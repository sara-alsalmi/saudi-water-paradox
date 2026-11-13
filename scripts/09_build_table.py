"""Step 5 — Build data/processed/jeddah_flood_ml.{csv,parquet} and data_dictionary.csv.

Grids (all EPSG:32637, 60 m-snapped origins so they nest exactly):
  G30  reference grid over the AOI (output cells)
  G10  10 m sub-grid = G30 + clip pad, 3x3 subpixels per G30 cell (labels, S2)

Labels (GFM 2022-11-29, encoding: see scripts/01_probe_labels.py):
  each GFM scene is reprojected (nearest) to G10; a subpixel is VALID if observed (!=255) and
  not excluded (exclusion_mask==0) -- GFM codes excluded pixels as 0 in the flood layer, so
  the exclusion mask is applied explicitly. Overlapping scenes: flooded if ANY valid scene is
  flooded. Reference water (1 permanent, 2 seasonal) is never dry/negative.
  label=1  flood_fraction >= labels.flood_fraction_threshold
  label=0  all 9 subpixels valid & dry & not reference water, and >= negative_buffer_m from
           any flooded subpixel (distance computed on the padded grid, so edge-correct);
           stratified random sample (landcover x elev_bin, proportional allocation),
           neg_per_pos_ratio x n_pos, fixed seed.
Only cells whose centre lies inside the AOI WGS84 bbox are used.
"""
from __future__ import annotations

import csv
import json
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import h5py
import rasterio
from affine import Affine
from pyproj import Transformer
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT
from rasterio.windows import from_bounds
from scipy import ndimage
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, SNAP_M, load_config, make_grid  # noqa: E402
import schema  # noqa: E402

GFM_ND = 255
INTERIM = ROOT / "data" / "interim"


def warp(path_or_href, transform, width, height, crs, resampling, band=1, src_nodata=None, init=np.nan, dtype="float32"):
    with rasterio.open(path_or_href) as src:
        nd = src.nodata if src_nodata is None else src_nodata
        with WarpedVRT(src, crs=crs, transform=transform, width=width, height=height, resampling=resampling,
                       src_nodata=nd, nodata=nd) as vrt:
            a = vrt.read(band).astype(dtype)
    if nd is not None and not (isinstance(nd, float) and math.isnan(nd)):
        a[a == nd] = init
    return a


def blocks(a, off, H, W):
    """(H*3, W*3) inner region of a padded 10 m array -> (H, 3, W, 3) view."""
    return a[off:off + 3 * H, off:off + 3 * W].reshape(H, 3, W, 3)


def manifest_rows(source):
    with open(ROOT / "data" / "raw" / "manifest.csv", encoding="utf-8") as f:
        return [r for r in csv.DictReader(f) if r["source"] == source]


# ------------------------------------------------------------------------------ labels
def gfm_on_g10(cfg, T10, W10, H10):
    crs = cfg["grid"]["crs"]
    d = cfg["labels"]["post_event_date"]
    files = sorted((ROOT / cfg["paths"]["raw"] / "gfm" / d).glob("*.tif"))
    items = sorted({f.name.split("__")[1][:-4] for f in files})
    shape = (H10, W10)
    flood, dry, refw = np.zeros(shape, bool), np.zeros(shape, bool), np.zeros(shape, bool)
    obs_item = np.full(shape, -1, np.int8)
    item_times = []
    for k, it in enumerate(items):
        lay = {a: warp(ROOT / cfg["paths"]["raw"] / "gfm" / d / f"{a}__{it}.tif", T10, W10, H10, crs,
                       Resampling.nearest, src_nodata=GFM_ND, init=GFM_ND, dtype="uint8")
               for a in ("ensemble_flood_extent", "exclusion_mask", "reference_water_mask")}
        f, e, r = lay["ensemble_flood_extent"], lay["exclusion_mask"], lay["reference_water_mask"]
        observed = f != GFM_ND
        valid = observed & (e == 0)
        flood |= valid & (f == 1)
        dry |= valid & (f == 0)
        refw |= np.isin(r, (1, 2))
        obs_item[valid & (obs_item < 0)] = k
        m = re.search(r"(\d{8}T\d{6})", it)
        item_times.append(datetime.strptime(m.group(1), "%Y%m%dT%H%M%S").replace(tzinfo=timezone.utc))
    dry &= ~flood & ~refw
    flood &= ~refw          # GFM already removes reference water from flood; enforce
    return flood, dry, refw, obs_item, item_times


# --------------------------------------------------------------------------------- S2
def s2_cell_indices(cfg, T10i, W10i, H10i, H, W):
    crs = cfg["grid"]["crs"]
    meta = {}
    for r in manifest_rows("sentinel2_l2a"):
        m = re.search(r"scale=([-\d.e]+) offset=([-\d.e]+) nodata=([-\d.e]+)", r["notes"])
        meta[Path(r["local_path"]).name] = tuple(float(v) for v in m.groups()) if m else (1e-4, 0.0, 0.0)
    out = {}
    for tag, d in (("pre", cfg["verification"]["s2_pre_date"]), ("post", cfg["verification"]["s2_post_date"])):
        folder = ROOT / cfg["paths"]["raw"] / "s2" / d
        items = sorted({f.name.split("__")[1][:-4] for f in folder.glob("*.tif")})
        bands = {b: np.full((H10i, W10i), np.nan, np.float32) for b in ("green", "nir", "swir16")}
        for it in items:
            scl = warp(folder / f"scl__{it}.tif", T10i, W10i, H10i, crs, Resampling.nearest, src_nodata=0, init=0,
                       dtype="uint8")
            ok = ~np.isin(scl, (0, 1, 3, 8, 9, 10))        # nodata, saturated, shadow, cloud med/high, cirrus
            for b in bands:
                name = f"{b}__{it}.tif"
                sc, of, nd = meta[name]
                dn = warp(folder / name, T10i, W10i, H10i, crs, Resampling.nearest, src_nodata=nd, init=np.nan)
                refl = np.where(ok, dn * sc + of, np.nan).astype(np.float32)
                fill = np.isnan(bands[b]) & ~np.isnan(refl)
                bands[b][fill] = refl[fill]
            del scl, ok
        with np.errstate(invalid="ignore", divide="ignore"):
            ndwi = (bands["green"] - bands["nir"]) / (bands["green"] + bands["nir"])
            mndwi = (bands["green"] - bands["swir16"]) / (bands["green"] + bands["swir16"])
        both = ~np.isnan(ndwi) & ~np.isnan(mndwi)
        nd4, md4, v4 = (a.reshape(H, 3, W, 3) for a in (np.where(both, ndwi, np.nan), np.where(both, mndwi, np.nan), both))
        with np.errstate(invalid="ignore"):
            out[tag] = {"ndwi": np.nanmean(nd4, axis=(1, 3)), "mndwi": np.nanmean(md4, axis=(1, 3)),
                        "n_valid": v4.sum(axis=(1, 3))}
        del bands, ndwi, mndwi
    return out


# ------------------------------------------------------------------------------- main
def main():
    cfg = load_config()
    lab, ver = cfg["labels"], cfg["verification"]
    crs = cfg["grid"]["crs"]
    T30, W, H, b30 = make_grid(cfg, cfg["grid"]["resolution_m"])
    pad = math.ceil(lab["clip_pad_m"] / SNAP_M) * SNAP_M
    off = int(pad // 10)
    T10 = Affine(10, 0, b30[0] - pad, 0, -10, b30[3] + pad)
    W10, H10 = 3 * W + 2 * off, 3 * H + 2 * off
    print(f"G30 {W}x{H}, G10 {W10}x{H10}")

    # cell centres + AOI membership
    cols, rows = np.meshgrid(np.arange(W), np.arange(H))
    X = b30[0] + (cols + 0.5) * 30.0
    Y = b30[3] - (rows + 0.5) * 30.0
    to_wgs = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
    LON, LAT = to_wgs.transform(X, Y)
    w_, s_, e_, n_ = cfg["aoi"]["bbox_wgs84"]
    in_aoi = (LON >= w_) & (LON <= e_) & (LAT >= s_) & (LAT <= n_)

    # ---- labels
    flood10, dry10, refw10, obs_item10, item_times = gfm_on_g10(cfg, T10, W10, H10)
    flood_frac = blocks(flood10, off, H, W).mean(axis=(1, 3))
    all_dry = blocks(dry10, off, H, W).all(axis=(1, 3))
    dist_flood10 = ndimage.distance_transform_edt(~flood10).astype(np.float32) * 10.0 if flood10.any() \
        else np.full(flood10.shape, np.inf, np.float32)
    min_dist_flood = blocks(dist_flood10, off, H, W).min(axis=(1, 3))
    del dist_flood10
    dist_ref = (ndimage.distance_transform_edt(~refw10).astype(np.float32) * 10.0)[off + 1:off + 3 * H:3, off + 1:off + 3 * W:3]
    obs_item_c = blocks(obs_item10, off, H, W).max(axis=(1, 3))   # any valid subpixel's scene (-1 = none)
    del dry10, refw10

    pos = in_aoi & (flood_frac >= lab["flood_fraction_threshold"])
    neg_cand = in_aoi & all_dry & (min_dist_flood >= lab["negative_buffer_m"]) & ~pos

    # ---- terrain (interim grid nests with G30)
    with rasterio.open(INTERIM / "terrain_30m.tif") as t:
        win = from_bounds(*b30, t.transform)
        assert abs(win.col_off - round(win.col_off)) < 1e-6 and abs(win.row_off - round(win.row_off)) < 1e-6
        win = win.round_offsets().round_lengths()
        terr = {t.descriptions[i]: t.read(i + 1, window=win) for i in range(t.count)}
        # coast distance on a padded window so the nearest sea is found even near AOI edges
        cp = 3000
        wpad = from_bounds(b30[0] - cp, b30[1] - cp, b30[2] + cp, b30[3] + cp, t.transform).round_offsets().round_lengths()
        sea = t.read(t.descriptions.index("sea_mask") + 1, window=wpad, boundless=True, fill_value=1) == 1
    lab_sea, n_sea = ndimage.label(sea, structure=np.ones((3, 3)))
    sizes = ndimage.sum(sea, lab_sea, range(1, n_sea + 1))
    big = np.isin(lab_sea, 1 + np.nonzero(sizes >= 1_000_000 / 900)[0])
    dco = ndimage.distance_transform_edt(~big) * 30.0
    pc = cp // 30
    dist_coast = dco[pc:pc + H, pc:pc + W]

    # ---- land cover (mode of 10 m -> 30 m)
    wc_files = sorted((ROOT / cfg["paths"]["raw"] / "worldcover").glob("*.tif"))
    lc = np.zeros((H, W), np.float32)
    for f in wc_files:
        a = warp(f, T30, W, H, crs, Resampling.mode, src_nodata=0, init=0, dtype="float32")
        lc = np.where((lc == 0) & (a > 0), a, lc)
    lc[lc == 0] = np.nan

    # ---- rain
    rain = {}
    with rasterio.open(INTERIM / "rain_0p1deg.tif") as r:
        descs = list(r.descriptions)
    for name in ("rain_event_24h_mm", "rain_event_max30_mmhr", "rain_antecedent_7d_mm"):
        rain[name] = warp(INTERIM / "rain_0p1deg.tif", T30, W, H, crs, Resampling.nearest, band=descs.index(name) + 1)

    # ---- SMAP (nearest EASE-Grid 2.0 cell; coastal no-data -> nearest valid cell, flagged)
    smf = sorted((ROOT / cfg["paths"]["raw"] / "smap").glob("*.subset.nc4"))
    assert len(smf) == 1, smf
    with h5py.File(smf[0]) as h:
        ex, ey = h["x"][:], h["y"][:]
        sm = {}
        for v in cfg["soil_moisture"]["variables"]:
            ds = h[f"Geophysical_Data/{v}"]
            a = np.squeeze(ds[...]).astype(np.float64)
            fv = ds.attrs.get("_FillValue", -9999.0)
            a[(a == np.squeeze(fv)) | (a < 0)] = np.nan
            sm[v] = a
    to_ease = Transformer.from_crs(crs, "EPSG:6933", always_xy=True)
    EX, EY = to_ease.transform(X, Y)
    ix = np.clip(np.rint((EX - ex[0]) / (ex[1] - ex[0])).astype(int), 0, ex.size - 1)
    iy = np.clip(np.rint((EY - ey[0]) / (ey[1] - ey[0])).astype(int), 0, ey.size - 1)
    gx, gy = np.meshgrid(ex, ey)
    valid_e = ~np.isnan(sm["sm_surface"]) & ~np.isnan(sm["sm_rootzone"])
    tree = cKDTree(np.c_[gx[valid_e], gy[valid_e]])
    vy, vx = np.nonzero(valid_e)
    need = ~valid_e[iy, ix]
    _, nn = tree.query(np.c_[EX[need], EY[need]])
    iy2, ix2 = iy.copy(), ix.copy()
    iy2[need], ix2[need] = vy[nn], vx[nn]
    sm_surface = sm["sm_surface"][iy2, ix2]
    sm_rootzone = sm["sm_rootzone"][iy2, ix2]
    sm_filled = need.astype(np.int8)

    # ---- S2 label-quality flag (AOI-only 10 m grid)
    T10i = Affine(10, 0, b30[0], 0, -10, b30[3])
    s2 = s2_cell_indices(cfg, T10i, 3 * W, 3 * H, H, W)
    ok = (s2["pre"]["n_valid"] >= 5) & (s2["post"]["n_valid"] >= 5)
    s2_flag = np.where(ok, (((s2["post"]["ndwi"] - s2["pre"]["ndwi"]) > ver["s2_dndwi_min"]) |
                            (s2["post"]["mndwi"] > ver["s2_mndwi_post_min"])).astype(np.float32), np.nan)

    # ---- elevation bins & strata
    elev = terr["elevation_m"]
    elev_bin = pd.cut(pd.Series(elev.ravel()), bins=schema.ELEV_BINS, labels=schema.ELEV_LABELS, right=False)
    elev_bin = elev_bin.astype(object).where(elev_bin.notna(), "missing").to_numpy().reshape(H, W)
    lc_name = np.vectorize(lambda v: schema.WORLDCOVER.get(int(v), "missing") if not np.isnan(v) else "missing")(lc)
    stratum = np.char.add(np.char.add(lc_name.astype(str), "|"), elev_bin.astype(str))

    # ---- negative sampling (fixed seed); see labels.negative_strategy / extra_variants
    n_pos = int(pos.sum())
    n_neg_target = int(round(lab["neg_per_pos_ratio"] * n_pos))
    cand_idx = np.flatnonzero(neg_cand)
    elev_flat = elev.ravel()

    def largest_remainder(quota, total):
        alloc = np.floor(quota).astype(int)
        alloc[np.argsort(-(quota - alloc), kind="stable")[:total - alloc.sum()]] += 1
        return alloc

    def sample_stratified():
        """Proportional stratified random sample over landcover x elev_bin."""
        rng = np.random.default_rng(lab["seed"])
        cand_strata = stratum.ravel()[cand_idx]
        uniq, counts = np.unique(cand_strata, return_counts=True)
        info = {"candidates_per_stratum": {s: int(c) for s, c in zip(uniq, counts)}, "shortfall": {}}
        if n_neg_target >= cand_idx.size:
            return cand_idx, stratum.ravel().astype(object), info
        alloc = largest_remainder(counts / counts.sum() * n_neg_target, n_neg_target)
        chosen = np.concatenate([rng.choice(cand_idx[cand_strata == s], size=k, replace=False)
                                 for s, k in zip(uniq, alloc) if k > 0])
        return chosen, stratum.ravel().astype(object), info

    def sample_elevation_matched():
        """Negatives follow the positives' elevation distribution: quantile bins of POSITIVE
        elevations; per bin, quota = ratio x positives in bin; random draw without replacement.
        Candidates outside [min, max] of positive elevation are never drawn."""
        rng = np.random.default_rng(lab["seed"])
        nb = int(lab["elev_match_quantile_bins"])
        pidx = np.flatnonzero(pos)
        pe = elev_flat[pidx]
        q = np.unique(np.quantile(pe[~np.isnan(pe)], np.linspace(0, 1, nb + 1)))
        nbin = len(q) - 1
        ce = elev_flat[cand_idx]
        in_rng = (ce >= q[0]) & (ce <= q[-1])
        bc = np.clip(np.searchsorted(q, ce, side="right") - 1, 0, nbin - 1)
        bp = np.clip(np.searchsorted(q, pe, side="right") - 1, 0, nbin - 1)
        alloc = largest_remainder(np.bincount(bp, minlength=nbin) * lab["neg_per_pos_ratio"], n_neg_target)
        names = [f"elevq{k:02d}[{q[k]:.2f},{q[k + 1]:.2f})" for k in range(nbin)]
        chosen, shortfall, cand_per = [], {}, {}
        lab_arr = np.full(H * W, "", dtype=object)
        lab_arr[pidx] = np.array(names, dtype=object)[bp]
        for k in range(nbin):
            pool = cand_idx[in_rng & (bc == k)]
            cand_per[names[k]] = int(pool.size)
            m = min(int(alloc[k]), pool.size)
            if m < alloc[k]:
                shortfall[names[k]] = int(alloc[k] - m)
            if m:
                pick = rng.choice(pool, size=m, replace=False)
                chosen.append(pick)
                lab_arr[pick] = names[k]
        return np.concatenate(chosen), lab_arr, {"candidates_per_stratum": cand_per, "shortfall": shortfall,
                                                  "positive_elev_quantile_edges_m": [round(float(v), 3) for v in q]}

    samplers = {"stratified": sample_stratified, "elevation_matched": sample_elevation_matched}
    suffix = {"stratified": "_stratified", "elevation_matched": "_elev_matched"}
    default = lab["negative_strategy"]
    variants = [(default, "jeddah_flood_ml")] + [(v, "jeddah_flood_ml" + suffix[v])
                                                for v in lab.get("extra_variants", []) if v != default]

    gfm_ver = sorted({r["version"] for r in manifest_rows("copernicus_gfm")})
    label_source = f"Copernicus GFM ensemble_flood_extent ({', '.join(gfm_ver)}), Sentinel-1 {lab['post_event_date']}"
    obs_t = np.array([t.strftime("%Y-%m-%dT%H:%M:%SZ") for t in item_times] + [""], dtype=object)
    out = ROOT / cfg["paths"]["processed"]
    out.mkdir(parents=True, exist_ok=True)
    cand_below10 = float(np.mean(elev_flat[cand_idx] < 10))
    summary = {
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "grid": {"crs": crs, "res_m": 30, "width": W, "height": H, "bounds": b30},
        "aoi_cells": int(in_aoi.sum()),
        "n_pos": n_pos, "n_neg_candidates": int(cand_idx.size), "n_neg_target": n_neg_target,
        "candidates_below_10m_share": round(cand_below10, 4),
        "pos_cells_flood_frac_between_0_and_thr": int((in_aoi & (flood_frac > 0) & (flood_frac < lab["flood_fraction_threshold"])).sum()),
        "default_strategy": default, "variants": {},
    }
    for strategy, fname in variants:
        chosen, strat_flat, info = samplers[strategy]()
        df = assemble(cfg, pos, chosen, X, Y, LON, LAT, elev, terr, lc, lc_name, elev_bin, dist_coast, dist_ref, rain,
                      sm_surface, sm_rootzone, sm_filled, flood_frac, s2_flag, strat_flat, obs_t, obs_item_c, label_source)
        df.to_csv(out / f"{fname}.csv", index=False, float_format="%.6g")
        df.to_parquet(out / f"{fname}.parquet", index=False)
        neg = df[df.label == 0]
        summary["variants"][fname] = {
            "strategy": strategy, "rows": len(df), "n_pos": int((df.label == 1).sum()), "n_neg": len(neg),
            "neg_below_10m_share": round(float((neg.elevation_m < 10).mean()), 4),
            "pos_below_10m_share": round(float((df[df.label == 1].elevation_m < 10).mean()), 4),
            "label_obs_times": sorted(set(df["label_obs_utc"])),
            "elevation_below_0_rows": int((df["elevation_m"] < 0).sum()),
            **info,
        }

    schema.write_dictionary(ROOT / "data_dictionary.csv")

    (out / "build_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({k: ({n: {kk: vv for kk, vv in v.items() if kk != "candidates_per_stratum"}
                           for n, v in summary["variants"].items()} if k == "variants" else v)
                      for k, v in summary.items()}, indent=1))


def assemble(cfg, pos, chosen, X, Y, LON, LAT, elev, terr, lc, lc_name, elev_bin, dist_coast, dist_ref, rain,
             sm_surface, sm_rootzone, sm_filled, flood_frac, s2_flag, strat_flat, obs_t, obs_item_c, label_source):
    sel = np.zeros(pos.size, bool)
    sel[np.flatnonzero(pos)] = True
    sel[chosen] = True
    idx = np.flatnonzero(sel)
    flat = lambda a: np.asarray(a).ravel()[idx]  # noqa: E731
    df = pd.DataFrame({
        "cell_id": idx.astype(np.int64),
        "event_id": cfg["project"]["event_id"],
        "x": flat(X), "y": flat(Y), "lon": np.round(flat(LON), 6), "lat": np.round(flat(LAT), 6),
        "block_id": [f"{int(a // 1000)}_{int(b // 1000)}" for a, b in zip(flat(X), flat(Y))],
        "t0_utc": cfg["time"]["t0_utc"],
        "elevation_m": flat(elev),
        "slope_deg": flat(terr["slope_deg"]),
        "flowacc_km2": flat(terr["flowacc_cells"]) * 900 / 1e6,
        "twi": flat(terr["twi"]),
        "hand_m": flat(terr["hand_m"]),
        "hand_ref_is_sea": flat(terr["hand_ref_is_sea"]),
        "landcover": pd.array(np.where(np.isnan(flat(lc)), pd.NA, flat(lc)).astype(object), dtype="Int16"),
        "landcover_name": flat(lc_name),
        "elev_bin": flat(elev_bin),
        "dist_to_coast_m": flat(dist_coast),
        "dist_to_ref_water_m": flat(dist_ref),
        "rain_antecedent_7d_mm": flat(rain["rain_antecedent_7d_mm"]),
        "sm_surface_m3m3": flat(sm_surface),
        "sm_rootzone_m3m3": flat(sm_rootzone),
        "sm_filled_from_neighbor": flat(sm_filled),
        "rain_event_24h_mm": flat(rain["rain_event_24h_mm"]),
        "rain_event_max30_mmhr": flat(rain["rain_event_max30_mmhr"]),
        "rain_upstream_event_mm": flat(terr["rain_upstream_event_mm"]),
        "label": flat(pos).astype(np.int8),
        "label_source": label_source,
        "label_obs_utc": obs_t[flat(obs_item_c)],
        "flood_fraction": np.round(flat(flood_frac), 4),
        "s2_water_confirmed": flat(s2_flag),
        "sampling_stratum": flat(strat_flat),
    })
    assert list(df.columns) == schema.NAMES, set(df.columns) ^ set(schema.NAMES)
    df = df.sort_values(["label", "cell_id"], ascending=[False, True]).reset_index(drop=True)
    for c in df.select_dtypes("float").columns:
        df[c] = df[c].astype("float64").round(6)
    return df


if __name__ == "__main__":
    main()
