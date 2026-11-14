"""Step 5 — Build data/processed/jeddah_flood_ml.{csv,parquet} (+ variants) for ALL events in
config.yaml, and data_dictionary.csv. One row per (cell_id, event_id).

Grids (all EPSG:32637, 60 m-snapped origins so they nest exactly):
  G30  reference grid over the AOI (output cells; same cell_id in every event)
  G10  10 m sub-grid = G30 + clip pad, 3x3 subpixels per G30 cell (labels, S2)

Static layers (computed once): terrain/hydrology, land cover (WorldCover 2021 for every event ->
landcover_year_gap), distance to coast. Per event: labels, distance to reference water, rain,
SMAP, Sentinel-2 flag, negative sampling.

Labels per event (GFM, encoding: see scripts/01_probe_labels.py):
  each GFM scene is reprojected (nearest) to G10; a subpixel is VALID if observed (!=255) and
  not excluded (exclusion_mask==0) -- GFM codes excluded pixels as 0 in the flood layer, so
  the exclusion mask is applied explicitly. Overlapping scenes: flooded if ANY valid scene is
  flooded. Reference water (1 permanent, 2 seasonal) is never flood/dry.
  "New flood vs baseline" (approved): subpixels flooded on the baseline scene (gfm_pre_date) are
  BASELINE WATER -> never flood, never dry, and they count in the negative buffer.
  label=1  new-flood fraction >= labels.flood_fraction_threshold
  label=0  all 9 subpixels valid & dry & not reference water & not baseline water, and
           >= negative_buffer_m from any new-flood or baseline-water subpixel (distance on the
           padded grid, so edge-correct); sampled per event (neg_per_pos_ratio x n_pos of that
           event, labels.seed for every event).
Only cells whose centre lies inside the AOI WGS84 bbox are used.
"""
from __future__ import annotations

import csv
import json
import math
import re
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import rasterio
from affine import Affine
from pyproj import Transformer
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT
from rasterio.windows import from_bounds
from scipy import ndimage
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, SNAP_M, events, interim_path, load_config, make_grid  # noqa: E402
import schema  # noqa: E402

GFM_ND = 255
INTERIM = ROOT / "data" / "interim"
SMAP_PAT = re.compile(r"SMAP_L4_SM_gph_(\d{8}T\d{6})_")


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
def gfm_layers(cfg, d, T10, W10, H10):
    """Merged GFM layers of all scenes acquired on date d (folder data/raw/gfm/<d>)."""
    crs = cfg["grid"]["crs"]
    folder = ROOT / cfg["paths"]["raw"] / "gfm" / d
    items = sorted({f.name.split("__")[1][:-4] for f in folder.glob("*.tif")})
    if not items:
        raise SystemExit(f"no GFM files for {d} — run 06_download_labels.py")
    shape = (H10, W10)
    flood, dry, refw, valid_any = (np.zeros(shape, bool) for _ in range(4))
    obs_item = np.full(shape, -1, np.int8)
    item_times = []
    for k, it in enumerate(items):
        lay = {a: warp(folder / f"{a}__{it}.tif", T10, W10, H10, crs,
                       Resampling.nearest, src_nodata=GFM_ND, init=GFM_ND, dtype="uint8")
               for a in ("ensemble_flood_extent", "exclusion_mask", "reference_water_mask")}
        f, e, r = lay["ensemble_flood_extent"], lay["exclusion_mask"], lay["reference_water_mask"]
        observed = f != GFM_ND
        valid = observed & (e == 0)
        flood |= valid & (f == 1)
        dry |= valid & (f == 0)
        refw |= np.isin(r, (1, 2))
        valid_any |= valid
        obs_item[valid & (obs_item < 0)] = k
        m = re.search(r"(\d{8}T\d{6})", it)
        item_times.append(datetime.strptime(m.group(1), "%Y%m%dT%H%M%S").replace(tzinfo=timezone.utc))
    return flood, dry, refw, valid_any, obs_item, item_times


def gfm_on_g10(cfg, ev, T10, W10, H10):
    flood, dry, refw, _, obs_item, item_times = gfm_layers(cfg, ev["gfm_post_date"], T10, W10, H10)
    pre_flood, _, _, pre_valid, _, _ = gfm_layers(cfg, ev["gfm_pre_date"], T10, W10, H10)
    base = pre_flood & ~refw                    # baseline water: neither class, inside negative buffer
    stats = {"label_scene_flood_subpx": int((flood & ~refw).sum()),
             "baseline_water_subpx": int(base.sum()),
             "label_flood_on_baseline_water_subpx": int((flood & ~refw & base).sum()),
             "label_flood_where_baseline_not_valid_subpx": int((flood & ~refw & ~pre_valid).sum())}
    dry &= ~flood & ~refw & ~base
    flood &= ~refw & ~base                      # GFM already removes reference water from flood; enforce
    stats["new_flood_subpx"] = int(flood.sum())
    return flood, dry, refw, base, obs_item, item_times, stats


# --------------------------------------------------------------------------------- S2
def s2_cell_indices(cfg, ev, T10i, W10i, H10i, H, W):
    """Per tag (pre/post): dates in priority order, items sorted by id; for every band the first
    clear (SCL) observation per pixel wins."""
    crs = cfg["grid"]["crs"]
    meta = {}
    for r in manifest_rows("sentinel2_l2a"):
        m = re.search(r"scale=([-\d.e]+) offset=([-\d.e]+) nodata=([-\d.e]+)", r["notes"])
        meta[Path(r["local_path"]).name] = tuple(float(v) for v in m.groups()) if m else (1e-4, 0.0, 0.0)
    out = {}
    for tag, dates in (("pre", ev["s2_pre_dates"]), ("post", ev["s2_post_dates"])):
        bands = {b: np.full((H10i, W10i), np.nan, np.float32) for b in ("green", "nir", "swir16")}
        for d in dates:
            folder = ROOT / cfg["paths"]["raw"] / "s2" / d
            items = sorted({f.name.split("__")[1][:-4] for f in folder.glob("*.tif")})
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


# ------------------------------------------------------------------------------- SMAP
def smap_file(cfg, t0):
    """Latest SPL4SMGP granule whose 3-h window [centre - 1.5 h, centre + 1.5 h) ENDS <= t0."""
    best = None
    for f in sorted((ROOT / cfg["paths"]["raw"] / "smap").glob("*.subset.nc4")):
        c = datetime.strptime(SMAP_PAT.search(f.name).group(1), "%Y%m%dT%H%M%S").replace(tzinfo=timezone.utc)
        end = c + timedelta(minutes=90)
        if end <= t0 and t0 - end < timedelta(hours=3) and (best is None or end > best[0]):
            best = (end, f)
    if best is None:
        raise SystemExit(f"no SMAP granule ending <= {t0} within 3 h — run 05_download_smap.py")
    return best[1]


def smap_on_cells(cfg, t0, X, Y, crs):
    """Nearest EASE-Grid 2.0 cell; coastal no-data -> nearest valid cell (flagged)."""
    with h5py.File(smap_file(cfg, t0)) as h:
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
    return sm["sm_surface"][iy2, ix2], sm["sm_rootzone"][iy2, ix2], need.astype(np.int8)


# --------------------------------------------------------------------------- sampling
def largest_remainder(quota, total):
    alloc = np.floor(quota).astype(int)
    alloc[np.argsort(-(quota - alloc), kind="stable")[:total - alloc.sum()]] += 1
    return alloc


def sample_stratified(lab, cand_idx, n_neg_target, stratum, **_):
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


def sample_elevation_matched(lab, cand_idx, n_neg_target, pos, elev_flat, **_):
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
    lab_arr = np.full(pos.size, "", dtype=object)
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


SAMPLERS = {"stratified": sample_stratified, "elevation_matched": sample_elevation_matched}
SUFFIX = {"stratified": "_stratified", "elevation_matched": "_elev_matched"}


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
    evs = events(cfg)
    print(f"G30 {W}x{H}, G10 {W10}x{H10}; events {[e['event_id'] for e in evs]}")

    # cell centres + AOI membership
    cols, rows = np.meshgrid(np.arange(W), np.arange(H))
    X = b30[0] + (cols + 0.5) * 30.0
    Y = b30[3] - (rows + 0.5) * 30.0
    to_wgs = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
    LON, LAT = to_wgs.transform(X, Y)
    w_, s_, e_, n_ = cfg["aoi"]["bbox_wgs84"]
    in_aoi = (LON >= w_) & (LON <= e_) & (LAT >= s_) & (LAT <= n_)

    # ================= static layers
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

    # ---- land cover (mode of 10 m -> 30 m); WorldCover 2021 for every event
    wc_files = sorted((ROOT / cfg["paths"]["raw"] / "worldcover").glob("*.tif"))
    lc = np.zeros((H, W), np.float32)
    for f in wc_files:
        a = warp(f, T30, W, H, crs, Resampling.mode, src_nodata=0, init=0, dtype="float32")
        lc = np.where((lc == 0) & (a > 0), a, lc)
    lc[lc == 0] = np.nan

    # ---- elevation bins & strata
    elev = terr["elevation_m"]
    elev_flat = elev.ravel()
    elev_bin = pd.cut(pd.Series(elev_flat), bins=schema.ELEV_BINS, labels=schema.ELEV_LABELS, right=False)
    elev_bin = elev_bin.astype(object).where(elev_bin.notna(), "missing").to_numpy().reshape(H, W)
    lc_name = np.vectorize(lambda v: schema.WORLDCOVER.get(int(v), "missing") if not np.isnan(v) else "missing")(lc)
    stratum = np.char.add(np.char.add(lc_name.astype(str), "|"), elev_bin.astype(str))
    T10i = Affine(10, 0, b30[0], 0, -10, b30[3])

    default = lab["negative_strategy"]
    variants = [(default, "jeddah_flood_ml")] + [(v, "jeddah_flood_ml" + SUFFIX[v])
                                                for v in lab.get("extra_variants", []) if v != default]
    gfm_rows = manifest_rows("copernicus_gfm")
    frames = {fname: [] for _, fname in variants}
    summary = {
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "dataset_version": cfg["project"].get("version", 1),
        "grid": {"crs": crs, "res_m": 30, "width": W, "height": H, "bounds": b30},
        "aoi_cells": int(in_aoi.sum()),
        "default_strategy": default, "events": {}, "variants": {},
    }

    # ================= per event
    for ev in evs:
        eid = ev["event_id"]
        print(f"[{eid}] labels ...")
        flood10, dry10, refw10, base10, obs_item10, item_times, gstats = gfm_on_g10(cfg, ev, T10, W10, H10)
        flood_frac = blocks(flood10, off, H, W).mean(axis=(1, 3))
        all_dry = blocks(dry10, off, H, W).all(axis=(1, 3))
        block10 = flood10 | base10              # negatives keep the buffer from new flood AND baseline water
        dist_flood10 = ndimage.distance_transform_edt(~block10).astype(np.float32) * 10.0 if block10.any() \
            else np.full(block10.shape, np.inf, np.float32)
        min_dist_flood = blocks(dist_flood10, off, H, W).min(axis=(1, 3))
        del dist_flood10, block10
        dist_ref = (ndimage.distance_transform_edt(~refw10).astype(np.float32) * 10.0)[off + 1:off + 3 * H:3, off + 1:off + 3 * W:3]
        obs_item_c = blocks(obs_item10, off, H, W).max(axis=(1, 3))   # any valid subpixel's scene (-1 = none)
        base_c = blocks(base10, off, H, W).any(axis=(1, 3))
        del dry10, refw10, base10, flood10

        pos = in_aoi & (flood_frac >= lab["flood_fraction_threshold"])
        neg_cand = in_aoi & all_dry & (min_dist_flood >= lab["negative_buffer_m"]) & ~pos

        # ---- rain (per event)
        rpath = interim_path("rain_0p1deg.tif", eid)
        with rasterio.open(rpath) as r:
            descs = list(r.descriptions)
        rain = {name: warp(rpath, T30, W, H, crs, Resampling.nearest, band=descs.index(name) + 1)
                for name in ("rain_event_24h_mm", "rain_event_max30_mmhr", "rain_antecedent_7d_mm")}
        rain["rain_upstream_event_mm"] = terr[f"rain_upstream_event_mm__{eid}"]

        # ---- SMAP (per event)
        sm_surface, sm_rootzone, sm_filled = smap_on_cells(cfg, ev["t0"], X, Y, crs)

        # ---- S2 label-quality flag (AOI-only 10 m grid)
        print(f"[{eid}] Sentinel-2 ...")
        s2 = s2_cell_indices(cfg, ev, T10i, 3 * W, 3 * H, H, W)
        ok = (s2["pre"]["n_valid"] >= 5) & (s2["post"]["n_valid"] >= 5)
        s2_flag = np.where(ok, (((s2["post"]["ndwi"] - s2["pre"]["ndwi"]) > ver["s2_dndwi_min"]) |
                                (s2["post"]["mndwi"] > ver["s2_mndwi_post_min"])).astype(np.float32), np.nan)

        # ---- per-event metadata
        gfm_ver = sorted({r["version"] for r in gfm_rows
                          if f"/gfm/{ev['gfm_post_date']}/" in r["local_path"] or f"/gfm/{ev['gfm_pre_date']}/" in r["local_path"]})
        label_source = f"Copernicus GFM ensemble_flood_extent ({', '.join(gfm_ver)}), Sentinel-1 {ev['gfm_post_date']}"
        meta = {"event_id": eid, "t0_utc": ev["t0_utc"], "label_source": label_source,
                "gfm_version": ", ".join(v.replace("GFM ", "") for v in gfm_ver),
                "imerg_run": ev["imerg_run"],
                "landcover_year_gap": int(ev["year"] - int(cfg["landcover"]["year"]))}
        obs_t = np.array([t.strftime("%Y-%m-%dT%H:%M:%SZ") for t in item_times] + [""], dtype=object)

        n_pos = int(pos.sum())
        n_neg_target = int(round(lab["neg_per_pos_ratio"] * n_pos))
        cand_idx = np.flatnonzero(neg_cand)
        summary["events"][eid] = {
            "t0_utc": ev["t0_utc"], "gfm_pre_date": ev["gfm_pre_date"], "gfm_post_date": ev["gfm_post_date"],
            "gfm_version": meta["gfm_version"], "imerg_run": ev["imerg_run"],
            "landcover_year_gap": meta["landcover_year_gap"],
            "s2_pre_dates": ev["s2_pre_dates"], "s2_post_dates": ev["s2_post_dates"],
            **gstats,
            "aoi_cells_touching_baseline_water": int((in_aoi & base_c).sum()),
            "n_pos": n_pos, "n_neg_candidates": int(cand_idx.size), "n_neg_target": n_neg_target,
            "candidates_below_10m_share": round(float(np.mean(elev_flat[cand_idx] < 10)), 4),
            "pos_cells_flood_frac_between_0_and_thr": int((in_aoi & (flood_frac > 0) & (flood_frac < lab["flood_fraction_threshold"])).sum()),
            "s2_flag_valid_share_aoi": round(float(ok[in_aoi].mean()), 4),
        }
        for strategy, fname in variants:
            chosen, strat_flat, info = SAMPLERS[strategy](lab=lab, cand_idx=cand_idx, n_neg_target=n_neg_target,
                                                          stratum=stratum, pos=pos.ravel(), elev_flat=elev_flat)
            df = assemble(meta, pos, chosen, X, Y, LON, LAT, elev, terr, lc, lc_name, elev_bin, dist_coast, dist_ref,
                          rain, sm_surface, sm_rootzone, sm_filled, flood_frac, s2_flag, strat_flat, obs_t, obs_item_c)
            frames[fname].append(df)
            neg = df[df.label == 0]
            v = summary["variants"].setdefault(fname, {"strategy": strategy, "events": {}})
            v["events"][eid] = {
                "rows": len(df), "n_pos": int((df.label == 1).sum()), "n_neg": len(neg),
                "neg_below_10m_share": round(float((neg.elevation_m < 10).mean()), 4),
                "pos_below_10m_share": round(float((df[df.label == 1].elevation_m < 10).mean()), 4),
                "label_obs_times": sorted(set(df["label_obs_utc"])),
                "elevation_below_0_rows": int((df["elevation_m"] < 0).sum()),
                **info,
            }
        del s2, flood_frac, all_dry, min_dist_flood

    out = ROOT / cfg["paths"]["processed"]
    out.mkdir(parents=True, exist_ok=True)
    for strategy, fname in variants:
        df = pd.concat(frames[fname], ignore_index=True)
        df.to_csv(out / f"{fname}.csv", index=False, float_format="%.6g")
        df.to_parquet(out / f"{fname}.parquet", index=False)
        v = summary["variants"][fname]
        v.update({"rows": len(df), "n_pos": int((df.label == 1).sum()), "n_neg": int((df.label == 0).sum())})

    schema.write_dictionary(ROOT / "data_dictionary.csv")
    (out / "build_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({"events": summary["events"],
                      "variants": {n: {k: (vv if k != "events" else
                                           {e: {kk: x for kk, x in d.items() if kk != "candidates_per_stratum"}
                                            for e, d in vv.items()})
                                       for k, vv in v.items()} for n, v in summary["variants"].items()}},
                     indent=1))


def assemble(meta, pos, chosen, X, Y, LON, LAT, elev, terr, lc, lc_name, elev_bin, dist_coast, dist_ref, rain,
             sm_surface, sm_rootzone, sm_filled, flood_frac, s2_flag, strat_flat, obs_t, obs_item_c):
    sel = np.zeros(pos.size, bool)
    sel[np.flatnonzero(pos)] = True
    sel[chosen] = True
    idx = np.flatnonzero(sel)
    flat = lambda a: np.asarray(a).ravel()[idx]  # noqa: E731
    df = pd.DataFrame({
        "cell_id": idx.astype(np.int64),
        "event_id": meta["event_id"],
        "x": flat(X), "y": flat(Y), "lon": np.round(flat(LON), 6), "lat": np.round(flat(LAT), 6),
        "block_id": [f"{int(a // 1000)}_{int(b // 1000)}" for a, b in zip(flat(X), flat(Y))],
        "t0_utc": meta["t0_utc"],
        "elevation_m": flat(elev),
        "slope_deg": flat(terr["slope_deg"]),
        "flowacc_km2": flat(terr["flowacc_cells"]) * 900 / 1e6,
        "twi": flat(terr["twi"]),
        "hand_m": flat(terr["hand_m"]),
        "hand_ref_is_sea": flat(terr["hand_ref_is_sea"]),
        "landcover": pd.array(np.where(np.isnan(flat(lc)), pd.NA, flat(lc)).astype(object), dtype="Int16"),
        "landcover_name": flat(lc_name),
        "landcover_year_gap": np.full(idx.size, meta["landcover_year_gap"], np.int16),
        "elev_bin": flat(elev_bin),
        "dist_to_coast_m": flat(dist_coast),
        "dist_to_ref_water_m": flat(dist_ref),
        "rain_antecedent_7d_mm": flat(rain["rain_antecedent_7d_mm"]),
        "sm_surface_m3m3": flat(sm_surface),
        "sm_rootzone_m3m3": flat(sm_rootzone),
        "sm_filled_from_neighbor": flat(sm_filled),
        "rain_event_24h_mm": flat(rain["rain_event_24h_mm"]),
        "rain_event_max30_mmhr": flat(rain["rain_event_max30_mmhr"]),
        "rain_upstream_event_mm": flat(rain["rain_upstream_event_mm"]),
        "imerg_run": meta["imerg_run"],
        "label": flat(pos).astype(np.int8),
        "label_source": meta["label_source"],
        "gfm_version": meta["gfm_version"],
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
