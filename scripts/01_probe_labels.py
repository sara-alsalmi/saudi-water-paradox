"""Step 2 — Label verification (review only; outputs are NOT part of the dataset).

Clips GFM (pre/post event) and Sentinel-2 L2A (pre/post event) to the AOI on a
common 20 m grid in the reference CRS, and writes:
  docs/checks/label_check_overview.png   full-AOI comparison
  docs/checks/label_check_zooms.png      zooms on the densest GFM flood areas
  docs/checks/label_check_stats.md       numbers behind the maps

Nothing is written to data/. Remote COGs are read window-only (no full downloads).

GFM encoding (GFM PDD v4.0, Tables 19-20 and §4.4-4.5, plus values read from files):
  ensemble_flood_extent : 0 no flood, 1 flood, 255 no data
                          NOTE: excluded pixels are 0 here, so the exclusion mask
                          must be applied explicitly before calling a pixel "dry".
  exclusion_mask        : 0 not excluded, 1 excluded, 255 no data
  reference_water_mask  : 0 no water, 1 permanent, 2 seasonal, 255 no data
  advisory_flags        : 0 none, 1 low backscatter, 2 rough water, 3 both
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch, Rectangle
from pystac_client import Client
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT
from rasterio.warp import transform as warp_transform

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, load_config, make_grid  # noqa: E402

NODATA = 255
# class codes for the merged GFM map
C_NODATA, C_EXCL, C_DRY, C_FLOOD, C_REFW = 0, 1, 2, 3, 4
CLASS_NAMES = {C_NODATA: "not observed", C_EXCL: "excluded", C_DRY: "observed, no water",
               C_FLOOD: "flood", C_REFW: "reference water"}
CLASS_CMAP = ListedColormap(["#ffffff", "#9e9e9e", "#e8e2c8", "#d7191c", "#2c7bb6"])


def day_range(d: str) -> str:
    start = datetime.fromisoformat(d).replace(tzinfo=timezone.utc)
    return f"{start:%Y-%m-%dT%H:%M:%SZ}/{start + timedelta(days=1) - timedelta(seconds=1):%Y-%m-%dT%H:%M:%SZ}"


def read_to_grid(href: str, transform, width, height, crs, resampling, src_nodata=None, dtype=None):
    with rasterio.open(href) as src:
        nd = src.nodata if src_nodata is None else src_nodata
        with WarpedVRT(src, crs=crs, transform=transform, width=width, height=height,
                       resampling=resampling, src_nodata=nd, nodata=nd) as vrt:
            arr = vrt.read(1)
        meta = {"nodata": src.nodata, "dtype": src.dtypes[0], "crs": str(src.crs)}
    return (arr.astype(dtype) if dtype else arr), meta


# ----------------------------------------------------------------------------- GFM
def load_gfm_date(cfg, date, grid):
    transform, width, height, _ = grid
    crs = cfg["grid"]["crs"]
    cat = Client.open(cfg["labels"]["stac_api"])
    items = list(cat.search(collections=cfg["labels"]["collection"], bbox=cfg["aoi"]["bbox_wgs84"],
                            datetime=day_range(date), max_items=200).items())
    shape = (height, width)
    flood_any = np.zeros(shape, bool)
    dry_any = np.zeros(shape, bool)
    excl_any = np.zeros(shape, bool)
    refw = np.zeros(shape, bool)
    n_valid = np.zeros(shape, np.uint8)
    advisory = np.zeros(shape, np.uint8)
    item_log, uniques = [], defaultdict(set)
    for it in sorted(items, key=lambda i: i.id):
        layers = {}
        for k in cfg["labels"]["assets"]:
            arr, meta = read_to_grid(it.assets[k].href, transform, width, height, crs, Resampling.nearest,
                                     src_nodata=NODATA)
            layers[k] = arr
            uniques[k] |= set(np.unique(arr).tolist())
            uniques[k + "__nodata_tag"].add(meta["nodata"])
        f, e = layers["ensemble_flood_extent"], layers["exclusion_mask"]
        r, a = layers["reference_water_mask"], layers["advisory_flags"]
        observed = f != NODATA
        valid = observed & (e == 0)                      # observed AND not excluded
        if not observed.any():
            item_log.append((it.id, it.datetime.isoformat(), 0, 0, 0))
            continue
        flood_any |= valid & (f == 1)
        dry_any |= valid & (f == 0)
        excl_any |= observed & (e == 1)
        refw |= observed & np.isin(r, (1, 2))
        n_valid += valid.astype(np.uint8)
        advisory = np.maximum(advisory, np.where(observed & (a != NODATA), a, 0).astype(np.uint8))
        item_log.append((it.id, it.datetime.isoformat(), int(observed.sum()), int(valid.sum()),
                         int((valid & (f == 1)).sum())))
    conflict = flood_any & dry_any                       # overlapping scenes disagree
    cls = np.full(shape, C_NODATA, np.uint8)
    cls[excl_any] = C_EXCL
    cls[dry_any] = C_DRY
    cls[flood_any] = C_FLOOD                             # rule: flooded if ANY valid scene says flooded
    cls[refw & (cls != C_NODATA)] = C_REFW               # permanent/seasonal water overrides
    return {"cls": cls, "conflict": conflict, "advisory": advisory, "n_valid": n_valid,
            "items": item_log, "uniques": {k: sorted(v, key=str) for k, v in uniques.items()}}


# ----------------------------------------------------------------------- Sentinel-2
def load_s2_date(cfg, date, grid):
    transform, width, height, _ = grid
    crs = cfg["grid"]["crs"]
    v = cfg["verification"]
    cat = Client.open(v["s2_stac_api"])
    items = list(cat.search(collections=v["s2_collection"], bbox=cfg["aoi"]["bbox_wgs84"],
                            datetime=day_range(date), max_items=50).items())
    items = [i for i in items if i.properties.get("eo:cloud_cover", 100) <= v["s2_max_cloud_pct"]]
    bands = {"blue": None, "green": None, "red": None, "nir": None, "swir16": None}
    out = {b: np.full((height, width), np.nan, np.float32) for b in bands}
    scl_all = np.zeros((height, width), np.uint8)
    used = []
    for it in sorted(items, key=lambda i: i.id):
        scl, _ = read_to_grid(it.assets["scl"].href, transform, width, height, crs, Resampling.nearest)
        ok = ~np.isin(scl, (0, 1, 3, 8, 9, 10))           # drop nodata/saturated/shadow/cloud/cirrus
        for b in bands:
            a = it.assets[b]
            rb = (a.extra_fields.get("raster:bands") or [{}])[0]
            scale, offset, nd = rb.get("scale", 1e-4), rb.get("offset", 0.0), rb.get("nodata", 0)
            dn, _ = read_to_grid(a.href, transform, width, height, crs, Resampling.average, src_nodata=nd,
                                 dtype=np.float32)
            refl = np.where((dn != nd) & ok, dn * scale + offset, np.nan).astype(np.float32)
            fill = np.isnan(out[b]) & ~np.isnan(refl)
            out[b][fill] = refl[fill]
        scl_all = np.where(scl_all == 0, scl, scl_all)
        used.append((it.id, it.properties.get("eo:cloud_cover")))
    with np.errstate(invalid="ignore", divide="ignore"):
        out["ndwi"] = (out["green"] - out["nir"]) / (out["green"] + out["nir"])      # McFeeters 1996
        out["mndwi"] = (out["green"] - out["swir16"]) / (out["green"] + out["swir16"])  # Xu 2006
    out["scl"] = scl_all
    out["items"] = used
    return out


def rgb(s2):
    stack = np.dstack([s2["red"], s2["green"], s2["blue"]])
    lo, hi = np.nanpercentile(stack, 2), np.nanpercentile(stack, 98)
    img = np.clip((stack - lo) / (hi - lo), 0, 1) ** 0.8
    return np.nan_to_num(img, nan=1.0)


# --------------------------------------------------------------------------- main
def main():
    cfg = load_config()
    v = cfg["verification"]
    res = v["check_resolution_m"]
    grid = make_grid(cfg, res)
    transform, width, height, bounds = grid
    outdir = ROOT / cfg["paths"]["checks"]
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"check grid {width}x{height} @ {res} m, bounds {bounds}")

    gfm = {}
    for tag, d in (("pre", cfg["labels"]["pre_event_date"]), ("post", cfg["labels"]["post_event_date"])):
        print(f"GFM {tag} {d} ..."); gfm[tag] = load_gfm_date(cfg, d, grid)
    s2 = {}
    for tag, d in (("pre", v["s2_pre_date"]), ("post", v["s2_post_date"])):
        print(f"S2 {tag} {d} ..."); s2[tag] = load_s2_date(cfg, d, grid)
    dndwi = s2["post"]["ndwi"] - s2["pre"]["ndwi"]

    post = gfm["post"]["cls"]
    flood, dry = post == C_FLOOD, post == C_DRY
    px_km2 = res * res / 1e6

    # ---------------- statistics: does S2 see water where GFM says flood?
    def describe(mask):
        n = int(mask.sum())
        if n == 0:
            return {"n_px": 0}
        return {
            "n_px": n,
            "area_km2": round(n * px_km2, 3),
            "median_dNDWI": round(float(np.nanmedian(dndwi[mask])), 4),
            "pct_NDWI_post_gt0": round(100 * float(np.nanmean(s2["post"]["ndwi"][mask] > 0)), 2),
            "pct_MNDWI_post_gt0": round(100 * float(np.nanmean(s2["post"]["mndwi"][mask] > 0)), 2),
            "pct_dNDWI_gt_0.1": round(100 * float(np.nanmean(dndwi[mask] > 0.1)), 2),
        }

    land_obs = post != C_NODATA
    summary = {
        "grid": {"crs": cfg["grid"]["crs"], "res_m": res, "width": width, "height": height},
        "gfm_value_check": {t: gfm[t]["uniques"] for t in gfm},
        "gfm_items": {t: gfm[t]["items"] for t in gfm},
        "s2_items": {t: s2[t]["items"] for t in s2},
        "gfm_class_area_km2": {t: {CLASS_NAMES[c]: round(int((gfm[t]["cls"] == c).sum()) * px_km2, 2)
                                   for c in CLASS_NAMES} for t in gfm},
        "pre_flood_px": int((gfm["pre"]["cls"] == C_FLOOD).sum()),
        "post_overlap_conflict_px": int(gfm["post"]["conflict"].sum()),
        "post_excluded_share_of_observed_pct": round(100 * float((post == C_EXCL).sum() / max(land_obs.sum(), 1)), 1),
        "post_advisory_flag_px_on_flood": {int(k): int(((gfm["post"]["advisory"] == k) & flood).sum()) for k in range(4)},
        "s2_vs_gfm": {"GFM flood (29 Nov)": describe(flood), "GFM observed dry (29 Nov)": describe(dry)},
    }

    # ---------------- zoom selection: densest non-overlapping windows of GFM flood
    half = v["zoom_half_width_m"]
    win = int(round(2 * half / res))
    counts = []
    for r0 in range(0, height - win + 1, win // 2):
        for c0 in range(0, width - win + 1, win // 2):
            counts.append((int(flood[r0:r0 + win, c0:c0 + win].sum()), r0, c0))
    counts.sort(reverse=True)
    zooms = []
    for n, r0, c0 in counts:
        if n == 0 or len(zooms) >= v["n_zoom_clusters"]:
            break
        if all(abs(r0 - zr) >= win or abs(c0 - zc) >= win for _, zr, zc in zooms):
            zooms.append((n, r0, c0))
    zinfo = []
    for n, r0, c0 in zooms:
        x, y = transform * (c0 + win / 2, r0 + win / 2)
        lon, lat = warp_transform(cfg["grid"]["crs"], "EPSG:4326", [x], [y])
        sub = (slice(r0, r0 + win), slice(c0, c0 + win))
        zinfo.append({"center_lon": round(lon[0], 4), "center_lat": round(lat[0], 4), "flood_px": n,
                      "stats": describe(np.zeros_like(flood) | _sub_mask(flood, sub))})
    summary["zooms"] = zinfo

    # ---------------- figures
    ext = (bounds[0] / 1e3, bounds[2] / 1e3, bounds[1] / 1e3, bounds[3] / 1e3)
    fig, axs = plt.subplots(1, 5, figsize=(26, 10), constrained_layout=True)
    axs[0].imshow(rgb(s2["pre"]), extent=ext); axs[0].set_title(f"S2 true colour {v['s2_pre_date']} (pre)")
    axs[1].imshow(rgb(s2["post"]), extent=ext); axs[1].set_title(f"S2 true colour {v['s2_post_date']} (post)")
    im = axs[2].imshow(dndwi, extent=ext, cmap="RdBu", vmin=-0.3, vmax=0.3)
    axs[2].set_title("ΔNDWI (post − pre); blue = wetter")
    fig.colorbar(im, ax=axs[2], shrink=0.6)
    for ax, t, d in ((axs[3], "post", cfg["labels"]["post_event_date"]), (axs[4], "pre", cfg["labels"]["pre_event_date"])):
        ax.imshow(gfm[t]["cls"], extent=ext, cmap=CLASS_CMAP, vmin=0, vmax=4, interpolation="nearest")
        ax.set_title(f"GFM {d} ({'label date' if t == 'post' else 'baseline'})")
    for n, r0, c0 in zooms:
        x0, y1 = transform * (c0, r0)
        for ax in axs:
            ax.add_patch(Rectangle((x0 / 1e3, (y1 - 2 * half) / 1e3), 2 * half / 1e3, 2 * half / 1e3,
                                   fill=False, ec="k", lw=1.5))
    axs[4].legend(handles=[Patch(fc=CLASS_CMAP(i), ec="k", label=CLASS_NAMES[i]) for i in CLASS_NAMES],
                  loc="lower left", fontsize=9)
    for ax in axs:
        ax.set_xlabel("UTM 37N easting (km)")
    axs[0].set_ylabel("northing (km)")
    fig.suptitle("Jeddah 24 Nov 2022 — GFM label check vs Sentinel-2 (review only, not dataset)", fontsize=14)
    fig.savefig(outdir / "label_check_overview.png", dpi=110)
    plt.close(fig)

    if zooms:
        fig, axs = plt.subplots(len(zooms), 4, figsize=(20, 5 * len(zooms)), constrained_layout=True, squeeze=False)
        for i, (n, r0, c0) in enumerate(zooms):
            sl = (slice(r0, r0 + win), slice(c0, c0 + win))
            fl = flood[sl]
            z = zinfo[i]
            axs[i, 0].imshow(rgb({k: s2["pre"][k][sl] for k in ("red", "green", "blue")}))
            axs[i, 0].set_title(f"S2 pre {v['s2_pre_date']}  @ {z['center_lat']}N {z['center_lon']}E")
            axs[i, 1].imshow(rgb({k: s2["post"][k][sl] for k in ("red", "green", "blue")}))
            axs[i, 1].contour(fl, levels=[0.5], colors="red", linewidths=0.8)
            axs[i, 1].set_title(f"S2 post {v['s2_post_date']} + GFM flood outline (red)")
            im = axs[i, 2].imshow(dndwi[sl], cmap="RdBu", vmin=-0.3, vmax=0.3)
            axs[i, 2].contour(fl, levels=[0.5], colors="k", linewidths=0.6)
            axs[i, 2].set_title("ΔNDWI + GFM flood outline (black)")
            axs[i, 3].imshow(post[sl], cmap=CLASS_CMAP, vmin=0, vmax=4, interpolation="nearest")
            axs[i, 3].set_title(f"GFM {cfg['labels']['post_event_date']} classes ({n} flood px)")
            for ax in axs[i]:
                ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(im, ax=axs[:, 2], shrink=0.5)
        axs[0, 3].legend(handles=[Patch(fc=CLASS_CMAP(i), ec="k", label=CLASS_NAMES[i]) for i in CLASS_NAMES],
                         loc="lower left", fontsize=8)
        fig.savefig(outdir / "label_check_zooms.png", dpi=100)
        plt.close(fig)

    (outdir / "label_check_stats.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    write_md(summary, outdir / "label_check_stats.md", cfg)
    print(json.dumps({k: summary[k] for k in ("gfm_class_area_km2", "pre_flood_px", "post_overlap_conflict_px",
                                              "post_excluded_share_of_observed_pct", "s2_vs_gfm", "zooms",
                                              "gfm_value_check")}, indent=1, default=str))


def _sub_mask(mask, sub):
    m = np.zeros_like(mask)
    m[sub] = mask[sub]
    return m


def write_md(s, path, cfg):
    L = ["# GFM label check — statistics (review only)", "",
         f"Generated by `scripts/01_probe_labels.py` on {datetime.now(timezone.utc):%Y-%m-%d %H:%M} UTC. "
         f"Common grid: {s['grid']['crs']} @ {s['grid']['res_m']} m.", "",
         "Thresholds below (NDWI>0, MNDWI>0, ΔNDWI>0.1) are descriptive only, not labelling rules.", "",
         "## GFM class areas (km²)", "", "| date | " + " | ".join(CLASS_NAMES.values()) + " |",
         "|---|" + "---|" * len(CLASS_NAMES)]
    for t, d in s["gfm_class_area_km2"].items():
        L.append(f"| {t} | " + " | ".join(str(d[n]) for n in CLASS_NAMES.values()) + " |")
    L += ["", f"- Flood pixels on pre-event baseline: **{s['pre_flood_px']}**",
          f"- Overlapping post-event scenes disagreeing (flood vs dry): **{s['post_overlap_conflict_px']}** px",
          f"- Excluded share of observed area (post): **{s['post_excluded_share_of_observed_pct']}%**",
          f"- Advisory flags on flood pixels (0 none/1 low-bs/2 rough/3 both): {s['post_advisory_flag_px_on_flood']}",
          "", "## Sentinel-2 agreement", "", "| group | n px | km² | median ΔNDWI | % NDWI_post>0 | % MNDWI_post>0 | % ΔNDWI>0.1 |",
          "|---|---|---|---|---|---|---|"]
    for g, d in s["s2_vs_gfm"].items():
        L.append(f"| {g} | {d.get('n_px')} | {d.get('area_km2')} | {d.get('median_dNDWI')} | "
                 f"{d.get('pct_NDWI_post_gt0')} | {d.get('pct_MNDWI_post_gt0')} | {d.get('pct_dNDWI_gt_0.1')} |")
    L += ["", "## Zoom windows (densest GFM flood areas)", "", "| # | lat | lon | flood px | median ΔNDWI | % NDWI_post>0 | % MNDWI_post>0 |", "|---|---|---|---|---|---|---|"]
    for i, z in enumerate(s["zooms"], 1):
        st = z["stats"]
        L.append(f"| {i} | {z['center_lat']} | {z['center_lon']} | {z['flood_px']} | {st.get('median_dNDWI')} | "
                 f"{st.get('pct_NDWI_post_gt0')} | {st.get('pct_MNDWI_post_gt0')} |")
    L += ["", "## Source items", ""]
    for t, items in s["gfm_items"].items():
        L.append(f"**GFM {t}** (id, datetime, observed px, valid px, flood px):")
        L += [f"- `{i[0]}` {i[1]} obs={i[2]} valid={i[3]} flood={i[4]}" for i in items]
    for t, items in s["s2_items"].items():
        L.append(f"**Sentinel-2 {t}**: " + ", ".join(f"`{i}` (cloud {c}%)" for i, c in items))
    L += ["", "## Raw values found in GFM layers (incl. nodata tag)", "", "```", json.dumps(s["gfm_value_check"], indent=1, default=str), "```"]
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
