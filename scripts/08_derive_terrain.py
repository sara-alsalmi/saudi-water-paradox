"""Derive terrain/hydrology variables (flow accumulation, TWI, HAND, slope) on a DEM extent
that contains the FULL upstream catchments of the AOI, plus upstream-mean event rainfall.

Stage A (coarse, terrain.coarse_res_m): condition DEM over dem.hydro_extent_wgs84, D8 flow,
  find every cell whose flow path enters the AOI (reverse BFS on the flow graph). Abort if
  that catchment touches the extent edge (extent too small -> expand config).
Stage B (fine, grid.resolution_m): recompute on the catchment bounding box (+pad) only:
  fill pits -> fill depressions -> resolve flats -> D8 -> accumulation, slope (Horn),
  TWI = ln(a / tan b) with a = acc_cells * res (specific catchment area, m), tan b >= tan_min,
  channels = acc_cells * res^2 >= hand_channel_area_m2, HAND w.r.t. those channels,
  rain_upstream_event_mm = flow-weighted mean of event rain over each cell's upstream area.

Sea: Copernicus DEM ocean = 0 m exactly; missing (ocean-only) tiles -> nodata. Both are set
to nodata so that land drains to the sea boundary. Copernicus DEM is a DSM (buildings incl.).

Output (interim): data/interim/terrain_30m.tif on a 60 m-snapped grid (nests with AOI grid),
                  data/interim/terrain_qa.json
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import rasterio
from affine import Affine
from rasterio.enums import Resampling
from rasterio.warp import reproject, transform_bounds
from scipy import ndimage
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import breadth_first_order

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, SNAP_M, aoi_bounds_utm, load_config  # noqa: E402

DIRMAP = (64, 128, 1, 2, 4, 8, 16, 32)          # N, NE, E, SE, S, SW, W, NW (pysheds default)
OFFS = {64: (-1, 0), 128: (-1, 1), 1: (0, 1), 2: (1, 1), 4: (1, 0), 8: (1, -1), 16: (0, -1), 32: (-1, -1)}
NODATA = -9999.0
INTERIM = ROOT / "data" / "interim"


def snap_bounds(b, pad=0.0):
    return (math.floor((b[0] - pad) / SNAP_M) * SNAP_M, math.floor((b[1] - pad) / SNAP_M) * SNAP_M,
            math.ceil((b[2] + pad) / SNAP_M) * SNAP_M, math.ceil((b[3] + pad) / SNAP_M) * SNAP_M)


def dem_on_grid(cfg, bounds, res):
    crs = cfg["grid"]["crs"]
    w = int(round((bounds[2] - bounds[0]) / res)); h = int(round((bounds[3] - bounds[1]) / res))
    tr = Affine(res, 0, bounds[0], 0, -res, bounds[3])
    out = np.full((h, w), np.nan, np.float32)
    for f in sorted((ROOT / cfg["paths"]["raw"] / "dem").glob("*.tif")):
        tmp = np.full((h, w), np.nan, np.float32)
        with rasterio.open(f) as src:
            reproject(rasterio.band(src, 1), tmp, dst_transform=tr, dst_crs=crs, dst_nodata=np.nan,
                      resampling=Resampling.bilinear if res <= 35 else Resampling.average)
        fill = np.isnan(out) & ~np.isnan(tmp)
        out[fill] = tmp[fill]
    sea = np.isnan(out) | (out == 0)
    return out, sea, tr


def condition(dem, sea, tr, crs, tag):
    """Write DEM (sea = nodata) and run pysheds conditioning + D8 + accumulation."""
    from pysheds.grid import Grid
    INTERIM.mkdir(parents=True, exist_ok=True)
    path = INTERIM / f"_dem_{tag}.tif"
    d = np.where(sea, NODATA, dem).astype(np.float64)
    with rasterio.open(path, "w", driver="GTiff", width=d.shape[1], height=d.shape[0], count=1, dtype="float64",
                       crs=crs, transform=tr, nodata=NODATA, compress="deflate") as dst:
        dst.write(d, 1)
    grid = Grid.from_raster(str(path))
    r = grid.read_raster(str(path))
    r = grid.fill_pits(r)
    r = grid.fill_depressions(r)
    r = grid.resolve_flats(r)
    fdir = grid.flowdir(r, dirmap=DIRMAP)
    acc = grid.accumulation(fdir, dirmap=DIRMAP)
    return grid, r, fdir, acc


def downstream_index(fdir):
    h, w = fdir.shape
    f = np.asarray(fdir)
    rows, cols = np.indices((h, w))
    tgt = np.full(h * w, -1, np.int64)
    for code, (dr, dc) in OFFS.items():
        m = f == code
        rr, cc = rows[m] + dr, cols[m] + dc
        ok = (rr >= 0) & (rr < h) & (cc >= 0) & (cc < w)
        src = (rows[m] * w + cols[m])[ok]
        tgt[src] = (rr * w + cc)[ok]
    return tgt


def upstream_of(mask_start, fdir):
    """Boolean mask of all cells whose D8 flow path reaches any cell in mask_start."""
    h, w = mask_start.shape
    n = h * w
    tgt = downstream_index(fdir)
    src = np.nonzero(tgt >= 0)[0]
    # reversed graph: downstream -> upstream, plus super-source node n -> start cells
    starts = np.flatnonzero(mask_start)
    rows = np.concatenate([tgt[src], np.full(starts.size, n)])
    cols = np.concatenate([src, starts])
    g = csr_matrix((np.ones(rows.size, np.int8), (rows, cols)), shape=(n + 1, n + 1))
    order = breadth_first_order(g, n, directed=True, return_predecessors=False)
    out = np.zeros(n + 1, bool)
    out[order] = True
    return out[:n].reshape(h, w)


def horn_slope_deg(z, res):
    zp = np.pad(z, 1, mode="edge")
    a, b, c = zp[:-2, :-2], zp[:-2, 1:-1], zp[:-2, 2:]
    d, f = zp[1:-1, :-2], zp[1:-1, 2:]
    g, h, i = zp[2:, :-2], zp[2:, 1:-1], zp[2:, 2:]
    dzdx = ((c + 2 * f + i) - (a + 2 * d + g)) / (8 * res)
    dzdy = ((g + 2 * h + i) - (a + 2 * b + c)) / (8 * res)
    return np.degrees(np.arctan(np.hypot(dzdx, dzdy)))


def main():
    cfg = load_config()
    tc = cfg["terrain"]
    crs = cfg["grid"]["crs"]
    res = cfg["grid"]["resolution_m"]
    aoi = aoi_bounds_utm(cfg)

    # ---------------- Stage A: catchment delineation (coarse)
    ext = snap_bounds(transform_bounds("EPSG:4326", crs, *cfg["dem"]["hydro_extent_wgs84"], densify_pts=101))
    cres = tc["coarse_res_m"]
    ext = (ext[0] + 2 * cres, ext[1] + 2 * cres, ext[2] - 2 * cres, ext[3] - 2 * cres)  # stay inside the tiles
    ext = snap_bounds(ext)
    dem_c, sea_c, tr_c = dem_on_grid(cfg, ext, cres)
    print(f"stage A grid {dem_c.shape} @ {cres} m")
    _, _, fdir_c, _ = condition(dem_c, sea_c, tr_c, crs, "coarse")
    inv = ~tr_c
    c0, r0 = inv * (aoi[0], aoi[3]); c1, r1 = inv * (aoi[2], aoi[1])
    aoi_mask = np.zeros(dem_c.shape, bool)
    aoi_mask[max(int(r0), 0):int(math.ceil(r1)), max(int(c0), 0):int(math.ceil(c1))] = True
    aoi_mask &= ~sea_c
    up = upstream_of(aoi_mask, fdir_c) & ~sea_c
    edge = np.zeros_like(up); edge[[0, -1], :] = True; edge[:, [0, -1]] = True
    touches = bool((up & edge).any())
    rr, cc = np.nonzero(up)
    xs, ys = tr_c * (np.array([cc.min(), cc.max() + 1]), np.array([rr.max() + 1, rr.min()]))
    catch_bounds = (float(xs[0]), float(ys[0]), float(xs[1]), float(ys[1]))
    catch_wgs = transform_bounds(crs, "EPSG:4326", *catch_bounds, densify_pts=101)
    print(f"catchment+AOI: {up.sum() * cres * cres / 1e6:.0f} km2, bbox WGS84 {np.round(catch_wgs, 3)}, touches edge: {touches}")
    if touches:
        raise SystemExit("Upstream catchment touches the DEM extent edge: enlarge dem.hydro_extent_wgs84 and re-run "
                         "02_download_dem.py (ask before >500 MB).")

    # ---------------- Stage B: fine terrain on catchment bbox
    fb = snap_bounds(catch_bounds, pad=tc["fine_pad_m"])
    dem, sea, tr = dem_on_grid(cfg, fb, res)
    print(f"stage B grid {dem.shape} @ {res} m")
    grid, cond, fdir, acc = condition(dem, sea, tr, crs, "fine")
    accn = np.asarray(acc, dtype=np.float64)
    land = ~sea
    slope = horn_slope_deg(np.where(sea, np.nan, dem).astype(np.float64), res)
    tanb = np.maximum(np.tan(np.radians(slope)), tc["twi_min_tan_slope"])
    with np.errstate(divide="ignore", invalid="ignore"):
        twi = np.log(accn * res / tanb)        # sea cells (acc=0) are set to NaN below
    channels = (accn * res * res) >= tc["hand_channel_area_m2"]
    # Sea as drainage outlet (approved): coastal outlet cells = land cells whose D8 receiver is sea,
    # or that have no receiver but touch the sea. HAND is measured to the first drainage cell
    # (channel OR coastal outlet) reached along the D8 path: HAND = z_cond(cell) - z_cond(drain).
    tgt = downstream_index(fdir)
    seaf = sea.ravel()
    to_sea = (tgt >= 0) & seaf[np.maximum(tgt, 0)]
    sea_adj = ndimage.binary_dilation(sea, structure=np.ones((3, 3), bool)).ravel()
    outlet = land.ravel() & (to_sea | ((tgt < 0) & sea_adj))
    del tgt
    drain = (channels & land).ravel() | outlet
    from pysheds.sview import Raster
    dr_r = Raster(drain.reshape(dem.shape), viewfinder=fdir.viewfinder)
    idx = np.asarray(grid.compute_hand(fdir, cond, dr_r, dirmap=DIRMAP, return_index=True)).ravel().astype(np.int64)
    zc = np.asarray(cond, dtype=np.float64).ravel()
    ok = idx >= 0
    hand = np.full(zc.size, np.nan)
    hand[ok] = zc[ok] - zc[idx[ok]]
    hand_ref_is_sea = np.zeros(zc.size, np.float32)
    hand_ref_is_sea[ok] = (outlet[idx[ok]] & ~(channels & land).ravel()[idx[ok]]).astype(np.float32)
    hand_ref_is_sea[~ok] = np.nan
    hand, hand_ref_is_sea = hand.reshape(dem.shape), hand_ref_is_sea.reshape(dem.shape)
    n_outlets = int(outlet.sum())
    del idx, zc, outlet, drain

    # upstream-mean event rain (flow-weighted accumulation / cell count)
    rain = np.full(dem.shape, np.nan, np.float32)
    with rasterio.open(INTERIM / "rain_0p1deg.tif") as src:
        assert src.descriptions[0] == "rain_event_24h_mm"
        reproject(rasterio.band(src, 1), rain, dst_transform=tr, dst_crs=crs, dst_nodata=np.nan,
                  resampling=Resampling.nearest)
    w_r = Raster(np.where(land, np.nan_to_num(rain, nan=0.0), 0.0).astype(np.float64), viewfinder=fdir.viewfinder)
    racc = np.asarray(grid.accumulation(fdir, weights=w_r, dirmap=DIRMAP), dtype=np.float64)
    rain_up = racc / accn
    rain_missing = int((land & np.isnan(rain)).sum())

    for a in (slope, twi, hand, rain_up, accn, hand_ref_is_sea):
        a[sea] = np.nan
    names = ["elevation_m", "slope_deg", "flowacc_cells", "twi", "hand_m", "rain_upstream_event_mm", "sea_mask",
             "hand_ref_is_sea"]
    arrs = [np.where(sea, np.nan, dem), slope, accn, twi, hand, rain_up, sea.astype(np.float32), hand_ref_is_sea]
    with rasterio.open(INTERIM / "terrain_30m.tif", "w", driver="GTiff", width=dem.shape[1], height=dem.shape[0],
                       count=len(arrs), dtype="float32", crs=crs, transform=tr, nodata=np.nan, compress="deflate",
                       tiled=True, predictor=3) as dst:
        for i, (n, a) in enumerate(zip(names, arrs), 1):
            dst.write(a.astype("float32"), i)
            dst.set_band_description(i, n)
    qa = {
        "coarse_res_m": cres, "fine_res_m": res, "fine_bounds_utm": fb,
        "catchment_plus_aoi_km2": round(float(up.sum() * cres * cres / 1e6), 1),
        "catchment_bbox_wgs84": [round(v, 4) for v in catch_wgs], "catchment_touches_extent_edge": touches,
        "acc_min_cells": float(np.nanmin(accn[land])), "hand_channel_area_m2": tc["hand_channel_area_m2"],
        "channel_cells": int((channels & land).sum()), "coastal_outlet_cells": n_outlets,
        "hand_nan_on_land": int(np.isnan(hand[land]).sum()),
        "hand_ref_is_sea_on_land": int(np.nansum(hand_ref_is_sea[land])),
        "hand_negative_on_land": int((hand[land] < -1e-6).sum()),
        "land_cells_without_event_rain": rain_missing,
        "unresolved_fdir_on_land": int(np.isin(np.asarray(fdir)[land], (-1, -2, 0)).sum()),
    }
    (INTERIM / "terrain_qa.json").write_text(json.dumps(qa, indent=2), encoding="utf-8")
    for p in INTERIM.glob("_dem_*.tif"):
        p.unlink()
    print(json.dumps(qa, indent=1))


if __name__ == "__main__":
    main()
