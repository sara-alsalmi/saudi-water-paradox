"""Shared helpers: config loading and reference-grid construction.

All grids are snapped to a 60 m origin so the 20 m check grid and the 30 m
reference grid nest exactly.
"""
from __future__ import annotations

import math
from datetime import datetime, timedelta, timezone
from pathlib import Path

import yaml
from affine import Affine
from rasterio.warp import transform_bounds

ROOT = Path(__file__).resolve().parents[1]
SNAP_M = 60.0


def load_config(path: Path | str = ROOT / "config.yaml") -> dict:
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def aoi_bounds_utm(cfg: dict) -> tuple[float, float, float, float]:
    """AOI bbox (WGS84) -> bounds in grid CRS, expanded outward to the 60 m snap."""
    b = transform_bounds("EPSG:4326", cfg["grid"]["crs"], *cfg["aoi"]["bbox_wgs84"], densify_pts=101)
    return (
        math.floor(b[0] / SNAP_M) * SNAP_M,
        math.floor(b[1] / SNAP_M) * SNAP_M,
        math.ceil(b[2] / SNAP_M) * SNAP_M,
        math.ceil(b[3] / SNAP_M) * SNAP_M,
    )


def clip_cog(href: str, bbox_wgs84, pad_m: float, out_path: Path) -> tuple[int, int]:
    """Save the window of a (remote) raster covering bbox_wgs84 (+pad, in metres) in its NATIVE
    CRS/resolution/values (no resampling). Returns (width, height); (0, 0) if no overlap."""
    import rasterio
    from rasterio.windows import Window, from_bounds

    with rasterio.open(href) as src:
        b = transform_bounds("EPSG:4326", src.crs, *bbox_wgs84, densify_pts=101)
        b = (b[0] - pad_m, b[1] - pad_m, b[2] + pad_m, b[3] + pad_m)
        w = from_bounds(*b, src.transform).round_offsets().round_lengths()
        try:
            w = w.intersection(Window(0, 0, src.width, src.height))
        except Exception:
            return 0, 0
        if w.width <= 0 or w.height <= 0:
            return 0, 0
        data = src.read(window=w)
        prof = src.profile.copy()
        prof.update(driver="GTiff", width=w.width, height=w.height, transform=src.window_transform(w),
                    compress="deflate", tiled=True, blockxsize=256, blockysize=256)
        for k in ("photometric",):
            prof.pop(k, None)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with rasterio.open(out_path, "w", **prof) as dst:
            dst.write(data)
            dst.update_tags(source_href=href)
        return int(w.width), int(w.height)


def iso(s: str) -> datetime:
    return datetime.fromisoformat(s.replace("Z", "+00:00"))


def events(cfg: dict, argv: list[str] | None = None) -> list[dict]:
    """Events from config (in config order) with derived time windows (all UTC datetimes):
      t0, event_end = t0 + rain.event_hours, ante_start = t0 - rain.antecedent_days,
      dl_end = gfm_post_date + 1 day (end of IMERG download), label_acq = label scene time.
    Optional argv: event_ids to keep (scripts accept them as command-line arguments)."""
    keep = set(argv or [])
    out = []
    for e in cfg["events"]:
        if keep and e["event_id"] not in keep:
            continue
        e = dict(e)
        e["t0"] = iso(e["t0_utc"])
        e["event_end"] = e["t0"] + timedelta(hours=cfg["rain"]["event_hours"])
        e["ante_start"] = e["t0"] - timedelta(days=cfg["rain"]["antecedent_days"])
        e["dl_end"] = datetime.fromisoformat(e["gfm_post_date"]).replace(tzinfo=timezone.utc) + timedelta(days=1)
        e["label_acq"] = iso(e["label_acq_utc"])
        e["year"] = e["t0"].year
        out.append(e)
    if keep and len(out) != len(keep):
        raise SystemExit(f"unknown event_id(s): {keep - {e['event_id'] for e in out}}")
    return out


def interim_path(name: str, event_id: str) -> Path:
    """Per-event interim file, e.g. data/interim/rain_0p1deg_JED_20221124.tif"""
    stem, dot, ext = name.rpartition(".")
    return ROOT / "data" / "interim" / f"{stem}_{event_id}.{ext}"


def make_grid(cfg: dict, res: float) -> tuple[Affine, int, int, tuple]:
    """Return (transform, width, height, bounds) for a north-up grid at `res` metres."""
    xmin, ymin, xmax, ymax = aoi_bounds_utm(cfg)
    width = int(round((xmax - xmin) / res))
    height = int(round((ymax - ymin) / res))
    return Affine(res, 0, xmin, 0, -res, ymax), width, height, (xmin, ymin, xmax, ymax)
