"""Download label sources per event as native-resolution clips (no resampling, values unchanged):

  data/raw/gfm/<date>/<asset>__<item_id>.tif   GFM baseline (gfm_pre_date) and label (gfm_post_date) scenes
  data/raw/s2/<date>/<band>__<item_id>.tif     Sentinel-2 L2A (verification only -> s2_water_confirmed)

Usage: 06_download_labels.py [EVENT_ID ...] [--s2-estimate]
  --s2-estimate   only print the uncompressed size of the S2 windows that would be downloaded

GFM clips cover the AOI plus labels.clip_pad_m so the negative buffer is correct at AOI edges.
S2 (minimised): only scenes with tile cloud <= verification.s2_max_cloud_pct, only
verification.s2_bands, windows clipped to the AOI + verification.s2_clip_pad_m.
Public sources; no account needed.
"""
from __future__ import annotations

import sys
from pathlib import Path

import rasterio
from pystac_client import Client
from rasterio.warp import transform_bounds
from rasterio.windows import Window, from_bounds

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, clip_cog, events, load_config  # noqa: E402
from manifest import already_have, record  # noqa: E402


def day(d):
    return f"{d}T00:00:00Z/{d}T23:59:59Z"


def s2_items(cat, ver, bbox, d):
    out = []
    for it in sorted(cat.search(collections=ver["s2_collection"], bbox=bbox, datetime=day(d), max_items=50).items(),
                     key=lambda i: i.id):
        cc = it.properties.get("eo:cloud_cover", 100)
        if cc > ver["s2_max_cloud_pct"]:
            print(f"S2 {it.id}: cloud {cc}% > max, skipped")
            continue
        out.append(it)
    return out


def window_bytes(href, bbox, pad):
    with rasterio.open(href) as src:
        b = transform_bounds("EPSG:4326", src.crs, *bbox, densify_pts=101)
        w = from_bounds(b[0] - pad, b[1] - pad, b[2] + pad, b[3] + pad, src.transform).round_offsets().round_lengths()
        try:
            w = w.intersection(Window(0, 0, src.width, src.height))
        except Exception:
            return 0
        return int(w.width * w.height * src.count * {"uint8": 1, "uint16": 2}.get(src.dtypes[0], 4))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    estimate = "--s2-estimate" in sys.argv
    cfg = load_config()
    lab, ver = cfg["labels"], cfg["verification"]
    bbox, pad = cfg["aoi"]["bbox_wgs84"], lab["clip_pad_m"]
    raw = ROOT / cfg["paths"]["raw"]
    gfm = Client.open(lab["stac_api"])
    s2 = Client.open(ver["s2_stac_api"])
    total = 0

    for ev in events(cfg, args):
        eid = ev["event_id"]
        if not estimate:
            for d in (ev["gfm_pre_date"], ev["gfm_post_date"]):
                for it in gfm.search(collections=lab["collection"], bbox=bbox, datetime=day(d), max_items=200).items():
                    proc_ver = it.properties.get("processing:version", "")
                    for a in lab["assets"]:
                        href = it.assets[a].href
                        path = raw / "gfm" / d / f"{a}__{it.id}.tif"
                        if already_have(path):
                            continue
                        w, h = clip_cog(href, bbox, pad, path)
                        if w == 0:
                            print(f"GFM {it.id} {a}: no overlap"); continue
                        record("copernicus_gfm", path, href, proc_ver,
                               f"window clip {w}x{h}, native Equi7Grid 20 m; acquired {it.datetime.isoformat()}; {eid}")
                    print(f"[{eid}] GFM {d} {it.id} ({proc_ver})")

        for d in ev["s2_pre_dates"] + ev["s2_post_dates"]:
            for it in s2_items(s2, ver, bbox, d):
                for b in ver["s2_bands"]:
                    href = it.assets[b].href
                    path = raw / "s2" / d / f"{b}__{it.id}.tif"
                    if estimate:
                        if not already_have(path):
                            total += window_bytes(href, bbox, ver["s2_clip_pad_m"])
                        continue
                    if already_have(path):
                        continue
                    w, h = clip_cog(href, bbox, ver["s2_clip_pad_m"], path)
                    if w == 0:
                        continue
                    cc = it.properties.get("eo:cloud_cover", 100)
                    rb = (it.assets[b].extra_fields.get("raster:bands") or [{}])[0]
                    record("sentinel2_l2a", path, href, it.properties.get("s2:processing_baseline", ""),
                           f"window clip {w}x{h}; cloud {cc:.3f}%; scale={rb.get('scale')} offset={rb.get('offset')} "
                           f"nodata={rb.get('nodata')}; VERIFICATION ONLY; {eid}")
                print(f"[{eid}] S2 {d} {it.id}")
    if estimate:
        print(f"S2 new windows, uncompressed: {total / 1e6:.0f} MB (files are deflate-compressed, typically smaller)")


if __name__ == "__main__":
    main()
