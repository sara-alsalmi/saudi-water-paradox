"""Download label sources as native-resolution clips (no resampling, values unchanged):

  data/raw/gfm/<date>/<asset>__<item_id>.tif   GFM pre (baseline) and post (label) dates
  data/raw/s2/<date>/<band>__<item_id>.tif     Sentinel-2 L2A (verification only -> s2_water_confirmed)

Clips cover the AOI plus labels.clip_pad_m so the negative buffer is correct at AOI edges.
Public sources; no account needed.
"""
from __future__ import annotations

import sys
from pathlib import Path

from pystac_client import Client

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, clip_cog, load_config  # noqa: E402
from manifest import already_have, record  # noqa: E402

S2_BANDS = ["green", "nir", "swir16", "scl"]


def day(d):
    return f"{d}T00:00:00Z/{d}T23:59:59Z"


def main():
    cfg = load_config()
    lab, ver = cfg["labels"], cfg["verification"]
    bbox, pad = cfg["aoi"]["bbox_wgs84"], lab["clip_pad_m"]
    raw = ROOT / cfg["paths"]["raw"]

    gfm = Client.open(lab["stac_api"])
    for d in (lab["pre_event_date"], lab["post_event_date"]):
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
                       f"window clip {w}x{h}, native Equi7Grid 20 m; acquired {it.datetime.isoformat()}")
            print(f"GFM {d} {it.id} ({proc_ver})")

    s2 = Client.open(ver["s2_stac_api"])
    for d in (ver["s2_pre_date"], ver["s2_post_date"]):
        for it in s2.search(collections=ver["s2_collection"], bbox=bbox, datetime=day(d), max_items=50).items():
            cc = it.properties.get("eo:cloud_cover", 100)
            if cc > ver["s2_max_cloud_pct"]:
                print(f"S2 {it.id}: cloud {cc}% > max, skipped"); continue
            for b in S2_BANDS:
                href = it.assets[b].href
                path = raw / "s2" / d / f"{b}__{it.id}.tif"
                if already_have(path):
                    continue
                w, h = clip_cog(href, bbox, pad, path)
                if w == 0:
                    continue
                rb = (it.assets[b].extra_fields.get("raster:bands") or [{}])[0]
                record("sentinel2_l2a", path, href, it.properties.get("s2:processing_baseline", ""),
                       f"window clip {w}x{h}; cloud {cc:.3f}%; scale={rb.get('scale')} offset={rb.get('offset')} "
                       f"nodata={rb.get('nodata')}; VERIFICATION ONLY")
            print(f"S2 {d} {it.id}")


if __name__ == "__main__":
    main()
