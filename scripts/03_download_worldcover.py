"""Download ESA WorldCover 10 m tiles (3x3 deg) covering the AOI. Public, CC-BY 4.0."""
from __future__ import annotations

import math
import sys
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, load_config  # noqa: E402
from manifest import already_have, record  # noqa: E402


def tiles(bbox):
    w, s, e, n = bbox
    for lat in range(math.floor(s / 3) * 3, math.ceil(n / 3) * 3, 3):
        for lon in range(math.floor(w / 3) * 3, math.ceil(e / 3) * 3, 3):
            yield f"{'N' if lat >= 0 else 'S'}{abs(lat):02d}{'E' if lon >= 0 else 'W'}{abs(lon):03d}"


def main():
    cfg = load_config()
    lc = cfg["landcover"]
    out = ROOT / cfg["paths"]["raw"] / "worldcover"
    out.mkdir(parents=True, exist_ok=True)
    for t in tiles(cfg["aoi"]["bbox_wgs84"]):
        url = lc["url_template"].format(version=lc["version"], year=lc["year"], tile=t)
        path = out / Path(url).name
        if already_have(path):
            print(f"{t}: up to date"); continue
        with requests.get(url, stream=True, timeout=300) as r:
            r.raise_for_status()
            with open(path, "wb") as f:
                for chunk in r.iter_content(1 << 20):
                    f.write(chunk)
        record("esa_worldcover", path, url, f"{lc['year']} {lc['version']}", "CC-BY 4.0; map year 2021 <= t0")
        print(f"{t}: {path.stat().st_size/1e6:.1f} MB")


if __name__ == "__main__":
    main()
