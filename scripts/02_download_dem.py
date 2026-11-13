"""Download Copernicus DEM GLO-30 1x1 deg tiles covering dem.hydro_extent_wgs84.

Copernicus DEM is a DSM (includes buildings/vegetation). Ocean-only tiles do not
exist in the bucket (HTTP 404) and are logged, not treated as errors.
Public bucket, no account needed. Abort if the total would exceed 500 MB.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, load_config  # noqa: E402
from manifest import already_have, record  # noqa: E402

CAP_BYTES = 500_000_000


def tile_names(bbox):
    w, s, e, n = bbox
    for lat in range(math.floor(s), math.ceil(n)):
        for lon in range(math.floor(w), math.ceil(e)):
            ns, ew = ("N" if lat >= 0 else "S"), ("E" if lon >= 0 else "W")
            yield f"{ns}{abs(lat):02d}_00_{ew}{abs(lon):03d}_00"


def main():
    cfg = load_config()
    d = cfg["dem"]
    out = ROOT / cfg["paths"]["raw"] / "dem"
    out.mkdir(parents=True, exist_ok=True)
    plan = []
    for t in tile_names(d["hydro_extent_wgs84"]):
        url = d["url_template"].format(tile=t)
        h = requests.head(url, timeout=60)
        if h.status_code == 404 or h.status_code == 403:
            print(f"{t}: no tile (ocean) — skipped")
            continue
        h.raise_for_status()
        plan.append((t, url, int(h.headers["Content-Length"])))
    total = sum(p[2] for p in plan)
    print(f"{len(plan)} tiles, {total/1e6:.0f} MB")
    if total > CAP_BYTES:
        sys.exit(f"ABORT: {total/1e6:.0f} MB exceeds 500 MB cap — ask before downloading.")
    for t, url, size in plan:
        path = out / Path(url).name
        if already_have(path):
            print(f"{t}: up to date"); continue
        with requests.get(url, stream=True, timeout=300) as r:
            r.raise_for_status()
            with open(path, "wb") as f:
                for chunk in r.iter_content(1 << 20):
                    f.write(chunk)
        assert path.stat().st_size == size, f"size mismatch for {t}"
        record("copernicus_dem_glo30", path, url, d["version"], "DSM (includes buildings); EGM2008 heights")
        print(f"{t}: {size/1e6:.1f} MB")


if __name__ == "__main__":
    main()
