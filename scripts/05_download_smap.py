"""Download SMAP L4 SPL4SMGP (3-hourly, 9 km EASE-Grid 2.0) per event — the latest granule whose
time window ENDS at or before the event's t0 — subset to soil_moisture.bbox_wgs84 via OPeNDAP.

Usage: 05_download_smap.py [EVENT_ID ...]   (default: all events)
"""
from __future__ import annotations

import sys
from datetime import timedelta
from pathlib import Path

import earthaccess
import h5py
import numpy as np
from pyproj import Transformer

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, events, iso, load_config  # noqa: E402
from manifest import already_have, record  # noqa: E402
from nasa_opendap import get, login, opendap_url, subset_url, var_dims  # noqa: E402


def main():
    cfg = load_config()
    session = login()
    for ev in events(cfg, sys.argv[1:]):
        print(f"[{ev['event_id']}]")
        download_event(cfg, session, ev["t0"])


def download_event(cfg, session, t0):
    sm = cfg["soil_moisture"]
    out = ROOT / cfg["paths"]["raw"] / "smap"
    out.mkdir(parents=True, exist_ok=True)

    grans = earthaccess.search_data(short_name=sm["short_name"], version=sm["version"],
                                    temporal=((t0 - timedelta(hours=12)).isoformat(), t0.isoformat()), count=-1)
    cands = []
    for g in grans:
        rng = g["umm"]["TemporalExtent"]["RangeDateTime"]
        end = iso(rng["EndingDateTime"])
        if end <= t0:
            cands.append((end, rng["BeginningDateTime"], g))
    end, begin, g = max(cands, key=lambda c: c[0])
    print(f"selected {g['umm']['GranuleUR']}  window {begin} -> {end.isoformat()} (<= t0 {t0.isoformat()})")

    base = opendap_url(g)
    dims = var_dims(session, base)
    tmp = out / "_xy.nc4"
    tmp.write_bytes(get(session, subset_url(base, ["/x", "/y"])).content)
    with h5py.File(tmp) as f:
        x, y = f["x"][:], f["y"][:]
    tmp.unlink()
    w, s, e, n = sm["bbox_wgs84"]
    tr = Transformer.from_crs("EPSG:4326", "EPSG:6933", always_xy=True)   # EASE-Grid 2.0 global
    xs, ys = tr.transform([w, e, w, e], [s, s, n, n])
    half = abs(x[1] - x[0]) / 2
    ix = np.where((x >= min(xs) - half) & (x <= max(xs) + half))[0]
    iy = np.where((y >= min(ys) - half) & (y <= max(ys) + half))[0]
    xr, yr = f"[{ix.min()}:{ix.max()}]", f"[{iy.min()}:{iy.max()}]"
    sl = {"/x": xr, "/y": yr}
    ce = []
    for v in sm["variables"]:
        key = f"/Geophysical_Data/{v}"
        ce.append(key + "".join(sl.get(d, "") for d in dims[key]))
    for v in ("/cell_lat", "/cell_lon"):
        ce.append(v + "".join(sl.get(d, "") for d in dims[v]))
    ce += [f"/x{xr}", f"/y{yr}", "/time"]
    print("dims:", {k: dims[k] for k in [f"/Geophysical_Data/{v}" for v in sm["variables"]]})

    path = out / g["umm"]["GranuleUR"].replace(".h5", ".subset.nc4")
    if already_have(path):
        print("up to date"); return
    url = subset_url(base, ce)
    path.write_bytes(get(session, url).content)
    with h5py.File(path) as f:
        print({k: f[k].shape for k in f.keys() if isinstance(f[k], h5py.Dataset)},
              {f"Geophysical_Data/{v}": f[f"Geophysical_Data/{v}"].shape for v in sm["variables"]})
    record("nasa_smap_l4_spl4smgp", path, url, g["umm"]["GranuleUR"].split("_")[-2],
           f"OPeNDAP subset; time-average window {begin} -> {end.isoformat()}; m3/m3; EASE-Grid 2.0 (EPSG:6933)")
    print("saved", path.name, path.stat().st_size, "bytes")


if __name__ == "__main__":
    main()
