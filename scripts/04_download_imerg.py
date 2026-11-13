"""Download GPM IMERG V07 Final half-hourly precipitation, subset to rain.bbox_wgs84 via OPeNDAP.

One small netCDF4 per half-hour in data/raw/imerg/. Requires NASA Earthdata credentials
(see nasa_opendap.login) and the "NASA GESDISC DATA ARCHIVE" application approved in the
user's Earthdata profile.
"""
from __future__ import annotations

import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import earthaccess
import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, load_config  # noqa: E402
from manifest import already_have, record  # noqa: E402
from nasa_opendap import get, login, opendap_url, subset_url, var_dims  # noqa: E402

VARS = ["precipitation", "precipitationQualityIndex"]


def idx_range(lo, hi, origin, step=0.1):
    """Index range of 0.1-deg cell centres (origin = first centre) fully covering [lo, hi]."""
    return int(np.floor((lo - origin) / step + 0.5 - 1e-9)), int(np.ceil((hi - origin) / step - 0.5 + 1e-9))


def main():
    cfg = load_config()
    rc = cfg["rain"]
    out = ROOT / cfg["paths"]["raw"] / "imerg"
    out.mkdir(parents=True, exist_ok=True)
    session = login()
    t_start = datetime.fromisoformat(rc["start_utc"].replace("Z", "+00:00"))
    t_end = datetime.fromisoformat(rc["end_utc"].replace("Z", "+00:00"))
    grans = earthaccess.search_data(short_name=rc["short_name"], version=rc["version"],
                                    temporal=(rc["start_utc"], rc["end_utc"]), count=-1)
    keep = []
    for g in grans:
        b = datetime.fromisoformat(g["umm"]["TemporalExtent"]["RangeDateTime"]["BeginningDateTime"].replace("Z", "+00:00"))
        if t_start <= b < t_end:
            keep.append(g)
    print(f"{len(keep)} half-hourly granules in [{rc['start_utc']}, {rc['end_utc']})")

    w, s, e, n = rc["bbox_wgs84"]
    i0, i1 = idx_range(w, e, -179.95)   # lon
    j0, j1 = idx_range(s, n, -89.95)    # lat
    dims = var_dims(session, opendap_url(keep[0]))
    pdims = dims["/Grid/precipitation"]
    print("precipitation dims:", pdims)
    sl = {"/Grid/time": "[0:0]", "/Grid/lon": f"[{i0}:{i1}]", "/Grid/lat": f"[{j0}:{j1}]"}
    ce = [f"/Grid/{v}" + "".join(sl[d] for d in pdims) for v in VARS]
    ce += [f"/Grid/lon[{i0}:{i1}]", f"/Grid/lat[{j0}:{j1}]", "/Grid/time"]

    def fetch(g):
        name = g["umm"]["GranuleUR"].split(":")[-1].replace(".HDF5", "") + ".subset.nc4"
        path = out / name
        if already_have(path):
            return name, "up to date"
        url = subset_url(opendap_url(g), ce)
        path.write_bytes(get(session, url).content)
        with h5py.File(path) as f:   # verify the subset really covers the bbox
            lon, lat = f["Grid/lon"][:], f["Grid/lat"][:]
            assert lon.min() <= w + 0.05 and lon.max() >= e - 0.05 and lat.min() <= s + 0.05 and lat.max() >= n - 0.05, name
        ver = g["umm"]["GranuleUR"].split(".")[-2] if "V07" in g["umm"]["GranuleUR"] else rc["version"]
        record("gpm_imerg_v07_final_hh", path, url, ver, "OPeNDAP spatial subset; precipitation in mm/hr")
        return name, "ok"

    with ThreadPoolExecutor(6) as ex:
        futs = [ex.submit(fetch, g) for g in keep]
        for k, f in enumerate(as_completed(futs), 1):
            name, st = f.result()
            if k % 48 == 0 or k == len(futs):
                print(f"{k}/{len(futs)} {name} {st}")


if __name__ == "__main__":
    main()
