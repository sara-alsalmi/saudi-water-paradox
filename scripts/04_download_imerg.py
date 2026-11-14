"""Download GPM IMERG V07 half-hourly precipitation per event, subset to rain.bbox_wgs84 via OPeNDAP.

Usage: 04_download_imerg.py [EVENT_ID ...]   (default: all events)

Window per event: [t0 - rain.antecedent_days, gfm_post_date + 1 day). Product per event:
events[].imerg_short_name — GPM_3IMERGHH (Final run) or GPM_3IMERGHHL (Late run, used only
where the Final run is not yet published). One small netCDF4 per half-hour in data/raw/imerg/
(file names carry the run: 3B-HHR = Final, 3B-HHR-L = Late). Requires NASA Earthdata credentials
(see nasa_opendap.login) and the "NASA GESDISC DATA ARCHIVE" application approved in the
user's Earthdata profile.
"""
from __future__ import annotations

import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import earthaccess
import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, events, iso, load_config  # noqa: E402
from manifest import already_have, record  # noqa: E402
from nasa_opendap import get, login, opendap_url, subset_url, var_dims  # noqa: E402

VARS = ["precipitation", "precipitationQualityIndex"]
SOURCE = {"GPM_3IMERGHH": "gpm_imerg_v07_final_hh", "GPM_3IMERGHHL": "gpm_imerg_v07_late_hh"}


def idx_range(lo, hi, origin, step=0.1):
    """Index range of 0.1-deg cell centres (origin = first centre) fully covering [lo, hi]."""
    return int(np.floor((lo - origin) / step + 0.5 - 1e-9)), int(np.ceil((hi - origin) / step - 0.5 + 1e-9))


def main():
    cfg = load_config()
    rc = cfg["rain"]
    out = ROOT / cfg["paths"]["raw"] / "imerg"
    out.mkdir(parents=True, exist_ok=True)
    session = login()
    w, s, e, n = rc["bbox_wgs84"]
    i0, i1 = idx_range(w, e, -179.95)   # lon
    j0, j1 = idx_range(s, n, -89.95)    # lat

    for ev in events(cfg, sys.argv[1:]):
        sn = ev["imerg_short_name"]
        t_start, t_end = ev["ante_start"], ev["dl_end"]
        grans = earthaccess.search_data(short_name=sn, version=rc["version"],
                                        temporal=(t_start.isoformat(), t_end.isoformat()), count=-1)
        keep = [g for g in grans
                if t_start <= iso(g["umm"]["TemporalExtent"]["RangeDateTime"]["BeginningDateTime"]) < t_end]
        expected = int((t_end - t_start).total_seconds() // 1800)
        print(f"[{ev['event_id']}] {sn}: {len(keep)} of {expected} half-hourly granules in [{t_start}, {t_end})")
        if len(keep) != expected:
            raise SystemExit("incomplete IMERG coverage for this window")

        dims = var_dims(session, opendap_url(keep[0]))
        pdims = dims["/Grid/precipitation"]
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
            ver = g["umm"]["GranuleUR"].split(".")[-2] if ".V07" in g["umm"]["GranuleUR"] else rc["version"]
            record(SOURCE[sn], path, url, ver, f"OPeNDAP spatial subset; precipitation in mm/hr; {ev['imerg_run']} run")
            return name, "ok"

        with ThreadPoolExecutor(6) as ex:
            futs = [ex.submit(fetch, g) for g in keep]
            for k, f in enumerate(as_completed(futs), 1):
                name, st = f.result()
                if k % 48 == 0 or k == len(futs):
                    print(f"{k}/{len(futs)} {name} {st}")


if __name__ == "__main__":
    main()
