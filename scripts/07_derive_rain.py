"""Derive rainfall variables per event from IMERG V07 half-hourly subsets (0.1 deg grid).

Usage: 07_derive_rain.py [EVENT_ID ...]   (default: all events)

Only files of the event's run are used (events[].imerg_short_name: GPM_3IMERGHH = Final
'3B-HHR.*', GPM_3IMERGHHL = Late '3B-HHR-L.*').

Outputs per event (interim, not deliverables):
  data/interim/rain_0p1deg_<event_id>.tif   EPSG:4326, bands:
      1 rain_event_24h_mm        sum over [t0, t0 + event_hours)            (trigger)
      2 rain_event_max30_mmhr    max half-hourly rate in the same window    (trigger)
      3 rain_antecedent_7d_mm    sum over [t0 - antecedent_days, t0)        (<= t0)
      4 qa_rain_after_event_mm   sum over [event end, label acquisition)    (QA ONLY, never a feature)
  data/interim/rain_qa_<event_id>.json      completeness + daily rain between event and label acquisition

Units: IMERG 'precipitation' is a rate in mm/hr per 30-min interval -> depth = rate * 0.5 mm.
Negative values (fill -9999.9) -> NaN. Any missing half-hour in a window aborts (no silent gaps).
"""
from __future__ import annotations

import json
import re
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import h5py
import numpy as np
import rasterio
from rasterio.transform import from_origin

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, events, interim_path, load_config  # noqa: E402

PAT = re.compile(r"3IMERG\.(\d{8})-S(\d{6})-E(\d{6})")
PREFIX = {"GPM_3IMERGHH": "3B-HHR.", "GPM_3IMERGHHL": "3B-HHR-L."}


def main():
    cfg = load_config()
    for ev in events(cfg, sys.argv[1:]):
        derive_event(cfg, ev)


def derive_event(cfg, ev):
    t0 = ev["t0"]
    eid = ev["event_id"]
    files = []
    for f in sorted((ROOT / cfg["paths"]["raw"] / "imerg").glob(PREFIX[ev["imerg_short_name"]] + "*.subset.nc4")):
        m = PAT.search(f.name)
        t = datetime.strptime(m.group(1) + m.group(2), "%Y%m%d%H%M%S").replace(tzinfo=timezone.utc)
        if ev["ante_start"] <= t < ev["dl_end"]:
            files.append(f)
    times, stack, lon, lat = [], [], None, None
    for f in files:
        m = PAT.search(f.name)
        times.append(datetime.strptime(m.group(1) + m.group(2), "%Y%m%d%H%M%S").replace(tzinfo=timezone.utc))
        with h5py.File(f) as h:
            p = h["Grid/precipitation"][...].astype(np.float64)
            dims = [d.decode() if isinstance(d, bytes) else d for d in h["Grid/precipitation"].attrs.get("DIMENSION_LIST_NAMES", [])]
            lo, la = h["Grid/lon"][:], h["Grid/lat"][:]
        p = np.squeeze(p, axis=0)                 # (lon, lat) for V07
        if p.shape == (lo.size, la.size):
            p = p.T                               # -> (lat, lon)
        assert p.shape == (la.size, lo.size), (f.name, p.shape)
        if lon is None:
            lon, lat = lo, la
        assert np.allclose(lo, lon) and np.allclose(la, lat)
        p[p < 0] = np.nan                         # fill value -9999.9
        stack.append(p)
    times = np.array(times)
    P = np.stack(stack)                           # (t, lat, lon) mm/hr
    order = np.argsort(times)
    times, P = times[order], P[order]

    def window(a, b, name):
        sel = (times >= a) & (times < b)
        expected = int((b - a) / timedelta(minutes=30))
        if sel.sum() != expected:
            raise SystemExit(f"[{eid}] {name}: {sel.sum()} of {expected} half-hours present — re-run 04_download_imerg.py")
        return sel

    # QA window ends at the half-hour boundary that CONTAINS the label acquisition (ceil)
    la = ev["label_acq"]
    qa_end = la.replace(minute=0, second=0, microsecond=0) + timedelta(minutes=30 if la.minute < 30 else 60)
    if la.minute in (0, 30) and la.second == 0:
        qa_end = la
    evw = window(t0, ev["event_end"], "event")
    an = window(ev["ante_start"], t0, "antecedent")
    qa = window(ev["event_end"], qa_end, "qa")

    bands = {
        "rain_event_24h_mm": np.nansum(P[evw] * 0.5, axis=0),
        "rain_event_max30_mmhr": np.nanmax(P[evw], axis=0),
        "rain_antecedent_7d_mm": np.nansum(P[an] * 0.5, axis=0),
        "qa_rain_after_event_mm": np.nansum(P[qa] * 0.5, axis=0),
    }
    nan_counts = {k: int(np.isnan(P[s]).sum()) for k, s in (("event", evw), ("antecedent", an), ("qa", qa))}

    # north-up raster
    flip = lat[0] < lat[-1]
    tr = from_origin(lon.min() - 0.05, lat.max() + 0.05, 0.1, 0.1)
    out = ROOT / "data" / "interim"
    out.mkdir(parents=True, exist_ok=True)
    with rasterio.open(interim_path("rain_0p1deg.tif", eid), "w", driver="GTiff", width=lon.size, height=lat.size, count=4,
                       dtype="float32", crs="EPSG:4326", transform=tr, nodata=np.nan) as dst:
        for i, (k, a) in enumerate(bands.items(), 1):
            dst.write((a[::-1] if flip else a).astype("float32"), i)
            dst.set_band_description(i, k)

    # QA: daily rain over cells intersecting the AOI, between event end and GFM acquisition
    w, s, e, n = cfg["aoi"]["bbox_wgs84"]
    li = (lon + 0.05 > w) & (lon - 0.05 < e)
    lj = (lat + 0.05 > s) & (lat - 0.05 < n)
    daily = []
    d = ev["event_end"]
    while d < qa_end:
        sel = (times >= d) & (times < min(d + timedelta(days=1), qa_end))
        depth = np.nansum(P[sel][:, lj][:, :, li] * 0.5, axis=0)
        daily.append({"from_utc": d.isoformat(), "to_utc": min(d + timedelta(days=1), qa_end).isoformat(),
                      "aoi_mean_mm": round(float(depth.mean()), 2), "aoi_max_mm": round(float(depth.max()), 2)})
        d += timedelta(days=1)
    ev_aoi = bands["rain_event_24h_mm"][np.ix_(lj, li)]
    qa_json = {
        "event_id": eid, "imerg_short_name": ev["imerg_short_name"], "imerg_run": ev["imerg_run"],
        "imerg_versions": sorted({f.name.split(".")[-3] for f in files}),
        "label_acq_utc": ev["label_acq_utc"], "qa_window_end_utc": qa_end.isoformat(),
        "n_files": len(files), "first": times[0].isoformat(), "last": times[-1].isoformat(),
        "nan_values_in_windows": nan_counts,
        "aoi_cells": int(lj.sum() * li.sum()),
        "event_24h_aoi_mean_mm": round(float(ev_aoi.mean()), 2), "event_24h_aoi_max_mm": round(float(ev_aoi.max()), 2),
        "event_max30_aoi_max_mmhr": round(float(bands["rain_event_max30_mmhr"][np.ix_(lj, li)].max()), 2),
        "antecedent_7d_aoi_max_mm": round(float(bands["rain_antecedent_7d_mm"][np.ix_(lj, li)].max()), 2),
        "rain_between_event_and_label": daily,
    }
    interim_path("rain_qa.json", eid).write_text(json.dumps(qa_json, indent=2), encoding="utf-8")
    print(json.dumps(qa_json, indent=1))


if __name__ == "__main__":
    main()
