"""Step 7 — Final QA of every table variant in data/processed -> docs/qa_report.md (fails loudly).

Checks run on the whole file AND per event (every event in config.yaml must be present).
"""
from __future__ import annotations

import csv
import hashlib
import json
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, events, interim_path, load_config  # noqa: E402
import schema  # noqa: E402


def md_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    out = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        out.append("| " + " | ".join("" if (isinstance(v, float) and np.isnan(v)) else str(v) for v in r) + " |")
    return "\n".join(out)


def smap_windows():
    """{file name: window end} from the manifest notes."""
    out = {}
    for r in csv.DictReader(open(ROOT / "data/raw/manifest.csv", encoding="utf-8")):
        if r["source"] == "nasa_smap_l4_spl4smgp":
            out[Path(r["local_path"]).name] = datetime.fromisoformat(r["notes"].split("-> ")[1].split(";")[0])
    return out


def file_checks(cfg, df, dfc, evs):
    checks = []

    def check(name, ok, detail=""):
        checks.append(("ALL", name, "PASS" if ok else "FAIL", detail))

    check("Columns match schema (order)", list(df.columns) == schema.NAMES)
    check("Key (cell_id, event_id) unique", not df.duplicated(["cell_id", "event_id"]).any(),
          f"{len(df)} rows, {df[['cell_id', 'event_id']].drop_duplicates().shape[0]} unique keys")
    check("CSV and Parquet identical rows/keys", len(dfc) == len(df) and (dfc["cell_id"].to_numpy() == df["cell_id"].to_numpy()).all()
          and (dfc["event_id"].to_numpy() == df["event_id"].to_numpy()).all())
    check("Event set = config events", sorted(df.event_id.unique()) == sorted(e["event_id"] for e in evs),
          ", ".join(sorted(df.event_id.unique())))
    check("Labels are only 0/1", set(df.label.unique()) <= {0, 1})
    check("label_obs_utc NOT in model inputs", "label_obs_utc" not in schema.MODEL_INPUTS)
    check("No label/label_meta/meta/id column in model inputs", not (set(schema.MODEL_INPUTS) & set(schema.NEVER_FEATURES)),
          "model inputs: " + ", ".join(schema.MODEL_INPUTS))
    check("s2_water_confirmed, imerg_run, gfm_version, landcover_year_gap not model inputs",
          not ({"s2_water_confirmed", "imerg_run", "gfm_version", "landcover_year_gap"} & set(schema.MODEL_INPUTS)))
    feat_time = {c[0]: c[7] for c in schema.COLUMNS if c[1] == "feature"}
    check("Every non-trigger feature <= t0 (schema)", all("<= t0" in v for v in feat_time.values()))
    rng = {"slope_deg": (0, 90), "sm_surface_m3m3": (0, 1), "sm_rootzone_m3m3": (0, 1), "rain_event_24h_mm": (0, 1000),
           "rain_event_max30_mmhr": (0, 500), "rain_antecedent_7d_mm": (0, 1000), "rain_upstream_event_mm": (0, 1000),
           "dist_to_coast_m": (0, 1e6), "dist_to_ref_water_m": (0, 1e6), "flowacc_km2": (0, 1e5), "hand_m": (0, 5000)}
    bad = {c: int(((df[c] < lo) | (df[c] > hi)).sum()) for c, (lo, hi) in rng.items()}
    check("Physical ranges (incl. hand_m >= 0)", not any(bad.values()), json.dumps({k: v for k, v in bad.items() if v}) or "all in range")
    check("hand_ref_is_sea is 0/1 where hand_m present",
          set(df.loc[df.hand_m.notna(), "hand_ref_is_sea"].dropna().unique()) <= {0, 1})
    static = ["elevation_m", "slope_deg", "flowacc_km2", "twi", "hand_m", "landcover", "dist_to_coast_m", "x", "y", "block_id"]
    g = df[df.duplicated("cell_id", keep=False)].groupby("cell_id")[static].nunique(dropna=False)
    check("Static columns identical for a cell across events", bool((g <= 1).all().all()),
          f"{len(g)} cells appear in > 1 event")
    return checks


def event_checks(cfg, ev, d, vinfo, smap_end):
    checks = []
    eid = ev["event_id"]

    def check(name, ok, detail=""):
        checks.append((eid, name, "PASS" if ok else "FAIL", detail))

    check("Single t0 = config t0", d.t0_utc.nunique() == 1 and d.t0_utc.iloc[0] == ev["t0_utc"], d.t0_utc.iloc[0])
    n1, n0 = int((d.label == 1).sum()), int((d.label == 0).sum())
    ratio = cfg["labels"]["neg_per_pos_ratio"]
    short = sum(vinfo.get("shortfall", {}).values())
    check("Has positives and negatives", n1 > 0 and n0 > 0)
    check("Negatives = ratio x positives - documented shortfall", n0 == int(round(ratio * n1)) - short,
          f"pos={n1}, neg={n0}, target={int(round(ratio * n1))}, shortfall={short}")
    thr = cfg["labels"]["flood_fraction_threshold"]
    check(f"All positives flood_fraction >= {thr}", (d.loc[d.label == 1, "flood_fraction"] >= thr).all())
    check("All negatives flood_fraction == 0", (d.loc[d.label == 0, "flood_fraction"] == 0).all())
    p, n = d[d.label == 1], d[d.label == 0]
    dmin = float(cKDTree(p[["x", "y"]].to_numpy()).query(n[["x", "y"]].to_numpy())[0].min())
    buf = cfg["labels"]["negative_buffer_m"]
    check(f"Negatives >= {buf} m from positives of the same event (centres, tol. 25 m)", dmin >= buf - 25, f"min = {dmin:.0f} m")
    obs = sorted(d.label_obs_utc.unique())
    check("label_obs_utc populated, on the label date", (d.label_obs_utc.astype(str).str.len() > 0).all()
          and all(o.startswith(ev["gfm_post_date"]) for o in obs), ", ".join(obs))
    check("imerg_run / gfm_version / landcover_year_gap constant and as configured",
          d.imerg_run.nunique() == 1 and d.imerg_run.iloc[0] == ev["imerg_run"] and d.gfm_version.nunique() == 1
          and d.landcover_year_gap.nunique() == 1 and int(d.landcover_year_gap.iloc[0]) == ev["year"] - int(cfg["landcover"]["year"]),
          f"{d.imerg_run.iloc[0]} / {d.gfm_version.iloc[0]} / {int(d.landcover_year_gap.iloc[0])}")
    check("SMAP averaging window ends <= t0 (and within 3 h)", smap_end is not None and smap_end <= ev["t0"]
          and ev["t0"] - smap_end < timedelta(hours=3), f"window end {smap_end}")
    check("Rain columns populated", d[["rain_event_24h_mm", "rain_antecedent_7d_mm", "rain_upstream_event_mm"]].notna().all().all())
    return checks


def qa_variant(cfg, evs, fname, v, smap_by_event):
    proc = ROOT / cfg["paths"]["processed"]
    df = pd.read_parquet(proc / f"{fname}.parquet")
    dfc = pd.read_csv(proc / f"{fname}.csv", low_memory=False)
    checks = file_checks(cfg, df, dfc, evs)
    for ev in evs:
        d = df[df.event_id == ev["event_id"]]
        if d.empty:
            checks.append((ev["event_id"], "Event present", "FAIL", "no rows"))
            continue
        checks += event_checks(cfg, ev, d, v["events"][ev["event_id"]], smap_by_event.get(ev["event_id"]))
    miss = pd.DataFrame({"column": df.columns, "role": [schema.ROLE[c] for c in df.columns],
                         "n_missing": [int(df[c].isna().sum() + (df[c].astype(str) == "missing").sum()) for c in df.columns]})
    miss["pct_missing"] = (100 * miss.n_missing / len(df)).round(2)
    for ev in evs:
        miss[f"pct_missing_{ev['event_id']}"] = [round(100 * float(df.loc[df.event_id == ev["event_id"], c].isna().mean()), 2)
                                                for c in df.columns]
    return df, checks, miss


def main():
    cfg = load_config()
    evs = events(cfg)
    proc = ROOT / cfg["paths"]["processed"]
    summ = json.loads((proc / "build_summary.json").read_text(encoding="utf-8"))
    rain_qa = {e["event_id"]: json.loads(interim_path("rain_qa.json", e["event_id"]).read_text(encoding="utf-8")) for e in evs}
    terr_qa = json.loads((ROOT / "data" / "interim" / "terrain_qa.json").read_text(encoding="utf-8"))
    buf, thr = cfg["labels"]["negative_buffer_m"], cfg["labels"]["flood_fraction_threshold"]

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    smap_bt = __import__("09_build_table")
    wins = smap_windows()
    smap_by_event = {e["event_id"]: wins.get(smap_bt.smap_file(cfg, e["t0"]).name) for e in evs}

    all_errors = []
    L = ["# QA report — jeddah_flood_ml (version 2, multi-event)", "",
         f"Generated by `scripts/10_qa.py` at {datetime.now(timezone.utc):%Y-%m-%d %H:%M} UTC; tables built {summ['built_utc']}.",
         f"Default negative strategy: **{summ['default_strategy']}**. AOI cells (30 m): {summ['aoi_cells']:,}. "
         f"Label rule: new flood vs baseline; negatives >= {buf} m from any new-flood or baseline-water pixel; "
         f"positives need new-flood fraction >= {thr}.", "",
         "## Events", "",
         "| event | t0 (UTC) | GFM baseline → label | GFM version | IMERG run | WorldCover gap (y) | new-flood subpx (10 m) | baseline-water subpx | label flood on baseline water | positives | negative candidates | candidates < 10 m | S2 flag valid (AOI) |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for eid, e in summ["events"].items():
        L.append(f"| {eid} | {e['t0_utc']} | {e['gfm_pre_date']} → {e['gfm_post_date']} | {e['gfm_version']} | {e['imerg_run']} | "
                 f"{e['landcover_year_gap']:+d} | {e['new_flood_subpx']:,} | {e['baseline_water_subpx']:,} | "
                 f"{e['label_flood_on_baseline_water_subpx']:,} | {e['n_pos']:,} | {e['n_neg_candidates']:,} | "
                 f"{100 * e['candidates_below_10m_share']:.1f} % | {100 * e['s2_flag_valid_share_aoi']:.1f} % |")
    L += ["", "## File overview", "", "| file | strategy | event | rows | positives | negatives | negatives < 10 m | positives < 10 m | shortfall |",
          "|---|---|---|---|---|---|---|---|---|"]
    for fname, v in summ["variants"].items():
        for eid, ve in v["events"].items():
            L.append(f"| `{fname}` | {v['strategy']} | {eid} | {ve['rows']} | {ve['n_pos']} | {ve['n_neg']} | "
                     f"{100 * ve['neg_below_10m_share']:.1f} % | {100 * ve['pos_below_10m_share']:.1f} % | {sum(ve['shortfall'].values())} |")
        L.append(f"| `{fname}` | {v['strategy']} | **all** | **{v['rows']}** | {v['n_pos']} | {v['n_neg']} | | | |")
    L.append("")
    console = []
    for k, (fname, v) in enumerate(summ["variants"].items(), 1):
        df, checks, miss = qa_variant(cfg, evs, fname, v, smap_by_event)
        errs = [f"{c[0]}: {c[1]}" for c in checks if c[2] == "FAIL"]
        all_errors += [f"{fname}: {e}" for e in errs]
        sha = hashlib.sha256((proc / f"{fname}.csv").read_bytes()).hexdigest()
        blocks = df.groupby("block_id").label.agg(["size", "sum"])
        L += [f"## {k}. `{fname}` ({v['strategy']})", "", f"CSV sha256 `{sha}`", "",
              f"**{'ALL CHECKS PASSED' if not errs else 'FAILED: ' + ', '.join(errs)}** ({len(checks)} checks)", "",
              "| scope | check | result | detail |", "|---|---|---|---|"] + [f"| {a} | {b} | {c} | {d} |" for a, b, c, d in checks]
        L += ["", f"- Spatial blocks (1 km): {len(blocks)}; with positives: {int((blocks['sum'] > 0).sum())}; "
                  f"with both classes: {int(((blocks['sum'] > 0) & (blocks['sum'] < blocks['size'])).sum())}",
              f"- Cells appearing in more than one event: {int(df.cell_id.duplicated(keep=False).sum())} rows; "
              f"cells positive in one event and negative in another: "
              f"{int((df.groupby('cell_id').label.nunique() > 1).sum())}", ""]
        for eid, ve in v["events"].items():
            if ve["shortfall"]:
                L += [f"**Shortfall {eid}** (positive-elevation bins with fewer eligible observed-dry candidates than the quota; "
                      "no sampling with replacement):", "", md_table(pd.DataFrame(
                          [{"bin": b, "missing_negatives": s, "candidates_available": ve["candidates_per_stratum"].get(b)}
                           for b, s in ve["shortfall"].items()])), ""]
        L += ["### Composition per event", "",
              md_table(pd.crosstab([df.event_id, df.elev_bin], df.label).rename(columns={0: "neg", 1: "pos"}).reset_index()), "",
              md_table(pd.crosstab([df.event_id, df.landcover_name], df.label).rename(columns={0: "neg", 1: "pos"}).reset_index()), "",
              "s2_water_confirmed (label-quality metadata, NOT a feature):", "",
              md_table(df.groupby(["event_id", "label"])["s2_water_confirmed"].agg(
                  n="size", confirmed=lambda s: int((s == 1).sum()), not_confirmed=lambda s: int((s == 0).sum()),
                  nan=lambda s: int(s.isna().sum())).reset_index()), "",
              "### Missing values per column", "", md_table(miss[miss.n_missing > 0]) if (miss.n_missing > 0).any() else "none", "",
              "`hand_m` NaN: low flats whose D8 path ends in an unresolved flat before reaching a channel or the sea. "
              "`slope_deg`/`twi` NaN: cells touching the sea. `s2_water_confirmed` NaN: < 5 of 9 S2 pixels valid.", ""]
        if k == 1:
            # Label check: very low positives vs the rest (0.6 m ~ median 2022 positive elevation). Copernicus
            # DEM absolute vertical accuracy is < 4 m LE90 (Product Handbook v5.0): descriptive only.
            rows = []
            for eid in summ["events"]:
                pp = df[(df.label == 1) & (df.event_id == eid)]
                for nm, g in (("elevation < 0.6 m", pp[pp.elevation_m < 0.6]), ("elevation >= 0.6 m", pp[pp.elevation_m >= 0.6]),
                              ("all positives", pp)):
                    s = g.s2_water_confirmed
                    rows.append({"event": eid, "positives": nm, "n": len(g), "s2_with_data": int(s.notna().sum()),
                                 "s2_water_confirmed_pct": round(100 * float((s == 1).sum() / max(s.notna().sum(), 1)), 1),
                                 "median_dist_to_coast_m": round(float(g.dist_to_coast_m.median())) if len(g) else None,
                                 "median_dist_to_ref_water_m": round(float(g.dist_to_ref_water_m.median())) if len(g) else None})
            L += ["### Label check: positives below vs above 0.6 m", "", md_table(pd.DataFrame(rows)), "",
                  "The 0.6 m split is well inside Copernicus DEM vertical error (absolute < 4 m LE90, global mean "
                  "2.57 m; relative < 2 m on slopes <= 20 %; Product Handbook v5.0, 29 Nov 2022), so elevation "
                  "differences of this size near the coast are not physically resolved.", ""]
            num = [c for c in schema.MODEL_INPUTS if c not in ("landcover", "landcover_name", "elev_bin")]
            med = df.groupby(["event_id", "label"])[num].median(numeric_only=True).T
            med.columns = [f"{e}_{'pos' if l else 'neg'}" for e, l in med.columns]
            L += ["### Medians by event and label (sanity only — EDA is the team's task)", "",
                  md_table(med.round(3).reset_index().rename(columns={"index": "column"})), ""]
        console += [f"[{fname}] {c}  {a}: {b}  {d}" for a, b, c, d in checks]

    L += ["## Rain between the event and the label acquisition (QA only, not a feature)", ""]
    for eid, rq in rain_qa.items():
        L += [f"**{eid}** — IMERG V07 {rq['imerg_run']} run ({', '.join(rq['imerg_versions'])}) over the {rq['aoi_cells']} IMERG "
              f"cells intersecting the AOI (event 24 h: AOI mean {rq['event_24h_aoi_mean_mm']} mm, max {rq['event_24h_aoi_max_mm']} mm; "
              f"antecedent 7 d max {rq['antecedent_7d_aoi_max_mm']} mm). Label acquisition {rq['label_acq_utc']}.", "",
              "| from (UTC) | to (UTC) | AOI mean (mm) | AOI max (mm) |", "|---|---|---|---|"]
        L += [f"| {d['from_utc'][:16]} | {d['to_utc'][:16]} | {d['aoi_mean_mm']} | {d['aoi_max_mm']} |"
              for d in rq["rain_between_event_and_label"]]
        L += ["", f"Half-hourly files: {rq['n_files']} ({rq['first'][:16]} → {rq['last'][:16]} UTC); "
                  f"missing values in windows: {rq['nan_values_in_windows']}.", ""]
    L += ["## Terrain / hydrology QA", "", "```", json.dumps(terr_qa, indent=1), "```", ""]
    (ROOT / "docs" / "qa_report.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(console))
    if all_errors:
        sys.exit(f"QA FAILED: {all_errors}")


if __name__ == "__main__":
    main()
