"""Step 7 — Final QA of every table variant in data/processed -> docs/qa_report.md (fails loudly)."""
from __future__ import annotations

import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, load_config  # noqa: E402
import schema  # noqa: E402


def md_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    out = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        out.append("| " + " | ".join("" if (isinstance(v, float) and np.isnan(v)) else str(v) for v in r) + " |")
    return "\n".join(out)


def qa_variant(cfg, fname, vinfo, summ):
    proc = ROOT / cfg["paths"]["processed"]
    df = pd.read_parquet(proc / f"{fname}.parquet")
    dfc = pd.read_csv(proc / f"{fname}.csv", low_memory=False)
    checks = []

    def check(name, ok, detail=""):
        checks.append((name, "PASS" if ok else "FAIL", detail))

    check("Columns match schema (order)", list(df.columns) == schema.NAMES)
    check("Key (cell_id, event_id) unique", not df.duplicated(["cell_id", "event_id"]).any(),
          f"{len(df)} rows, {df[['cell_id', 'event_id']].drop_duplicates().shape[0]} unique keys")
    check("CSV and Parquet identical rows/keys", len(dfc) == len(df) and (dfc["cell_id"].to_numpy() == df["cell_id"].to_numpy()).all())
    check("Single event_id / t0", df["event_id"].nunique() == 1 and df["t0_utc"].nunique() == 1,
          f"{df['event_id'].iloc[0]} / {df['t0_utc'].iloc[0]}")
    n1, n0 = int((df.label == 1).sum()), int((df.label == 0).sum())
    ratio = cfg["labels"]["neg_per_pos_ratio"]
    short = sum(vinfo.get("shortfall", {}).values())
    expected = int(round(ratio * n1)) - short
    check("Labels are only 0/1", set(df.label.unique()) <= {0, 1})
    check("Negatives = ratio x positives - documented shortfall", n0 == expected,
          f"pos={n1}, neg={n0}, target={int(round(ratio * n1))}, shortfall={short}")
    thr = cfg["labels"]["flood_fraction_threshold"]
    check(f"All positives flood_fraction >= {thr}", (df.loc[df.label == 1, "flood_fraction"] >= thr).all())
    check("All negatives flood_fraction == 0", (df.loc[df.label == 0, "flood_fraction"] == 0).all())
    p, n = df[df.label == 1], df[df.label == 0]
    dmin = float(cKDTree(p[["x", "y"]].to_numpy()).query(n[["x", "y"]].to_numpy())[0].min())
    buf = cfg["labels"]["negative_buffer_m"]
    check(f"Negatives >= {buf} m from positives (centres, tol. 25 m)", dmin >= buf - 25, f"min = {dmin:.0f} m")
    check("label_obs_utc populated", (df.label_obs_utc.astype(str).str.len() > 0).all(),
          ", ".join(sorted(df.label_obs_utc.unique())))
    check("label_obs_utc NOT in model inputs", "label_obs_utc" not in schema.MODEL_INPUTS)
    check("No label/label_meta/id column in model inputs", not (set(schema.MODEL_INPUTS) & set(schema.NEVER_FEATURES)),
          "model inputs: " + ", ".join(schema.MODEL_INPUTS))
    check("s2_water_confirmed not a model input", "s2_water_confirmed" not in schema.MODEL_INPUTS)
    feat_time = {c[0]: c[7] for c in schema.COLUMNS if c[1] == "feature"}
    check("Every non-trigger feature <= t0 (schema)", all("<= t0" in v for v in feat_time.values()))
    smap_rows = [r for r in csv.DictReader(open(ROOT / "data/raw/manifest.csv", encoding="utf-8"))
                 if r["source"] == "nasa_smap_l4_spl4smgp"]
    smap_end = smap_rows[0]["notes"].split("-> ")[1].split(";")[0] if smap_rows else ""
    check("SMAP averaging window ends <= t0", bool(smap_end) and datetime.fromisoformat(smap_end) <=
          datetime.fromisoformat(cfg["time"]["t0_utc"].replace("Z", "+00:00")), f"window end {smap_end}")
    rng = {"slope_deg": (0, 90), "sm_surface_m3m3": (0, 1), "sm_rootzone_m3m3": (0, 1), "rain_event_24h_mm": (0, 1000),
           "rain_event_max30_mmhr": (0, 500), "rain_antecedent_7d_mm": (0, 1000), "rain_upstream_event_mm": (0, 1000),
           "dist_to_coast_m": (0, 1e6), "dist_to_ref_water_m": (0, 1e6), "flowacc_km2": (0, 1e5), "hand_m": (0, 5000)}
    bad = {c: int(((df[c] < lo) | (df[c] > hi)).sum()) for c, (lo, hi) in rng.items()}
    check("Physical ranges (incl. hand_m >= 0)", not any(bad.values()), json.dumps({k: v for k, v in bad.items() if v}) or "all in range")
    check("hand_ref_is_sea is 0/1 where hand_m present",
          set(df.loc[df.hand_m.notna(), "hand_ref_is_sea"].dropna().unique()) <= {0, 1})

    miss = pd.DataFrame({"column": df.columns, "role": [schema.ROLE[c] for c in df.columns],
                         "n_missing": [int(df[c].isna().sum() + (df[c].astype(str) == "missing").sum()) for c in df.columns]})
    miss["pct_missing"] = (100 * miss.n_missing / len(df)).round(2)
    miss["pct_missing_pos"] = [round(100 * float(p[c].isna().mean()), 2) for c in df.columns]
    miss["pct_missing_neg"] = [round(100 * float(n[c].isna().mean()), 2) for c in df.columns]
    blocks = df.groupby("block_id").label.agg(["size", "sum"])
    return df, checks, miss, blocks


def main():
    cfg = load_config()
    proc = ROOT / cfg["paths"]["processed"]
    summ = json.loads((proc / "build_summary.json").read_text(encoding="utf-8"))
    rain_qa = json.loads((ROOT / "data" / "interim" / "rain_qa.json").read_text(encoding="utf-8"))
    terr_qa = json.loads((ROOT / "data" / "interim" / "terrain_qa.json").read_text(encoding="utf-8"))
    buf, thr = cfg["labels"]["negative_buffer_m"], cfg["labels"]["flood_fraction_threshold"]
    all_errors = []
    L = ["# QA report — jeddah_flood_ml", "",
         f"Generated by `scripts/10_qa.py` at {datetime.now(timezone.utc):%Y-%m-%d %H:%M} UTC; tables built {summ['built_utc']}.",
         f"Default negative strategy: **{summ['default_strategy']}**. AOI cells (30 m): {summ['aoi_cells']:,}. "
         f"Eligible negative candidates: {summ['n_neg_candidates']:,} (observed, not excluded, not reference water, "
         f">= {buf} m from flood); **{100 * summ['candidates_below_10m_share']:.1f} % of candidates are below 10 m**. "
         f"Cells with 0 < flood_fraction < {thr} (neither label): {summ['pos_cells_flood_frac_between_0_and_thr']:,}.", ""]
    L += ["## File overview", "", "| file | strategy | rows | positives | negatives | negatives < 10 m | positives < 10 m | shortfall |",
          "|---|---|---|---|---|---|---|---|"]
    for fname, v in summ["variants"].items():
        L.append(f"| `{fname}` | {v['strategy']} | {v['rows']} | {v['n_pos']} | {v['n_neg']} | "
                 f"{100 * v['neg_below_10m_share']:.1f} % | {100 * v['pos_below_10m_share']:.1f} % | {sum(v['shortfall'].values())} |")
    L.append("")
    console = []
    for k, (fname, v) in enumerate(summ["variants"].items(), 1):
        df, checks, miss, blocks = qa_variant(cfg, fname, v, summ)
        errs = [c[0] for c in checks if c[1] == "FAIL"]
        all_errors += [f"{fname}: {e}" for e in errs]
        sha = hashlib.sha256((proc / f"{fname}.csv").read_bytes()).hexdigest()
        L += [f"## {k}. `{fname}` ({v['strategy']})", "", f"CSV sha256 `{sha}`", "",
              f"**{'ALL CHECKS PASSED' if not errs else 'FAILED: ' + ', '.join(errs)}**", "",
              "| check | result | detail |", "|---|---|---|"] + [f"| {a} | {b} | {c} |" for a, b, c in checks]
        L += ["", f"- Spatial blocks (1 km): {len(blocks)}; with positives: {int((blocks['sum'] > 0).sum())}; "
                  f"with both classes: {int(((blocks['sum'] > 0) & (blocks['sum'] < blocks['size'])).sum())}",
              f"- DSM elevation < 0 m rows: {v['elevation_below_0_rows']}", ""]
        if v["shortfall"]:
            L += ["**Shortfall** (positive-elevation bins with fewer eligible observed-dry candidates than the quota; "
                  "no sampling with replacement):", "", md_table(pd.DataFrame(
                      [{"bin": b, "missing_negatives": s, "candidates_available": v["candidates_per_stratum"].get(b)}
                       for b, s in v["shortfall"].items()])), ""]
        L += ["### Composition", "", md_table(pd.crosstab(df.elev_bin, df.label).rename(columns={0: "neg", 1: "pos"}).reset_index()), "",
              md_table(pd.crosstab(df.landcover_name, df.label).rename(columns={0: "neg", 1: "pos"}).reset_index()), "",
              "s2_water_confirmed (label-quality metadata, NOT a feature):", "",
              md_table(df.groupby("label")["s2_water_confirmed"].agg(
                  n="size", confirmed=lambda s: int((s == 1).sum()), not_confirmed=lambda s: int((s == 0).sum()),
                  nan=lambda s: int(s.isna().sum())).reset_index()), "",
              "### Missing values per column", "", md_table(miss[miss.n_missing > 0]) if (miss.n_missing > 0).any() else "none", "",
              "`hand_m` NaN: low flats whose D8 path ends in an unresolved flat before reaching a channel or the sea. "
              "`slope_deg`/`twi` NaN: cells touching the sea. `s2_water_confirmed` NaN: < 5 of 9 S2 pixels valid.", ""]
        if k == 1:
            # Label check: very low positives vs the rest (0.6 m ~ median positive elevation). Copernicus
            # DEM absolute vertical accuracy is < 4 m LE90 (Product Handbook v5.0), so this split is
            # far below the DEM's error and is descriptive only.
            pp = df[df.label == 1]
            rows = []
            for nm, g in (("elevation < 0.6 m", pp[pp.elevation_m < 0.6]), ("elevation >= 0.6 m", pp[pp.elevation_m >= 0.6]),
                          ("all positives", pp)):
                s = g.s2_water_confirmed
                rows.append({"positives": nm, "n": len(g), "s2_with_data": int(s.notna().sum()),
                             "s2_water_confirmed_pct": round(100 * float((s == 1).sum() / max(s.notna().sum(), 1)), 1),
                             "median_dist_to_coast_m": round(float(g.dist_to_coast_m.median())),
                             "median_dist_to_ref_water_m": round(float(g.dist_to_ref_water_m.median()))})
            L += ["### Label check: positives below vs above 0.6 m", "",
                  md_table(pd.DataFrame(rows)), "",
                  "Very low positives are *more* often optically confirmed and closer to the sea; they are not a "
                  "weaker label subset. The 0.6 m split is well inside Copernicus DEM vertical error (absolute "
                  "< 4 m LE90, global mean 2.57 m; relative < 2 m on slopes <= 20 %; Product Handbook v5.0, "
                  "29 Nov 2022), so elevation differences of this size near the coast are not physically resolved.", ""]
            num = [c for c in schema.MODEL_INPUTS if c not in ("landcover", "landcover_name", "elev_bin")]
            med = df.groupby("label")[num].median(numeric_only=True).T
            med.columns = ["median_neg", "median_pos"]
            L += ["### Medians by label (sanity only — EDA is the team's task)", "",
                  md_table(med.round(3).reset_index().rename(columns={"index": "column"})), ""]
        console += [f"[{fname}] {b}  {a}  {c}" for a, b, c in checks]

    L += ["## Rain between the event and the label acquisition (QA only, not a feature)", "",
          f"IMERG V07 Final over the {rain_qa['aoi_cells']} IMERG cells intersecting the AOI "
          f"(event 24 Nov: AOI mean {rain_qa['event_24h_aoi_mean_mm']} mm, max {rain_qa['event_24h_aoi_max_mm']} mm):", "",
          "| from (UTC) | to (UTC) | AOI mean (mm) | AOI max (mm) |", "|---|---|---|---|"]
    L += [f"| {d['from_utc'][:16]} | {d['to_utc'][:16]} | {d['aoi_mean_mm']} | {d['aoi_max_mm']} |"
          for d in rain_qa["rain_between_event_and_label"]]
    L += ["", f"Half-hourly files: {rain_qa['n_files']} ({rain_qa['first'][:16]} → {rain_qa['last'][:16]} UTC); "
              f"missing values in windows: {rain_qa['nan_values_in_windows']}.", "",
          "## Terrain / hydrology QA", "", "```", json.dumps(terr_qa, indent=1), "```", ""]
    (ROOT / "docs" / "qa_report.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(console))
    if all_errors:
        sys.exit(f"QA FAILED: {all_errors}")


if __name__ == "__main__":
    main()
