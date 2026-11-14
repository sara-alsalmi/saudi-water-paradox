# Handover — `jeddah_flood_ml` (version 2: two events)

**Files** (`data/processed/`, each as `.parquet` [preferred, typed] and `.csv` [same rows]):

| file | negatives | event | rows | label=1 | label=0 | use |
|---|---|---|---|---|---|---|
| `jeddah_flood_ml` | stratified (land cover × elevation bin) | JED_20171121 | 15,092 | 7,546 | 7,546 | **default** |
| | | JED_20221124 | 24,962 | 12,481 | 12,481 | |
| | | **total** | **40,054** | 20,027 | 20,027 | |
| `jeddah_flood_ml_elev_matched` | follow each event's positive elevation distribution | JED_20171121 | 13,849 | 7,546 | 6,303 | optional robustness check |
| | | JED_20221124 | 21,995 | 12,481 | 9,514 | |
| | | **total** | **35,844** | 20,027 | 15,817 | |

Same columns in both. Definitions: `data_dictionary.csv` (see `use_as_model_input` and `caution`).
Checks: `docs/qa_report.md`. Sources, licences, attributions, limitations, rejected events:
`docs/data_sources.md`. The 2022 rows are identical to version 1 (git keeps v1).

## What one row is
One 30 m grid cell (EPSG:32637) in the Jeddah AOI **for one event**. Key: `cell_id` + `event_id`.
The same `cell_id` is the same place in both events (static columns are identical).

- `label = 1`: GFM (Sentinel-1) saw **new** flood water (not present on the pre-event baseline
  scene) in ≥ 50 % of the cell on the label date — 25 Nov 2017 or 29 Nov 2022.
- `label = 0`: GFM observed the cell on the label date and found **no water**, ≥ 250 m from any
  new flood or pre-existing water. It is **not** "safe" — just "observed dry on that date".

## Which columns to use
Only columns with `use_as_model_input = yes`:

- Static, ≤ t0: `elevation_m, slope_deg, flowacc_km2, twi, hand_m, landcover /
  landcover_name, elev_bin, dist_to_coast_m, dist_to_ref_water_m`
- Antecedent, ≤ t0: `rain_antecedent_7d_mm`*, `sm_surface_m3m3, sm_rootzone_m3m3`
- Trigger (event rain): `rain_event_24h_mm`*, `rain_event_max30_mmhr`*, `rain_upstream_event_mm`*

\* **All rain variables: use with caution.** Only 24 distinct values per event, and in **both**
events IMERG rain is *lower* at the coastal positives (spurious spatial correlation). With two
events the rain level also identifies the event (2022 was much wetter). `rain_upstream_event_mm`
is **nearly identical to `rain_event_24h_mm`** (equal in ~96–98 % of rows) — **do not treat it as
an independent variable.**

**Never use as features:** `hand_ref_is_sea` (coastal proxy), `sm_filled_from_neighbor` (constant
0), `landcover_year_gap`, `imerg_run`, `gfm_version` (constant per event — they would only encode
the event), `label_obs_utc`, `flood_fraction`, `label_source`, `sampling_stratum`, and especially
**`s2_water_confirmed` — NOT a feature, do not use for training.** It is post-event label-quality
metadata; use it only for a sensitivity check (e.g. evaluate on positives with
`s2_water_confirmed == 1`). For 2017 it is **empty in the north of the city** (cloudy post-event
images; 96 % NaN north of 21.70°N) and available for only 2,571 of 7,546 positives.
Do not feed `cell_id, event_id, t0_utc, x, y, lon, lat, block_id` as predictors.
`landcover`/`landcover_name` are duplicates (keep one); `elev_bin` is a binned `elevation_m`.

## How to split the data

### 1. By event (new in version 2) — train on one event, test on the other
This answers *"does a model trained on one storm predict a different storm?"*. Run it **in both
directions** and report both scores:

```python
df = pd.read_parquet("data/processed/jeddah_flood_ml.parquet")
for train_ev, test_ev in [("JED_20221124", "JED_20171121"), ("JED_20171121", "JED_20221124")]:
    tr, te = df[df.event_id == train_ev], df[df.event_id == test_ev]
    model.fit(tr[FEATURES], tr["label"])
    score = roc_auc_score(te["label"], model.predict_proba(te[FEATURES])[:, 1])
```

Be aware: **the same places flood in both storms** (68 % of the 2017 positive cells are also
positive in 2022). An event split therefore still shares locations between train and test. For a
stricter test, also hold out space: drop from the test event every `block_id` that is used in
training, or combine both (`groups = block_id`, test only on the other event). The rain columns
differ in level between events (58 mm vs 23 mm AOI mean), so compare scores with and without rain.

### 2. Spatial cross-validation with `block_id` (within or across events)
Neighbouring cells are near-duplicates; a random split leaks them into the test set. Group by
`block_id` (1 km × 1 km). A cell appears in both events with the same `block_id`, so grouping by
block also keeps the two copies of a cell on the same side:

```python
from sklearn.model_selection import StratifiedGroupKFold
cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
for tr, te in cv.split(X, y, groups=df["block_id"]):
    ...
```

Check per-fold class counts and per-fold event mix; as a robustness check, also group by coarser
blocks (e.g. 5 km).

## Missing values
`hand_m` 1.45 % (2022) / 2.48 % (2017) (low inland flats with unresolved flow), `slope_deg` / `twi`
~0.13–0.15 % (cells touching the sea), `s2_water_confirmed` 0.3 % (2022) / 53.7 % (2017,
clouds — see above). Impute inside the CV loop (never `s2_water_confirmed`, which is not a feature).

## Modeling warnings
1. **Label bias toward low coastal ground.** ~96 % of positives are below 10 m in both events
   (vs ~38 % of the default negatives); 50 % (2022) and 60 % (2017) of positives are below 0.6 m.
   Elevation, `dist_to_coast_m` and `dist_to_ref_water_m` will dominate.
   **Label check by elevation** (positives, `s2_water_confirmed` where available): 2022 — below
   0.6 m 73.0 %, at/above 0.6 m 63.0 %; 2017 — 98.2 % vs 84.7 %. The very low positives are *not*
   a weaker subset. But a 0.6 m difference is far inside the DEM's error: Copernicus DEM absolute
   vertical accuracy is < 4 m (90 % linear error; global mean 2.57 m) and relative accuracy < 2 m
   on slopes ≤ 20 % (*Copernicus DEM Product Handbook*, v5.0, Airbus Defence and Space,
   29 Nov 2022, §2.1 and p. 10,
   https://dataspace.copernicus.eu/sites/default/files/media/files/2024-06/geo1988-copernicusdem-spe-002_producthandbook_i5.0.pdf),
   and coastlines were edited during production (p. 8). Sub-metre elevation near the coast is not
   physically meaningful; a model splitting at such values is fitting DEM noise.
2. **Few urban positives.** Built-up positives vs negatives: 169 vs 3,545 (2022), 50 vs 2,125
   (2017). GFM cannot see water in the dense city (~66 % of the observed AOI is excluded).
   "Built-up ⇒ dry" is a sensor artefact, not a finding.
3. **Few independent positive areas.** Positives fall in only **114 of 1,823** one-km blocks (two
   coastal clusters, largely the same in both events); effective sample size is far smaller than
   20,027.
4. **Rain caveat** (above). Report results with and without rain variables.
5. **Land cover is from 2021 for both events** (`landcover_year_gap` = −4 for 2017). Districts
   demolished in 2021–2023 and areas built after 2017 are misrepresented, especially for 2017.
6. **Baselines are mandatory.** Compare every model against an **elevation-only baseline**
   (e.g. logistic regression on `elevation_m`), in the event split as well. If a model barely beats
   it, it has learned "low = flooded". Also evaluate on `jeddah_flood_ml_elev_matched`.
7. **Never use random K-fold** (scores will be inflated): use the event split and/or `block_id`
   grouped CV.

## Limitations to state in the report
1. Urban blind spot (above). 2. **Label lag:** labels are residual water ~4.4 days (2017) and
~5 days (2022) after the rain, not peak extent. 3. **Optical confirmation:** in 2022 ~40 % of GFM
flood pixels show no optical water signal (`s2_water_confirmed` = 1 for 68 % of positives); in 2017
95 % of positives with S2 data are confirmed, but the north is not covered. 4. Coastal
concentration; same locations in both events. 5. IMERG event rain ≤ 81 mm (2022) vs an
NCM-reported 179 mm (08:00–14:00 local, via FloodList, 25 Nov 2022). 6. Copernicus DEM is a **DSM**
(buildings included). 7. **Two events, 1:1 sampling per event** → probabilities are not calibrated;
report ROC-AUC / PR-AUC and the sampling design. 8. No rain fell between either event and its label
date (max 1.3 mm). 9. The elevation-matched file is short of 1,243 (2017) and 2,967 (2022)
negatives (not enough observed-dry cells below ~0.3–0.6 m). 10. WorldCover 2021 for both events;
2021–2023 demolitions not reflected. 11. Pre-existing water on the 2017 baseline (persistent GFM
water outside the reference mask) is excluded from both classes. 12. Candidate events evaluated and
rejected (Dec 2025, Jan 2025, Jan 2023, Nov 2015, 2009/2011): see `docs/data_sources.md`.
