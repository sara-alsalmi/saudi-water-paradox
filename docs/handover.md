# Handover — `jeddah_flood_ml` (pilot table)

**Files** (`data/processed/`, each as `.parquet` [preferred, typed] and `.csv` [same rows]):

| file | negatives | rows | label=1 | label=0 | use |
|---|---|---|---|---|---|
| `jeddah_flood_ml` | stratified (land cover × elevation bin) | 24,962 | 12,481 | 12,481 | **default** |
| `jeddah_flood_ml_elev_matched` | follow the positives' elevation distribution | 21,995 | 12,481 | 9,514 | optional robustness check |

Same columns in both. Definitions: `data_dictionary.csv` (see `use_as_model_input` and `caution`).
Checks: `docs/qa_report.md`. Sources, licences, attributions, limitations: `docs/data_sources.md`.

## What one row is
One 30 m grid cell (EPSG:32637) in the Jeddah AOI for the 24 Nov 2022 event. Key: `cell_id` + `event_id`.

- `label = 1`: GFM (Sentinel-1) saw flood water in ≥ 50 % of the cell on **29 Nov 2022**.
- `label = 0`: GFM observed the cell on 29 Nov and found **no water**, ≥ 250 m from any flood.
  It is **not** "safe" — just "observed dry on that date".

## Which columns to use
Only columns with `use_as_model_input = yes`:

- Static, ≤ t0: `elevation_m, slope_deg, flowacc_km2, twi, hand_m, landcover /
  landcover_name, elev_bin, dist_to_coast_m, dist_to_ref_water_m`
- Antecedent, ≤ t0: `rain_antecedent_7d_mm`*, `sm_surface_m3m3, sm_rootzone_m3m3`
- Trigger (event rain): `rain_event_24h_mm`*, `rain_event_max30_mmhr`*, `rain_upstream_event_mm`*

\* **All rain variables: use with caution — single-event, spurious spatial correlation with the
label (lower IMERG rain at coastal positives).** `rain_upstream_event_mm` is **nearly identical
to `rain_event_24h_mm` in this event** (equal in ~96 % of rows, because most cells' catchments lie
inside one IMERG cell) — **do not treat it as an independent variable.**

**Never use as features:** `hand_ref_is_sea` (metadata flag; strongly associated with the label —
coastal proxy), `sm_filled_from_neighbor` (constant 0 in this event — no information), `label_obs_utc`, `flood_fraction`, `label_source`, `sampling_stratum`,
and especially **`s2_water_confirmed` — NOT a feature, do not use for training.** It is post-event
label-quality metadata; use it only for a sensitivity check (e.g. evaluate on positives with
`s2_water_confirmed == 1`). Do not feed `cell_id, x, y, lon, lat, block_id` as predictors.
`landcover`/`landcover_name` are duplicates (keep one); `elev_bin` is a binned `elevation_m`.

## Spatial cross-validation with `block_id`
Neighbouring cells are near-duplicates; a random split leaks them into the test set. Group by
`block_id` (1 km × 1 km):

```python
from sklearn.model_selection import StratifiedGroupKFold
cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
for tr, te in cv.split(X, y, groups=df["block_id"]):
    ...
```

Check per-fold class counts; as a robustness check, also group by coarser blocks (e.g. 5 km).

## Missing values
`hand_m` 1.45 % (low inland flats with unresolved flow), `slope_deg` / `twi`
0.13 % (cells touching the sea), `s2_water_confirmed` 0.3 %. Impute inside the CV loop.

## Modeling warnings
1. **Label bias toward low coastal ground.** 96 % of positives are below 10 m (vs 37.5 % of the
   default negatives); half of the positives are below 0.6 m. Elevation, `dist_to_coast_m` and
   `dist_to_ref_water_m` will dominate.
   **Label check by elevation** (positives only): below 0.6 m — n = 6,271, `s2_water_confirmed`
   73.0 %, median distance to sea 967 m; at/above 0.6 m — n = 6,210, 63.0 %, 1,504 m. The very low
   positives are *not* a weaker subset (they are more often optically confirmed). But a 0.6 m
   difference is far inside the DEM's error: Copernicus DEM absolute vertical accuracy is < 4 m
   (90 % linear error; global mean 2.57 m) and relative accuracy < 2 m on slopes ≤ 20 %
   (*Copernicus DEM Product Handbook*, v5.0, Airbus Defence and Space, 29 Nov 2022, §2.1 and p. 10,
   https://dataspace.copernicus.eu/sites/default/files/media/files/2024-06/geo1988-copernicusdem-spe-002_producthandbook_i5.0.pdf),
   and coastlines were edited during production (p. 8). Sub-metre elevation near the coast is
   therefore not physically meaningful; a model splitting at such values is fitting DEM noise
   and coastal editing, not real terrain.
2. **Few urban positives.** Only 169 positives are built-up cells (vs 3,545 negatives): GFM
   cannot see water in the dense city (66 % of the observed AOI is excluded). "Built-up ⇒ dry" is
   a sensor artefact, not a finding.
3. **Few independent positive areas.** Positives fall in only **95 of 1,748** one-km blocks (two
   coastal clusters); effective sample size is far smaller than 12,481.
4. **Rain caveat.** All rain variables come from one event at ~11 km resolution (24 distinct
   values) and are *lower* at the coastal positives — any learned rain effect is likely spurious.
   `rain_upstream_event_mm` is nearly identical to `rain_event_24h_mm` (~96 % of rows equal); do
   not treat it as an independent variable.
   Report results with and without rain variables.
5. **Baselines are mandatory.** Compare every model against an **elevation-only baseline**
   (e.g. logistic regression on `elevation_m`). If a model barely beats it, it has learned
   "low = flooded". Also evaluate on `jeddah_flood_ml_elev_matched`, where elevation alone
   separates the classes much less.
6. **Always use `block_id` grouped CV** (above). Random K-fold scores will be inflated.

## Limitations to state in the report
1. Urban blind spot (above). 2. **5-day lag:** labels are residual water on 29 Nov, not peak extent.
3. **~40 % of GFM flood pixels show no optical water signal** (Sentinel-2; possibly wet soil /
sabkha); `s2_water_confirmed` = 1 for 68 % of positives. 4. Coastal concentration (above).
5. IMERG event rain ≤ 81 mm vs an NCM-reported 179 mm (08:00–14:00 local, via FloodList,
25 Nov 2022). 6. Copernicus DEM is a **DSM** (buildings included). 7. **Single event, 1:1
sampling** → probabilities are not calibrated; report ROC-AUC / PR-AUC and the sampling design.
8. No rain fell between the event and the label date (max 1.3 mm on 25 Nov).
9. The elevation-matched file is short of 2,967 negatives (not enough observed-dry cells below
0.6 m), so its class ratio is 1 : 0.76.
