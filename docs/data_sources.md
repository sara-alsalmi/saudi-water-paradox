# Data sources — Jeddah flood-susceptibility pilot (event of 24 Nov 2022)

Every file under `data/raw/` is listed in `data/raw/manifest.csv` with its exact URL (including
OPeNDAP subset expressions), product version, SHA-256 checksum and download time. All parameters
live in `config.yaml`; `notebooks/01_data_preparation.ipynb` re-runs the full pipeline.

| # | Source | Role | Version | Native resolution | Access | Licence |
|---|---|---|---|---|---|---|
| 1 | Copernicus GFM (CEMS Global Flood Monitoring) | **Labels** | GFM v3.1.0 archive (`processing:version`) | 20 m, Equi7Grid AF020M | EODC STAC `https://stac.eodc.eu/api/v1`, collection `GFM`, no account | Copernicus open licence (free, attribution) |
| 2 | Sentinel-2 L2A (Earth Search, AWS) | **Label verification only** | processing baseline 04.00 | 10/20 m | STAC `https://earth-search.aws.element84.com/v1`, no account | Copernicus Sentinel data terms (free, attribution) |
| 3 | Copernicus DEM GLO-30 | Terrain, hydrology | GLO-30 public, 2021 release (AWS COG) | 1 arcsec (~30 m) | `s3://copernicus-dem-30m` (public HTTPS) | Copernicus DEM licence (free, attribution) |
| 4 | ESA WorldCover | Land cover | 2021 v200 | 10 m | `s3://esa-worldcover` (public HTTPS) | CC-BY 4.0 |
| 5 | GPM IMERG Final half-hourly (`GPM_3IMERGHH`) | Rainfall | V07 (V07B files) | 0.1° / 30 min | NASA GES DISC, Hyrax OPeNDAP subset, Earthdata login | NASA open data (cite) |
| 6 | SMAP L4 (`SPL4SMGP`) | Soil moisture | v008 (Vv8010) | 9 km EASE-Grid 2.0 / 3 h | NASA NSIDC DAAC, Hyrax OPeNDAP subset, Earthdata login | NASA open data, no restrictions (cite) |

**Not used:** FABDEM (only official access is a 10°×10° zip of 1.87 GB > 500 MB project cap) and
ERA5-Land (no CDS account). Event-level inventories (Global Flood Records, EM-DAT) and published
flood-point inventories are not used for labels.

## Required attributions and citations

Include these in the report / any figure or derived product:

- **GFM labels:** "Contains modified Copernicus Service information 2022" — Copernicus Emergency
  Management Service, Global Flood Monitoring (GFM) product. Also: "Contains modified Copernicus
  Sentinel data 2022" (GFM is derived from Sentinel-1).
- **Sentinel-2:** "Contains modified Copernicus Sentinel data 2022".
- **Copernicus DEM:** "Produced using Copernicus WorldDEM-30 © DLR e.V. 2010-2014 and © Airbus
  Defence and Space GmbH 2014-2018 provided under COPERNICUS by the European Union and ESA;
  all rights reserved."
- **ESA WorldCover:** "© ESA WorldCover project 2021 / Contains modified Copernicus Sentinel data
  (2021) processed by ESA WorldCover consortium." Zanaga, D. et al. (2022) ESA WorldCover 10 m
  2021 v200, doi:10.5281/zenodo.7254221.
- **IMERG:** Huffman, G.J., Stocker, E.F., Bolvin, D.T., Nelkin, E.J., Tan, J. (2023), GPM IMERG
  Final Precipitation L3 Half Hourly 0.1 degree x 0.1 degree V07, Greenbelt, MD, Goddard Earth
  Sciences Data and Information Services Center (GES DISC), doi:10.5067/GPM/IMERG/3B-HH/07.
- **SMAP L4:** SMAP L4 Global 3-hourly 9 km EASE-Grid Surface and Root Zone Soil Moisture
  Geophysical Data, Version 8, NASA NSIDC DAAC, doi:10.5067/T5RUATAQREF8 (use the full author
  citation from the NSIDC landing page).

## Processing summary

- **Reference grid:** 30 m, EPSG:32637 (UTM 37N), origin snapped to 60 m so the 10 m label/S2
  sub-grid nests exactly (3×3 subpixels per cell). AOI 39.05–39.40°E, 21.25–21.85°N; only cells
  whose centre lies inside the bbox are used.
- **Time:** all timestamps UTC. t0 = 2022-11-24T00:00Z. Non-trigger predictors are ≤ t0. Event
  rainfall (24 Nov 00:00–25 Nov 00:00 UTC) is used as the trigger (susceptibility design).
- **DEM / hydrology:** tiles mosaicked, sea (DEM = 0 m and missing ocean tiles) set to no-data.
  Stage A (90 m): D8 on the 38–41°E, 20–23°N extent, reverse flow-graph search for every cell
  draining into the AOI → **7,198 km²** catchment+AOI, bbox 39.02–40.52°E, 21.24–22.23°N (does
  not touch the extent edge). Stage B (30 m) on that bbox + 2 km: fill pits → fill depressions →
  resolve flats → D8 → accumulation (pysheds 0.5), Horn slope, TWI, HAND, flow-weighted upstream
  mean event rainfall. **HAND drainage = channel cells (upstream area ≥ 1 km²) OR coastal outlet
  cells** (land cells whose D8 receiver is the sea, or that have no receiver and touch the sea):
  the sea is treated as a drainage outlet. HAND = conditioned elevation of the cell minus that of
  the first drainage cell on its D8 path; `hand_ref_is_sea = 1` when that cell is a coastal outlet.
- **Rainfall:** IMERG `precipitation` is mm/hr per 30-min slot → depth = rate × 0.5. Fill
  (-9999.9) → NaN; any missing half-hour in a window aborts the build (none were missing).
- **Soil moisture:** SMAP L4 granule `SMAP_L4_SM_gph_20221123T223000_Vv8010_001` (3-h mean
  21:00–24:00 UTC 23 Nov, i.e. ends exactly at t0). Fill → NaN; cells without a value would take
  the nearest valid EASE cell (`sm_filled_from_neighbor = 1`) — in this build 0 rows needed it.
- **Land cover:** majority (mode) of 10 m WorldCover pixels per 30 m cell; class 0 → NaN.
- **Distances:** coast = connected sea body (DEM = 0 m, ≥ 1 km²); reference water = GFM reference
  water mask values 1 (permanent) and 2 (seasonal).

## Label construction and its limitations

### GFM encoding and version (important)
- The 2022 archive is processed with **GFM v3.1.0** (STAC `processing:version`), while the
  encoding documentation consulted is the GFM Product Definition Document **v4.0** (Oct 2025,
  Tables 19–20, §4.4–4.5). Values were therefore **verified directly in the files**:
  `ensemble_flood_extent` {0, 1, 255}, `exclusion_mask` {0, 1, 255}, `reference_water_mask`
  {0, 1, 2, 255}, `advisory_flags` {0, 255 here}; no-data tag 255 in all layers.
- **Excluded pixels are coded 0 ("no flood") in `ensemble_flood_extent`, not 255.** Read naively,
  1,316 km² of excluded area (mostly the urban core) would become "observed dry" negatives. The
  pipeline applies the exclusion mask explicitly: a subpixel is valid only if observed (≠ 255)
  AND `exclusion_mask == 0`.
- Overlapping scenes on 29 Nov (two Equi7 tiles, 15:29:57 UTC) never disagreed (0 conflicting
  pixels); merge rule: flooded if any valid scene is flooded. Advisory flags on flood pixels: none.
- Pre-event baseline (17 Nov, same orbit): **0** flood pixels in the AOI.

### Definitions
- `label = 1`: ≥ 50 % of the cell's 10 m subpixels flagged flooded on 2022-11-29.
- `label = 0`: **observed with no water on 2022-11-29** — all 9 subpixels observed, not excluded,
  not reference water, and ≥ 250 m from any flooded subpixel. **It does not mean "safe".**
- Negatives (`labels.negative_strategy` in `config.yaml`), seed 42, from 468,480 eligible cells:
  - **`stratified` (default → `jeddah_flood_ml.*`)**: proportional stratified random sample over
    WorldCover class × fixed elevation bin, 1:1 with positives → 12,481 negatives.
  - **`elevation_matched` (optional variant → `jeddah_flood_ml_elev_matched.*`)**: negatives follow
    the positives' elevation distribution (20 quantile bins of positive elevation, quota = positives
    per bin; candidates outside the positive elevation range are never drawn). Half of the positives
    lie below 0.6 m and there are not enough observed-dry cells that low, so the 10 lowest bins are
    short by **2,967** negatives in total (no sampling with replacement) → **9,514 negatives**.
- **Elevation bias of negatives:** only **37.5 %** of the default negatives (and 37.5 % of all
  eligible candidates) are below 10 m, versus **96.1 %** of positives. In the elevation-matched
  variant, 95.6 % of negatives are below 10 m.

### Known biases and limitations (must be stated in the report)
1. **Urban blind spot.** 66 % of the observed AOI is in the GFM exclusion mask (SAR has no
   sensitivity in dense urban areas). Labels describe where water *remained in open terrain*, not
   where neighbourhoods flooded. Only 169 positives are built-up cells vs 3,545 negatives — a model
   will learn "built-up ⇒ dry" as an artefact of the sensor, not of hydrology.
2. **Five-day lag.** The flood peaked on 24 Nov; the Sentinel-1 label scene is 29 Nov 15:29 UTC.
   Labels capture residual ponding in depressions/sabkha, not maximum flood extent.
3. **Optical confirmation is partial.** In the 20 m label check, only 59.7 % of GFM flood pixels
   showed an NDWI rise > 0.1 and 42.6 % had post-event MNDWI > 0 (vs 11.9 % / 4.6 % for observed-dry
   pixels), i.e. **~40 % of flood pixels have no optical water signal** — plausibly wet soil or
   sabkha that is dark/smooth in SAR. At cell level, `s2_water_confirmed` = 1 for 68 % of
   positives (32 % not confirmed) and 13 % of negatives. Sentinel-2 is 2 days *before* the label
   scene (27 vs 29 Nov), so disagreement is not proof of error.
4. **Coastal concentration.** 96 % of positives (vs 37.5 % of default negatives) are in the 0–10 m elevation bin, mostly on the
   southern coastal sabkha (≈21.33–21.38°N) and north of Obhur (≈21.83°N); positives fall in only
   95 of 1,748 one-km blocks. Tidal influence cannot be excluded, although the 17 Nov baseline
   (same orbit) shows no flooding. `dist_to_coast_m` and `dist_to_ref_water_m` are provided so the
   team can test sensitivity.
5. **No rain between event and label.** IMERG over the AOI, 25 Nov 00:00 – 29 Nov 15:30 UTC:
   25 Nov mean 0.1 mm (max 1.3 mm); 26–29 Nov 0.0 mm. Water seen on 29 Nov is therefore not from
   later rain.
6. **IMERG underestimates the local extreme.** AOI event totals are 41–81 mm (mean 58 mm) at
   0.1°. The reported point value is **179 mm between 08:00 and 14:00 (local time, = 05:00–11:00
   UTC) on 24 Nov 2022**, attributed to Saudi Arabia's National Center for Meteorology (NCM):
   FloodList (article by R. Davies on the Jeddah floods), 25 Nov 2022,
   https://floodlist.com/asia/saudi-arabia-floods-jeddah-november-2022 (The Watchers, 24 Nov 2022,
   https://watchers.news/2022/11/24/record-breaking-rains-cause-devastating-flash-floods-in-jeddah-saudi-arabia/
   gives the same value for "Jeddah's southern area"). Neither names the station, and the original
   NCM statement could not be retrieved, so treat it as an NCM-reported point value, not a verified
   gauge record.
   **Rain variables — USE WITH CAUTION: single-event, spurious spatial correlation with the label.**
   Rain varies over only 24 distinct values across the AOI and IMERG rain is *lower* at the coastal
   positive clusters (median 42 mm) than at negatives (62 mm) — a spatial confound, not a physical
   "less rain ⇒ more flooding" effect. Applies to `rain_event_24h_mm`, `rain_event_max30_mmhr`,
   `rain_upstream_event_mm` and `rain_antecedent_7d_mm`.
7. **DSM, not DTM, with metre-level vertical error.** Copernicus DEM includes buildings and
   vegetation; flow paths and HAND inside the city are affected. 561 rows have DSM elevation
   < 0 m (coastal; binned in "0–10"). Absolute vertical accuracy is specified as < 4 m (90 %
   linear error; global mean LE90 2.57 m) and relative accuracy < 2 m on slopes ≤ 20 %; local
   deviations occur, and shore-/coastlines were edited (Copernicus DEM Product Handbook v5.0,
   Airbus Defence and Space, 29 Nov 2022, p. 8, p. 10, §2.1 and Table 12:
   https://dataspace.copernicus.eu/sites/default/files/media/files/2024-06/geo1988-copernicusdem-spe-002_producthandbook_i5.0.pdf).
   Half the positives lie below 0.6 m, so their elevation differences are within DEM error.
   Positives < 0.6 m vs ≥ 0.6 m: `s2_water_confirmed` 73.0 % vs 63.0 %, median distance to sea
   967 m vs 1,504 m (see `docs/qa_report.md`).
8. **HAND near the coast.** With channels only, `hand_m` was NaN for 9.5 % of rows (13.8 % of
   positives vs 5.2 % of negatives) because flow reached the sea without crossing a channel. The
   sea is now a drainage outlet (`hand_ref_is_sea`, 2,009 default rows = 1, of which 1,515
   positives). Remaining NaN: 1.45 % (1.67 % pos / 1.24 % neg) — low flats ~2.5 km inland whose
   D8 path ends in an unresolved flat. `hand_ref_is_sea` is strongly label-associated (coastal
   proxy) and is therefore a **metadata flag, not a model input** (`use_as_model_input = NO`).
9. **Single event, artificial prevalence.** One event; the 1:1 class ratio is a sampling choice
   (true flooded share of observed area ≈ 1–2 %). Predicted probabilities are not calibrated
   flood probabilities.
