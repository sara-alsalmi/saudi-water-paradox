# Data sources — Jeddah flood-susceptibility dataset (version 2: events of 21 Nov 2017 and 24 Nov 2022)

Every file under `data/raw/` is listed in `data/raw/manifest.csv` with its exact URL (including
OPeNDAP subset expressions), product version, SHA-256 checksum and download time. All parameters
live in `config.yaml` (event list under `events:`); `notebooks/01_data_preparation.ipynb` re-runs
the full pipeline.

| # | Source | Role | Version | Native resolution | Access | Licence |
|---|---|---|---|---|---|---|
| 1 | Copernicus GFM (CEMS Global Flood Monitoring) | **Labels** | GFM v3.1.0 archive (`processing:version`) for both events | 20 m, Equi7Grid AF020M | EODC STAC `https://stac.eodc.eu/api/v1`, collection `GFM`, no account | Copernicus open licence (free, attribution) |
| 2 | Sentinel-2 L2A (Earth Search, AWS) | **Label verification only** | processing baselines as recorded per file in the manifest | 10/20 m | STAC `https://earth-search.aws.element84.com/v1`, no account | Copernicus Sentinel data terms (free, attribution) |
| 3 | Copernicus DEM GLO-30 | Terrain, hydrology | GLO-30 public, 2021 release (AWS COG) | 1 arcsec (~30 m) | `s3://copernicus-dem-30m` (public HTTPS) | Copernicus DEM licence (free, attribution) |
| 4 | ESA WorldCover | Land cover (**2021 map used for both events**) | 2021 v200 | 10 m | `s3://esa-worldcover` (public HTTPS) | CC-BY 4.0 |
| 5 | GPM IMERG Final half-hourly (`GPM_3IMERGHH`) | Rainfall | V07 (V07B files) | 0.1° / 30 min | NASA GES DISC, Hyrax OPeNDAP subset, Earthdata login | NASA open data (cite) |
| 6 | SMAP L4 (`SPL4SMGP`) | Soil moisture | v008 (Vv8010) | 9 km EASE-Grid 2.0 / 3 h | NASA NSIDC DAAC, Hyrax OPeNDAP subset, Earthdata login | NASA open data, no restrictions (cite) |

**Not used:** FABDEM (only official access is a 10°×10° zip of 1.87 GB > 500 MB project cap) and
ERA5-Land (no CDS account). Event-level inventories (Global Flood Records, EM-DAT) and published
flood-point inventories are used only to **find and date** candidate events, never for labels.

## Events

| event_id | t0 (UTC) | Rain (IMERG, AOI) | GFM baseline → label scene | Label lag | Sentinel-2 pre / post | Status |
|---|---|---|---|---|---|---|
| `JED_20171121` | 2017-11-21 00:00 | 04:00–07:30 UTC, peak 05:30; 24 h AOI mean 23 mm, max 45 mm | 2017-11-13 → **2017-11-25 15:28:43** | ~4.4 days | 18 Nov (fallback 13 Nov) / 23 Nov | **included** |
| `JED_20221124` | 2022-11-24 00:00 | 24 h AOI mean 58 mm, max 81 mm | 2022-11-17 → **2022-11-29 15:29:57** | ~5 days | 22 Nov / 27 Nov | **included** (unchanged from v1) |

Event dates were confirmed from news reports and, for 2017 and 2022, from the Global Flood Records
inventory (Dartmouth Flood Observatory archive, Zenodo v0.9.0, doi:10.5281/zenodo.19288171,
records #4542 and #5329, whose polygons intersect the AOI). 2017 sources: Al Jazeera, 22 Nov 2017
(https://www.aljazeera.com/amp/news/2017/11/22/floods-hit-saudi-arabia); The Watchers, 21 Nov 2017
(https://watchers.news/2017/11/21/heavy-rain-wreaks-havoc-across-jeddah-saudi-arabia/). Rain timing
is from IMERG half-hourly over the AOI. t0 is 00:00 UTC of the event day for every event, so the
trigger window (`rain_event_*`) is always 24 h.

### Evaluated and rejected events (2026-11-14)

| Event | Reason |
|---|---|
| **9 Dec 2025** (`JED_20251209`) | **Evaluated and rejected — labels indistinguishable from noise.** GFM (v4.1.0) label scene 2025-12-13, baseline 2025-12-01. (1) Between two *dry* scenes (2025-11-19 → 2025-12-01, no rain) GFM already produced **4,360** "new" flood pixels (20 m), against **17,258** new-flood pixels after the event, so a large share of the signal is likely noise. (2) Two of the three densest new-flood zones (≈21.35°N 39.24°E and ≈21.83°N 39.10°E) look like scattered speckle on dry sand next to large "baseline water" patches (21,915 px of pre-existing GFM water). (3) Sentinel-2 shows water (MNDWI > 0) in only **4.9 %** of new-flood pixels vs **1.8 %** of observed-dry pixels (post-event S2 is 7–9 days after the event). Also: IMERG Final not yet published (Late run only). Config entry kept commented out in `config.yaml`; maps in `docs/checks/JED_20251209/`. |
| Early Jan 2023 | Rain on several days (29 Dec, 1 Jan, 3 Jan); the baseline scene (23 Dec 2022) was still flooded from Nov 2022 (37,527 px), and only 4,609 new pixels appeared — close to the noise level. |
| 6 Jan 2025 | Clean baseline (5 Jan) but the next GFM scene over the AOI is 17 Jan, **11 days** after the rain: much of the water would have drained, turning truly flooded cells into false negatives. |
| 17 Nov 2015 | No GFM scene over the AOI between 12 Oct and 1 Dec 2015. |
| 7 Nov 2026 | Light rain (AOI mean 3 mm, Late run), no confirmed flood report; the storm of 10–14 Nov 2026 was still ongoing at the time of the review. Can be re-evaluated later. |
| 25 Nov 2009, 26 Jan 2011 | Before Sentinel-1 (no GFM). The Global Flood Database (MODIS, 913 events 2000–2018) contains **no** Saudi event; SMAP starts 2015-03-31. The only label sources would be published flood-point inventories (e.g. IKONOS-based, Youssef et al. 2016), a different and incompatible label source. |

## Required attributions and citations

Include these in the report / any figure or derived product:

- **GFM labels:** "Contains modified Copernicus Service information 2017, 2022" — Copernicus Emergency
  Management Service, Global Flood Monitoring (GFM) product. Also: "Contains modified Copernicus
  Sentinel data 2017, 2022" (GFM is derived from Sentinel-1).
- **Sentinel-2:** "Contains modified Copernicus Sentinel data 2017, 2022".
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
- **Event inventory (dating only):** Global Flood Records v0.9.0 (Dartmouth Flood Observatory
  archive), Zenodo, doi:10.5281/zenodo.19288171.

## Processing summary

- **Reference grid:** 30 m, EPSG:32637 (UTM 37N), origin snapped to 60 m so the 10 m label/S2
  sub-grid nests exactly (3×3 subpixels per cell). AOI 39.05–39.40°E, 21.25–21.85°N; only cells
  whose centre lies inside the bbox are used. The grid is the same for every event, so `cell_id`
  identifies the same place in both events.
- **Time:** all timestamps UTC. t0 per event (table above). Non-trigger predictors are ≤ t0. Event
  rainfall [t0, t0 + 24 h) is used as the trigger (susceptibility design).
- **Static vs per-event columns:** terrain/hydrology, land cover, `elev_bin` and `dist_to_coast_m`
  are identical for a cell in both events (checked in QA). Per event: label, `dist_to_ref_water_m`
  (from that event's GFM reference-water layer), rain, soil moisture, `s2_water_confirmed`.
- **DEM / hydrology:** tiles mosaicked, sea (DEM = 0 m and missing ocean tiles) set to no-data.
  Stage A (90 m): D8 on the 38–41°E, 20–23°N extent, reverse flow-graph search for every cell
  draining into the AOI → **7,198 km²** catchment+AOI, bbox 39.02–40.52°E, 21.24–22.23°N (does
  not touch the extent edge). Stage B (30 m) on that bbox + 2 km: fill pits → fill depressions →
  resolve flats → D8 → accumulation (pysheds 0.5), Horn slope, TWI, HAND, and — **per event** —
  flow-weighted upstream mean event rainfall. **HAND drainage = channel cells (upstream area ≥ 1 km²) OR coastal outlet
  cells** (land cells whose D8 receiver is the sea, or that have no receiver and touch the sea):
  the sea is treated as a drainage outlet. HAND = conditioned elevation of the cell minus that of
  the first drainage cell on its D8 path; `hand_ref_is_sea = 1` when that cell is a coastal outlet.
- **Rainfall:** IMERG `precipitation` is mm/hr per 30-min slot → depth = rate × 0.5. Fill
  (-9999.9) → NaN; any missing half-hour in a window aborts the build (none were missing).
  Download window per event: [t0 − 7 d, label date + 1 d). Both events use the **Final** run
  (`imerg_run = Final`).
- **Soil moisture:** per event, the SMAP L4 granule whose 3-h window ends exactly at t0:
  `SMAP_L4_SM_gph_20171120T223000_Vv8010_001` (2017) and `SMAP_L4_SM_gph_20221123T223000_Vv8010_001`
  (2022). Fill → NaN; cells without a value would take the nearest valid EASE cell
  (`sm_filled_from_neighbor = 1`) — 0 rows needed it in either event.
- **Land cover:** majority (mode) of 10 m WorldCover pixels per 30 m cell; class 0 → NaN.
  `landcover_year_gap` = event year − 2021 (−4 for 2017, +1 for 2022).
- **Distances:** coast = connected sea body (DEM = 0 m, ≥ 1 km²); reference water = GFM reference
  water mask values 1 (permanent) and 2 (seasonal).

## Label construction and its limitations

### GFM encoding and version (important)
- Both events' archive scenes are processed with **GFM v3.1.0** (STAC `processing:version`;
  column `gfm_version`), while the encoding documentation consulted is the GFM Product Definition
  Document **v4.0** (Oct 2025, Tables 19–20, §4.4–4.5). Values were therefore **verified directly
  in the files**: `ensemble_flood_extent` {0, 1, 255}, `exclusion_mask` {0, 1, 255},
  `reference_water_mask` {0, 1, 2, 255}, `advisory_flags` {0, 1, 255}; no-data tag 255 in all layers.
- **Excluded pixels are coded 0 ("no flood") in `ensemble_flood_extent`, not 255.** Read naively,
  1,316 km² of excluded area (mostly the urban core) would become "observed dry" negatives. The
  pipeline applies the exclusion mask explicitly: a subpixel is valid only if observed (≠ 255)
  AND `exclusion_mask == 0`. Excluded share of the observed AOI: 66.2 % (2022), 66.3 % (2017).
- Overlapping scenes on the label date (two Equi7 tiles, one Sentinel-1 pass) never disagreed
  (0 conflicting pixels) in either event; merge rule: flooded if any valid scene is flooded.
  Advisory flags on new-flood pixels: none.

### "New flood vs baseline" rule (version 2, approved 2026-11-14, all events)
Each event has a **baseline** GFM scene from the same orbit ~1–2 weeks before the event. A 10 m
subpixel flooded on the label scene **and also flooded on the baseline** is *baseline water*:
it is neither a positive nor a negative, and it counts inside the 250 m negative buffer.
- 2022: baseline (17 Nov) has **0** flood pixels → the rule changes nothing; the 2022 rows are
  identical to version 1 (all 32 v1 columns, same rows and order).
- 2017: baseline (13 Nov) has 3,854 flood pixels (20 m) — persistent GFM water outside the
  reference-water mask. 1,248 of them are still "flooded" on 25 Nov and are removed from the
  positives. Sentinel-2 shows **no** water in any of these 1,248 pixels (MNDWI > 0: 0 %), versus
  55.7 % of the new-flood pixels — the removal is justified. 2,323 AOI cells touch baseline water.

### Definitions
- `label = 1`: ≥ 50 % of the cell's 10 m subpixels are **new** flood on the event's label scene.
- `label = 0`: **observed with no water on the label scene** — all 9 subpixels observed, not
  excluded, not reference water, not baseline water, and ≥ 250 m from any new-flood or
  baseline-water subpixel. **It does not mean "safe".**
- Negatives are sampled **per event** (`labels.negative_strategy` in `config.yaml`), seed 42 for
  every event, ratio 1:1 with that event's positives:
  - **`stratified` (default → `jeddah_flood_ml.*`)**: proportional stratified random sample over
    WorldCover class × fixed elevation bin → 7,546 (2017) + 12,481 (2022) negatives.
  - **`elevation_matched` (variant → `jeddah_flood_ml_elev_matched.*`)**: negatives follow each
    event's positive elevation distribution (20 quantile bins; candidates outside the positive
    range are never drawn). Not enough observed-dry cells exist at the lowest elevations, so the
    lowest bins are short by **1,243 (2017)** and **2,967 (2022)** negatives (no sampling with
    replacement) → 6,303 + 9,514 negatives.
- **Elevation bias of negatives:** 38.4 % (2017) / 37.5 % (2022) of default negatives are below
  10 m, versus 95.8 % / 96.1 % of positives.

| | 2017 | 2022 | total |
|---|---|---|---|
| `jeddah_flood_ml` rows (pos / neg) | 15,092 (7,546 / 7,546) | 24,962 (12,481 / 12,481) | **40,054** |
| `jeddah_flood_ml_elev_matched` rows (pos / neg) | 13,849 (7,546 / 6,303) | 21,995 (12,481 / 9,514) | **35,844** |
| eligible negative candidates | 474,778 | 468,480 | |

### Known biases and limitations (must be stated in the report)
1. **Urban blind spot.** ~66 % of the observed AOI is in the GFM exclusion mask (SAR has no
   sensitivity in dense urban areas). Labels describe where water *remained in open terrain*, not
   where neighbourhoods flooded. Built-up positives vs negatives: 50 vs 2,125 (2017), 169 vs 3,545
   (2022) — a model will learn "built-up ⇒ dry" as an artefact of the sensor, not of hydrology.
2. **Label lag.** The label scenes are ~4.4 days (2017) and ~5 days (2022) after the rain. Labels
   capture residual ponding in depressions/sabkha, not maximum flood extent.
3. **Optical confirmation is partial.**
   - 2022: only 59.7 % of GFM flood pixels showed an NDWI rise > 0.1 and 42.6 % had post-event
     MNDWI > 0 (vs 11.9 % / 4.6 % for observed-dry pixels), i.e. **~40 % of flood pixels have no
     optical water signal** — plausibly wet soil or sabkha. `s2_water_confirmed` = 1 for 68 % of
     positives and 13 % of negatives. Sentinel-2 is 2 days *before* the label scene.
   - 2017: where Sentinel-2 is available, `s2_water_confirmed` = 1 for **95.3 %** of positives and
     1.0 % of negatives (S2 post scene 23 Nov, 2 days before the label scene).
   - **2017 — `s2_water_confirmed` is empty in the north of the city:** the post-event Sentinel-2
     tiles were cloudy there (23 Nov tile 37QDE > 10 % cloud; the 28 Nov fallback was dropped to stay
     within the download cap). NaN share: **96 % of 2017 rows north of 21.70°N**, 39 % south of it;
     overall 53.7 % of 2017 rows and 65.9 % of 2017 positives (only 2,571 of 7,546 positives have a
     value). Treat the 2017 flag as a southern sample, not a full check.
4. **Coastal concentration.** ~96 % of positives are in the 0–10 m elevation bin in both events,
   mostly on the southern coastal sabkha (≈21.27–21.38°N) and north of Obhur (≈21.78–21.83°N).
   Positives fall in only 114 of 1,823 one-km blocks of the default file. **The same places flood
   in both events:** 5,166 of the 7,546 positive cells of 2017 (68 %) are also positive in 2022;
   81 cells are positive in one event and negative in the other. Tidal influence cannot be
   excluded. `dist_to_coast_m` and `dist_to_ref_water_m` are provided for sensitivity tests.
5. **No rain between event and label.** IMERG over the AOI between the end of the event window and
   the label acquisition: 2017 — 22 Nov mean 0.02 mm (max 0.5 mm), 23–25 Nov 0.0 mm; 2022 — 25 Nov
   mean 0.1 mm (max 1.3 mm), 26–29 Nov 0.0 mm. Water seen on the label scenes is not from later rain.
6. **IMERG underestimates local extremes.** 2022 AOI event totals are 41–81 mm (mean 58 mm) at
   0.1°. The reported point value is **179 mm between 08:00 and 14:00 (local time, = 05:00–11:00
   UTC) on 24 Nov 2022**, attributed to Saudi Arabia's National Center for Meteorology (NCM):
   FloodList (article by R. Davies on the Jeddah floods), 25 Nov 2022,
   https://floodlist.com/asia/saudi-arabia-floods-jeddah-november-2022 (The Watchers, 24 Nov 2022,
   https://watchers.news/2022/11/24/record-breaking-rains-cause-devastating-flash-floods-in-jeddah-saudi-arabia/
   gives the same value for "Jeddah's southern area"). Neither names the station, and the original
   NCM statement could not be retrieved, so treat it as an NCM-reported point value, not a verified
   gauge record. 2017 AOI event totals are 14–45 mm (mean 23 mm).
   **Rain variables — USE WITH CAUTION.** Rain takes only 24 distinct values per event across the
   AOI, and in **both** events IMERG rain is *lower* at the coastal positive clusters (median
   positives vs negatives: 20 vs 24 mm in 2017, 42 vs 62 mm in 2022) — a spatial confound, not a
   physical "less rain ⇒ more flooding" effect. With two events the rain columns now also differ
   *between* events (2022 wetter), so a model can use them to tell the events apart. Applies to
   `rain_event_24h_mm`, `rain_event_max30_mmhr`, `rain_upstream_event_mm` (equal to the 24 h value
   in ~96–98 % of rows) and `rain_antecedent_7d_mm` (≤ 3.3 mm in both events).
7. **DSM, not DTM, with metre-level vertical error.** Copernicus DEM includes buildings and
   vegetation; flow paths and HAND inside the city are affected. 561 (2022) / 574 (2017) default rows
   have DSM elevation < 0 m (coastal; binned in "0–10"). Absolute vertical accuracy is specified as
   < 4 m (90 % linear error; global mean LE90 2.57 m) and relative accuracy < 2 m on slopes ≤ 20 %;
   local deviations occur, and shore-/coastlines were edited (Copernicus DEM Product Handbook v5.0,
   Airbus Defence and Space, 29 Nov 2022, p. 8, p. 10, §2.1 and Table 12:
   https://dataspace.copernicus.eu/sites/default/files/media/files/2024-06/geo1988-copernicusdem-spe-002_producthandbook_i5.0.pdf).
   About half the positives lie below 0.6 m (50 % in 2022, 60 % in 2017), so their elevation
   differences are within DEM error. Positives < 0.6 m vs ≥ 0.6 m, `s2_water_confirmed`:
   73.0 % vs 63.0 % (2022), 98.2 % vs 84.7 % (2017; see `docs/qa_report.md`).
8. **HAND near the coast.** The sea is a drainage outlet (`hand_ref_is_sea`). Remaining `hand_m`
   NaN: 1.45 % (2022) / 2.48 % (2017) of default rows — low flats whose D8 path ends in an
   unresolved flat. `hand_ref_is_sea` is strongly label-associated (coastal proxy) and is therefore
   a **metadata flag, not a model input** (`use_as_model_input = NO`).
9. **Land cover from 2021 for both events (`landcover_year_gap`).** WorldCover 2021 is 4 years after
   the 2017 event and 1 year before 2022. Jeddah's large demolition and redevelopment programme of
   unplanned districts (2021–2023) and new construction since 2017 are not reflected: some cells
   labelled built-up were not yet built in 2017, and some demolished areas still appear as built-up.
   No year-matched land-cover map was used.
10. **Only two events, artificial prevalence.** Two events with largely the same flood locations;
    the 1:1 class ratio is a sampling choice (true flooded share of observed area ≈ 1–2 %).
    Predicted probabilities are not calibrated flood probabilities.
