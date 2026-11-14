# The Saudi Water Paradox — Jeddah flood dataset

*How can a water-scarce country experience devastating floods?*
This folder holds the pilot dataset for that project, built around **two Jeddah floods: 21 November 2017
and 24 November 2022**.

> **Version 2 (14 Nov 2026).** The dataset now covers **two storms** instead of one. Same file names, new
> rows for 2017; the 2022 rows are exactly the same as in version 1 (the old version is kept in git).
> New columns: `landcover_year_gap`, `imerg_run`, `gfm_version` (information only, not for training).
> New rule: water that was **already there before the storm** is no longer counted as flood.
> Other candidate storms (2009, 2011, 2015, Jan 2023, Jan 2025, Dec 2025) were checked and rejected —
> reasons in [docs/data_sources.md](docs/data_sources.md).

---

## 1. What is this dataset?

We divided the Jeddah area into small squares, each **30 m × 30 m** (about the size of a small sports
field), and took a sample of them. For each square we collected facts such as: how high is it, how
steep is it, is it buildings or bare ground, how much rain fell, how wet was the soil before the storm.

Each square also has an **answer** (the *label*):

- **1** = a satellite saw flood water in this square.
- **0** = a satellite looked at this square and saw **no** water.

**The question it helps answer:** *what makes a place more likely to collect flood water?*
Teammates will use it to train a machine-learning model (a program that learns patterns from examples).

> ⚠️ **0 does not mean "safe".** It only means "seen dry on the satellite image a few days after the storm"
> (25 November 2017 or 29 November 2022).

---

## 2. Where the data comes from

| Source | What it adds |
|---|---|
| **Copernicus GFM** (European radar flood monitoring) | The answer: where the satellite saw flood water. Radar sees through clouds. |
| **Copernicus DEM** (elevation map) | The shape of the land: height, slope, where water flows and collects. |
| **ESA WorldCover** (land-cover map) | What the ground is: buildings, bare land, shrubs, farms, water. |
| **NASA GPM IMERG** (rain satellites) | How much rain fell: on flood day, in the heaviest half hour, and in the week before. |
| **NASA SMAP** (soil-moisture satellite) | How wet the soil was just before the rain started. |

A sixth source, **Sentinel-2** colour photos from space, was used **only to double-check** the answers,
never for training. All sources are free and open; required credits are in [docs/data_sources.md](docs/data_sources.md).

---

## 3. How the data was collected and cleaned

1. **Chose the area:** a rectangle around Jeddah (39.05–39.40° E, 21.25–21.85° N).
2. **Chose the storms:** looked for every Jeddah flood since 2009 and kept only storms with a clean radar
   image a few days after the rain (2017 and 2022).
3. **Checked the answers first:** compared the radar flood map with colour photos from before and after
   each storm. The team reviewed and approved this before anything else ([docs/checks/](docs/checks/),
   one folder per storm). Water already visible on the radar image *before* the storm is left out.
4. **Downloaded only what we need** from the official sources, cut to our area. Every file is listed in
   [data/raw/manifest.csv](data/raw/manifest.csv) with its link, date and a *checksum* (a fingerprint
   that proves the file has not changed).
5. **Cleaned it:**
   - error codes (like −9999) became "empty" — we never invented values;
   - units were made consistent (e.g. rain in millimetres);
   - all maps were lined up on the same 30 m grid so every square matches across sources;
   - all times are in UTC (Saudi time = UTC + 3 h).
6. **Traced the valleys:** found every valley (*wadi*) that drains into Jeddah — about **7,200 km²**,
   reaching the Taif mountains — because flood water comes from there.
7. **No peeking at the future:** everything used for prediction is from **before** the storm began
   (00:00 UTC on 21 Nov 2017 or 24 Nov 2022), except the storm's own rain, which is the cause.
8. **Picked the examples, storm by storm:** all squares where new water was seen (1) and an equal number
   of squares seen dry (0), chosen randomly but balanced, at least 250 m away from any water.
9. **Quality-checked it:** 32 automatic checks per file, for the whole file and for each storm (no
   duplicates, no future information, sensible values). Results: [docs/qa_report.md](docs/qa_report.md).

---

## 4. What one row and the main columns mean

**One row = one 30 m × 30 m square in Jeddah, for one storm.** The same square can appear twice (once per
storm); `event_id` says which storm.

| Column | Meaning |
|---|---|
| `label` | **The answer:** 1 = water seen, 0 = seen dry |
| `event_id` | Which storm: `JED_20171121` or `JED_20221124` |
| `elevation_m` | Height above sea level (m) |
| `slope_deg` | How steep the ground is (degrees) |
| `flowacc_km2` | Area of land that drains into this square (km²) |
| `twi` | Wetness index: higher = water tends to collect here |
| `hand_m` | Height above the nearest stream or the sea (m) |
| `landcover_name` | Type of surface (built-up, bare, shrubland…) |
| `dist_to_coast_m` | Distance to the sea (m) |
| `dist_to_ref_water_m` | Distance to the nearest permanent or seasonal water (m) |
| `sm_surface_m3m3`, `sm_rootzone_m3m3` | Soil wetness just before the storm |
| `rain_event_24h_mm` | Rain on flood day (mm) ⚠️ |
| `rain_event_max30_mmhr` | Heaviest half hour of rain (mm per hour) ⚠️ |
| `rain_upstream_event_mm` | Average rain over all valleys draining into this square (mm) ⚠️ — in this storm almost always the same as `rain_event_24h_mm`, so not a separate factor |
| `rain_antecedent_7d_mm` | Rain in the week before the storm (mm) ⚠️ |
| `block_id` | A 1 km × 1 km "neighbourhood" ID, used to split data fairly for testing |

⚠️ **Rain columns: use with caution.** They are coarse (about 11 km), and in both storms the rain values are
by chance *lower* where the flooding was seen, which can mislead a model.

### 🚫 Do NOT use these columns for training
| Column | Why |
|---|---|
| `s2_water_confirmed` | Comes from photos taken *after* the flood; only for checking answer quality (for 2017 mostly empty in the north of the city, because of clouds) |
| `flood_fraction`, `label_source`, `label_obs_utc` | Describe the answer itself — using them is cheating |
| `hand_ref_is_sea` | Almost the same as "near the coast", so it leaks the answer |
| `sm_filled_from_neighbor` | Always 0, so it carries no information |
| `landcover_year_gap`, `imerg_run`, `gfm_version` | Describe the data sources of each storm; the same for every row of a storm |
| `sampling_stratum` | Only records how the square was sampled |
| `cell_id`, `event_id`, `x`, `y`, `lon`, `lat`, `block_id`, `t0_utc` | IDs and locations, not causes (use `event_id` and `block_id` to split the data) |

Full details for every column: [data_dictionary.csv](data_dictionary.csv) (see `use_as_model_input`).

---

## 5. The two files

Both are in `data/processed/`, as `.csv` (opens in Excel) and `.parquet` (faster for code). Same columns.

- **`jeddah_flood_ml.csv` — the main file for modeling.**
- **`jeddah_flood_ml_elev_matched.csv` — a check file.**
  In the main file, most flooded places are low and most dry places are higher, so a model can score
  well just by learning "low = flooded". In the matched file, flooded and dry places have **similar
  elevations**, so the model must rely on other factors.
  **Train on both and compare:** a big drop in the matched file means the model mostly relies on elevation.

**New with two storms:** train on one storm and test on the other (both ways round) to see whether the
model works on a storm it has never seen — see [docs/handover.md](docs/handover.md).

## 6. Final numbers

| File | Storm | Rows | Flooded (1) | Dry (0) |
|---|---|---|---|---|
| `jeddah_flood_ml` | 2017 | 15,092 | 7,546 | 7,546 |
| | 2022 | 24,962 | 12,481 | 12,481 |
| | **both** | **40,054** | 20,027 | 20,027 |
| `jeddah_flood_ml_elev_matched` | 2017 | 13,849 | 7,546 | 6,303 |
| | 2022 | 21,995 | 12,481 | 9,514 |
| | **both** | **35,844** | 20,027 | 15,817 |

The matched file has fewer dry squares because there are not enough dry places as low as the flooded ones.

---

## 7. Main limitations

1. **Radar cannot see inside the dense city.** About two-thirds of the area it observed was left out,
   mostly neighbourhoods. So the data shows where water **stayed in open land**, not which streets flooded.
2. **The picture is from 4–5 days later.** Storms on 21 Nov 2017 and 24 Nov 2022; radar images from
   25 Nov 2017 and 29 Nov 2022. We see water that **remained**, not the peak flood.
3. **Not all "water" spots look like water in colour photos:** in 2022 about 40% do not (maybe wet soil or
   salt flats, *sabkha*). In 2017 almost all checked spots do, but the north of the city could not be
   checked because of clouds.
4. **Most flooded spots are low and near the coast** (96% below 10 m) — and **mostly the same places in both
   storms**. The elevation map is only accurate to about 2–4 m, so tiny height differences near the sea
   are not reliable.
5. **Satellite rain is lower than measured:** our highest 2022 value is 81 mm, while the national weather
   centre reported 179 mm in 6 hours (as reported by FloodList).
6. **The elevation map includes buildings and trees**, not just bare ground.
7. **The land-cover map is from 2021 for both storms.** Jeddah changed a lot (new building since 2017,
   large demolitions of old districts in 2021–2023), so for 2017 some squares have the wrong surface type.
   The column `landcover_year_gap` shows how many years apart the map and the storm are.
8. **Only two storms**, and we chose a 50/50 mix of flooded and dry squares, so model scores are
   *relative*, not real-world probabilities.
9. ✅ **No rain fell between each storm and its radar image** (at most 1.3 mm), so the water seen came
   from that storm.

More detail: [docs/data_sources.md](docs/data_sources.md) · training advice: [docs/handover.md](docs/handover.md).

---

## 8. Re-run everything from scratch

**You need:** Python 3.13, internet, a free **NASA Earthdata** account, and about **1 GB** of free space.

**Step 1 — set up** (once, in a terminal inside this folder):
```
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
```
(For exactly the same versions, use `requirements.lock.txt` instead.)

To use the notebook, register the environment as a notebook *kernel* (the Python the notebook runs on):
```
.venv\Scripts\python -m ipykernel install --sys-prefix --name saudi-water-paradox --display-name "Python (saudi-water-paradox .venv)"
```

**Step 2 — log in to NASA** (once; you type your own password, nothing is saved in the project):
```
.venv\Scripts\python -c "import earthaccess; earthaccess.login(strategy='interactive', persist=True)"
```
Then on the Earthdata website: Profile → Applications → Authorized Apps → approve **NASA GESDISC DATA ARCHIVE**.

**Step 3 — run all steps in order.** Either open
[notebooks/01_data_preparation.ipynb](notebooks/01_data_preparation.ipynb), choose the kernel
**Python (saudi-water-paradox .venv)** and click *Run All*, or run:
```
.venv\Scripts\python scripts\01_probe_labels.py
.venv\Scripts\python scripts\02_download_dem.py
.venv\Scripts\python scripts\03_download_worldcover.py
.venv\Scripts\python scripts\04_download_imerg.py
.venv\Scripts\python scripts\05_download_smap.py
.venv\Scripts\python scripts\06_download_labels.py
.venv\Scripts\python scripts\07_derive_rain.py
.venv\Scripts\python scripts\08_derive_terrain.py
.venv\Scripts\python scripts\09_build_table.py
.venv\Scripts\python scripts\10_qa.py
```

**Good to know**
- All settings live in one file: [config.yaml](config.yaml). Don't change numbers inside the scripts.
- Re-running is safe: files already downloaded are not downloaded again.
- Results are identical every time (same numbers, same fingerprints in the quality report).
- Total download is about 690 MB; a full run takes roughly 30–50 minutes (an estimate, depending on internet and computer).

---

## Folder map

```
README.md                  ← this guide
config.yaml                ← all settings
requirements.txt           ← required libraries and versions
data_dictionary.csv        ← every column explained
data/raw/                  ← original downloads (re-created by the scripts) + manifest.csv
data/interim/              ← in-between files (re-created automatically)
data/processed/            ← the final tables
scripts/                   ← steps 01 to 10
notebooks/                 ← notebook that runs every step
docs/data_sources.md       ← sources, licences, credits, limitations
docs/qa_report.md          ← quality-check report
docs/handover.md           ← training guide for the team
docs/checks/<event_id>/    ← maps used to check the answers, one folder per storm (review only;
                              JED_20251209 = rejected storm, kept as evidence)
```
