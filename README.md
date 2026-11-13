# The Saudi Water Paradox — Jeddah flood dataset

*How can a water-scarce country experience devastating floods?*
This folder holds the pilot dataset for that project, built around the **Jeddah flood of 24 November 2022**.

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

> ⚠️ **0 does not mean "safe".** It only means "seen dry on 29 November 2022".

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
2. **Checked the answers first:** compared the radar flood map with colour photos from before and after
   the storm. The team reviewed and approved this before anything else ([docs/checks/](docs/checks/)).
3. **Downloaded only what we need** from the official sources, cut to our area. Every file is listed in
   [data/raw/manifest.csv](data/raw/manifest.csv) with its link, date and a *checksum* (a fingerprint
   that proves the file has not changed).
4. **Cleaned it:**
   - error codes (like −9999) became "empty" — we never invented values;
   - units were made consistent (e.g. rain in millimetres);
   - all maps were lined up on the same 30 m grid so every square matches across sources;
   - all times are in UTC (Saudi time = UTC + 3 h).
5. **Traced the valleys:** found every valley (*wadi*) that drains into Jeddah — about **7,200 km²**,
   reaching the Taif mountains — because flood water comes from there.
6. **No peeking at the future:** everything used for prediction is from **before** the storm began
   (24 Nov 2022, 00:00 UTC), except the storm's own rain, which is the cause.
7. **Picked the examples:** all squares where water was seen (1) and an equal number of squares seen dry
   (0), chosen randomly but balanced, at least 250 m away from any water.
8. **Quality-checked it:** 17 automatic checks per file (no duplicates, no future information, sensible
   values). Results: [docs/qa_report.md](docs/qa_report.md).

---

## 4. What one row and the main columns mean

**One row = one 30 m × 30 m square in Jeddah.**

| Column | Meaning |
|---|---|
| `label` | **The answer:** 1 = water seen, 0 = seen dry |
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

⚠️ **Rain columns: use with caution.** They come from one storm at a coarse scale (about 11 km), and by
chance the rain values are *lower* where the flooding was seen, which can mislead a model.

### 🚫 Do NOT use these columns for training
| Column | Why |
|---|---|
| `s2_water_confirmed` | Comes from photos taken *after* the flood; only for checking answer quality |
| `flood_fraction`, `label_source`, `label_obs_utc` | Describe the answer itself — using them is cheating |
| `hand_ref_is_sea` | Almost the same as "near the coast", so it leaks the answer |
| `sm_filled_from_neighbor` | Always 0 in this storm, so it carries no information |
| `sampling_stratum` | Only records how the square was sampled |
| `cell_id`, `event_id`, `x`, `y`, `lon`, `lat`, `block_id`, `t0_utc` | IDs and locations, not causes |

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

## 6. Final numbers

| File | Rows | Flooded (1) | Dry (0) |
|---|---|---|---|
| `jeddah_flood_ml` | **24,962** | 12,481 | 12,481 |
| `jeddah_flood_ml_elev_matched` | **21,995** | 12,481 | 9,514 |

The matched file has fewer dry squares because there are not enough dry places as low as the flooded ones.

---

## 7. Main limitations

1. **Radar cannot see inside the dense city.** About two-thirds of the area it observed was left out,
   mostly neighbourhoods. So the data shows where water **stayed in open land**, not which streets flooded.
2. **The picture is from 5 days later.** The storm was on 24 Nov; the radar image is from 29 Nov. We see
   water that **remained**, not the peak flood.
3. **About 40% of the "water" spots do not look like water in colour photos** — they may be wet soil or
   salt flats (*sabkha*).
4. **Most flooded spots are low and near the coast** (96% below 10 m, half below 0.6 m). The elevation
   map is only accurate to about 2–4 m, so tiny height differences near the sea are not reliable.
5. **Satellite rain is lower than measured:** our highest value is 81 mm, while the national weather
   centre reported 179 mm in 6 hours (as reported by FloodList).
6. **The elevation map includes buildings and trees**, not just bare ground.
7. **Only one storm**, and we chose a 50/50 mix of flooded and dry squares, so model scores are
   *relative*, not real-world probabilities.
8. ✅ **No rain fell between the storm and the radar image** (at most 1.3 mm), so the water seen came
   from this storm.

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
- Total download is about 450 MB; a full run takes roughly 20–40 minutes (an estimate, depending on internet and computer).

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
docs/checks/               ← maps used to check the answers (review only)
```
