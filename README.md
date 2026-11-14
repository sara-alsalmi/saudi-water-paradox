# The Saudi Water Paradox — Jeddah flood dataset

How can a water-scarce country experience devastating floods? This repository contains a dataset for
studying which places in Jeddah are susceptible to flash-flood water, and why. It combines satellite flood
maps with terrain, land cover, rainfall and soil moisture for two storms: **21 November 2017** and
**24 November 2022**.

## Dataset at a glance

- **One row is one 30 m × 30 m cell in one storm.** The same cell can appear once per storm; `event_id`
  says which storm (`JED_20171121` or `JED_20221124`).
- **`label = 1`:** the radar saw new flood water in the cell a few days after the storm.
- **`label = 0`:** the radar observed the cell and saw no water. This means "observed dry on that date",
  not "safe".

| File | Storm | Rows | Flooded (1) | Dry (0) |
|---|---|---|---|---|
| `jeddah_flood_ml` | 2017 | 15,092 | 7,546 | 7,546 |
| | 2022 | 24,962 | 12,481 | 12,481 |
| | **both** | **40,054** | 20,027 | 20,027 |
| `jeddah_flood_ml_elev_matched` | 2017 | 13,849 | 7,546 | 6,303 |
| | 2022 | 21,995 | 12,481 | 9,514 |
| | **both** | **35,844** | 20,027 | 15,817 |

## Data sources

| Source | What it adds |
|---|---|
| Copernicus GFM (Sentinel-1 radar flood maps) | The label: where flood water was seen. Radar sees through clouds. |
| Copernicus DEM GLO-30 | Terrain: elevation, slope, flow paths and the wadis that drain into Jeddah. |
| ESA WorldCover 2021 | Land cover: built-up, bare ground, shrubland, cropland, water. |
| NASA GPM IMERG | Rainfall: on the storm day, the heaviest half hour, and the week before. |
| NASA SMAP | Soil moisture just before the storm. |

Sentinel-2 optical images were used only to verify the labels, never as model inputs.

## How it was built

1. Defined the study area around Jeddah (39.05–39.40° E, 21.25–21.85° N) on a 30 m grid.
2. Selected the storms: every Jeddah flood since 2009 was checked, and only storms with a clean radar
   image a few days after the rain were kept (2017 and 2022).
3. Built the labels from GFM, counting only *new* water: water already present before the storm is
   excluded. The labels were checked against Sentinel-2 images ([docs/checks/](docs/checks/)).
4. Downloaded the other sources, clipped to the area and listed in
   [data/raw/manifest.csv](data/raw/manifest.csv) with links and checksums.
5. Cleaned and aligned everything: fill codes became missing values, units were made consistent, all
   layers share the same 30 m grid, and all times are in UTC.
6. Traced the upstream wadis that drain into Jeddah (about 7,200 km²) for the terrain and rainfall features.
7. Used only information from before the storm began, except the storm's own rain.
8. Sampled, per storm, all flooded cells and an equal number of dry cells at least 250 m from any water.
9. Ran 32 automatic quality checks per file, for the whole file and for each storm
   ([docs/qa_report.md](docs/qa_report.md)).

## Files

| File | Description |
|---|---|
| `data/processed/jeddah_flood_ml.csv` / `.parquet` | Main table for modelling. |
| `data/processed/jeddah_flood_ml_elev_matched.csv` / `.parquet` | Check table: dry cells are sampled at the same elevations as flooded cells. |
| [data_dictionary.csv](data_dictionary.csv) | Description of every column, including `use_as_model_input`. |
| [notebooks/01_data_preparation.ipynb](notebooks/01_data_preparation.ipynb) | Rebuilds the data. |
| [notebooks/02_data_sanity_check.ipynb](notebooks/02_data_sanity_check.ipynb) | Quick check of the final data. |

Both tables have the same columns. CSV opens in Excel; Parquet is faster and keeps column types.

**Why two tables?** In the main table most flooded cells are low and most dry cells are higher, so a model
can score well by learning "low = flooded". In the elevation-matched table that shortcut no longer works.
It has fewer dry cells because there are not enough dry places as low as the flooded ones.

## Using the data

- Use only the columns with `use_as_model_input = yes` in the data dictionary. The others are IDs,
  label metadata or quality flags.
- Split by place, using `block_id` (1 km blocks), so that test cells are in areas the model has not seen.
- Random splits inflate scores, because neighbouring cells are nearly identical.
- Always compare against an elevation-only baseline, and report results on both tables.
- Use the rain columns with caution: they are coarse (about 11 km), and in both storms rain is by chance
  lower where flooding was seen.

More guidance: [docs/handover.md](docs/handover.md).

## Quick data check

From [notebooks/02_data_sanity_check.ipynb](notebooks/02_data_sanity_check.ipynb):

- The data is clean and balanced, and no leakage was found.
- The labels are consistent: 5,166 cells flooded in both storms, against about 35 expected by chance.
- Elevation alone reaches ROC-AUC 0.961; a Random Forest on place features reaches 0.978.
- On the elevation-matched table the model reaches 0.867, against 0.586 for elevation alone.

## Limitations

- The radar cannot see inside the dense city: about two-thirds of the observed area is excluded, so the
  labels show where water stayed in open land, not which streets flooded.
- The radar images are from 4–5 days after each storm (25 Nov 2017 and 29 Nov 2022), so they show
  remaining water, not the peak flood.
- Not all flood pixels look like water in optical images: in 2022 about 40% do not. In 2017 the north of
  the city could not be checked because of clouds.
- Most flooded cells are low and near the coast (96% below 10 m), and mostly the same places in both
  storms. The elevation data is accurate to about 2–4 m and includes buildings and trees.
- Satellite rain is lower than measured: the highest 2022 value is 81 mm, while 179 mm in 6 hours was
  reported by the national weather centre (via FloodList).
- Land cover is from 2021 for both storms, and there are only two storms. With a 50/50 sample, model
  scores are relative, not real-world probabilities.

No rain fell between each storm and its radar image (at most 1.3 mm). Details:
[docs/data_sources.md](docs/data_sources.md).

## Reproducing the dataset

You need Python 3.13, internet access, a free NASA Earthdata account and about 1 GB of free space.

**1. Set up the environment** (once, in a terminal in this folder):
```
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
```
Use `requirements.lock.txt` instead for exactly the same versions. To run the notebooks, register the
environment as a kernel:
```
.venv\Scripts\python -m ipykernel install --sys-prefix --name saudi-water-paradox --display-name "Python (saudi-water-paradox .venv)"
```

**2. Log in to NASA Earthdata** (once; you type your own password and nothing is saved in the project):
```
.venv\Scripts\python -c "import earthaccess; earthaccess.login(strategy='interactive', persist=True)"
```
Then, on the Earthdata website, go to Profile → Applications → Authorized Apps and approve
**NASA GESDISC DATA ARCHIVE**.

**3. Run all steps.** Open [notebooks/01_data_preparation.ipynb](notebooks/01_data_preparation.ipynb)
with the kernel **Python (saudi-water-paradox .venv)** and run all cells, or run:
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

All settings are in [config.yaml](config.yaml). Re-running is safe (files already downloaded are skipped)
and gives identical results. The total download is about 690 MB, and a full run takes roughly 30–50
minutes, depending on the connection and computer.

## Project structure

```
README.md                  this file
config.yaml                all settings
requirements.txt           required libraries
data_dictionary.csv        description of every column
data/raw/                  original downloads (re-created by the scripts) and manifest.csv
data/interim/              intermediate files (re-created automatically)
data/processed/            final tables
scripts/                   pipeline steps 01–10
notebooks/
  01_data_preparation.ipynb    rebuilds the data
  02_data_sanity_check.ipynb   quick check of the final data
docs/data_sources.md       sources, processing, limitations, rejected storms
docs/qa_report.md          quality-check report
docs/handover.md           guidance for modelling
docs/checks/<event_id>/    label-check maps, one folder per storm
```

## Data sources and credits

All sources are free and open. Licences, citations and required attributions are listed in
[docs/data_sources.md](docs/data_sources.md).
