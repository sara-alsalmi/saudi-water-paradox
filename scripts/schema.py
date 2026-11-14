"""Single source of truth for output columns -> used by 09_build_table, data_dictionary.csv and QA.

role:
  id             identifiers / coordinates / grouping (not predictors)
  feature        predictor known at or before t0
  trigger        event rainfall (allowed as trigger in the susceptibility design)
  label          target
  label_meta     describes the label; NOT a feature (leakage) — never train on it
"""
COLUMNS = [
    # name, role, unit, description, source, native_resolution, resampling_to_30m, time_reference
    ("cell_id", "id", "-", "Row-major index of the 30 m cell on the AOI reference grid (EPSG:32637, 60 m-snapped origin)", "derived", "-", "-", "-"),
    ("event_id", "id", "-", "Event identifier: JED_20171121 or JED_20221124. Key = cell_id + event_id. Use it to split train/test by event (train on one, test on the other)", "config", "-", "-", "-"),
    ("x", "id", "m", "Cell-centre easting, UTM 37N (EPSG:32637)", "derived", "30 m", "-", "-"),
    ("y", "id", "m", "Cell-centre northing, UTM 37N (EPSG:32637)", "derived", "30 m", "-", "-"),
    ("lon", "id", "deg", "Cell-centre longitude (WGS84)", "derived", "-", "-", "-"),
    ("lat", "id", "deg", "Cell-centre latitude (WGS84)", "derived", "-", "-", "-"),
    ("block_id", "id", "-", "1 km x 1 km spatial block (UTM km 'E_N' of lower-left corner) for spatial cross-validation", "derived", "1 km", "-", "-"),
    ("t0_utc", "id", "ISO-8601", "Prediction reference time = 00:00 UTC of the event day (2017-11-21, 2022-11-24)", "config", "-", "-", "-"),
    ("elevation_m", "feature", "m (EGM2008)", "Surface elevation. Copernicus DEM is a DSM: includes buildings/vegetation", "Copernicus DEM GLO-30 (2021 release)", "1 arcsec (~30 m)", "bilinear", "static (<= t0)"),
    ("slope_deg", "feature", "deg", "Slope (Horn 1981) of the 30 m DSM", "derived from Copernicus DEM", "30 m", "computed on 30 m grid", "static (<= t0)"),
    ("flowacc_km2", "feature", "km2", "D8 upstream contributing area incl. the cell (on catchment-extended, depression-filled DEM)", "derived from Copernicus DEM", "30 m", "computed on 30 m grid", "static (<= t0)"),
    ("twi", "feature", "-", "Topographic Wetness Index ln(a / tan b), a = acc_cells*30 m, tan b floored at 0.001", "derived from Copernicus DEM", "30 m", "computed on 30 m grid", "static (<= t0)"),
    ("hand_m", "feature", "m", "Height Above Nearest Drainage on the conditioned DSM. Drainage = channel cells (upstream area >= 1 km2) OR coastal outlet cells (sea treated as drainage outlet); see hand_ref_is_sea", "derived from Copernicus DEM", "30 m", "computed on 30 m grid", "static (<= t0)"),
    ("hand_ref_is_sea", "meta", "0/1", "1 if hand_m is measured to a coastal outlet cell (flow reaches the sea without crossing a >= 1 km2 channel), 0 if measured to a channel", "derived from Copernicus DEM", "30 m", "computed on 30 m grid", "static (<= t0)"),
    ("landcover", "feature", "class code", "ESA WorldCover class (10 tree,20 shrub,30 grass,40 crop,50 built-up,60 bare/sparse,80 water,90 wetland,95 mangrove)", "ESA WorldCover 2021 v200", "10 m", "mode (majority of 10 m pixels)", "2021 (<= t0)"),
    ("landcover_name", "feature", "-", "Name of landcover class (same information as landcover)", "ESA WorldCover 2021 v200", "10 m", "mode", "2021 (<= t0)"),
    ("landcover_year_gap", "meta", "years", "Event year minus WorldCover map year (2021): -4 (2017), +1 (2022). WorldCover 2021 is used for every event; large Jeddah demolitions/redevelopment 2021-2023 are NOT reflected", "derived", "-", "-", "-"),
    ("elev_bin", "feature", "-", "Fixed elevation bin: 0-10 (incl. <0), 10-25, 25-50, 50-100, >100 m; used for negative stratification", "derived", "30 m", "-", "static (<= t0)"),
    ("dist_to_coast_m", "feature", "m", "Euclidean distance from cell centre to nearest sea cell (DEM = 0 m, connected sea body)", "derived from Copernicus DEM", "30 m", "computed on 30 m grid", "static (<= t0)"),
    ("dist_to_ref_water_m", "feature", "m", "Distance from cell centre to nearest GFM reference (permanent/seasonal) water pixel", "Copernicus GFM reference water mask (2017-2021)", "20 m", "nearest to 10 m, distance at cell centre", "2017-2021 (<= t0)"),
    ("rain_antecedent_7d_mm", "feature", "mm", "IMERG precipitation total [t0 - 7 d, t0)", "GPM IMERG V07 Final half-hourly (see imerg_run)", "0.1 deg / 30 min", "nearest (cell containing centre)", "<= t0"),
    ("sm_surface_m3m3", "feature", "m3/m3", "SMAP L4 surface (0-5 cm) soil moisture, 3-h mean [t0 - 3 h, t0)", "NASA SMAP L4 SPL4SMGP v008", "9 km EASE-Grid 2.0", "nearest EASE cell", "<= t0"),
    ("sm_rootzone_m3m3", "feature", "m3/m3", "SMAP L4 root-zone (0-100 cm) soil moisture, same window", "NASA SMAP L4 SPL4SMGP v008", "9 km EASE-Grid 2.0", "nearest EASE cell", "<= t0"),
    ("sm_filled_from_neighbor", "meta", "0/1", "1 if SMAP had no value in the cell's EASE cell (coast) and the nearest valid EASE cell was used", "derived", "-", "-", "<= t0"),
    ("rain_event_24h_mm", "trigger", "mm", "IMERG total [t0, t0 + 24 h) (event trigger)", "GPM IMERG V07 Final half-hourly (see imerg_run)", "0.1 deg / 30 min", "nearest (cell containing centre)", "event window (trigger, after t0)"),
    ("rain_event_max30_mmhr", "trigger", "mm/hr", "Max half-hourly IMERG rate in the event window", "GPM IMERG V07 Final half-hourly (see imerg_run)", "0.1 deg / 30 min", "nearest", "event window (trigger, after t0)"),
    ("rain_upstream_event_mm", "trigger", "mm", "Mean event rainfall over the cell's full upstream catchment (flow-weighted accumulation / cell count)", "IMERG + Copernicus DEM", "0.1 deg routed on 30 m D8", "nearest to 30 m then routed", "event window (trigger, after t0)"),
    ("imerg_run", "meta", "-", "IMERG run used for this event's rain columns: Final (research quality, gauge-calibrated) or Late (near-real-time, no monthly gauge calibration). Both current events use Final; the column exists so that a future Late-run event (e.g. the rejected JED_20251209) is visible. Values are not strictly comparable across runs", "GPM IMERG V07", "-", "-", "-"),
    ("label", "label", "0/1", "1 = NEW flood on the event's GFM label scene (new-flood fraction >= 0.5; new = flooded on the label scene and not flooded on the baseline scene); 0 = observed with no water on the label scene (all 10 m subpixels observed, not excluded, not reference water, not baseline water, >= 250 m from any new-flood or baseline-water pixel). 0 does NOT mean safe", "Copernicus GFM ensemble_flood_extent (version: see gfm_version)", "20 m", "nearest to 10 m sub-grid, fraction per 30 m cell", "label scene 4-5 days after the event (after event)"),
    ("label_source", "label_meta", "-", "Label product, version and label-scene date", "Copernicus GFM", "-", "-", "-"),
    ("gfm_version", "label_meta", "-", "GFM processing version of the label and baseline scenes. Both current events: v3.1.0 (2017 and 2022 are archive reprocessing with the same version). Newer scenes (2025+) are v4.x with a different algorithm and exclusion mask", "Copernicus GFM", "-", "-", "-"),
    ("label_obs_utc", "label_meta", "ISO-8601", "Sentinel-1 acquisition time of the labelling scene — NOT a feature", "Copernicus GFM", "-", "-", "label scene (after event)"),
    ("flood_fraction", "label_meta", "0-1", "Fraction of the cell's 10 m subpixels flagged as NEW flood by GFM — derived from the label, NOT a feature", "Copernicus GFM", "20 m", "nearest to 10 m, mean per 30 m cell", "label scene (after event)"),
    ("s2_water_confirmed", "label_meta", "0/1/NaN", "Label-quality flag, NOT a feature — do not use for training. 1 if (NDWI post - NDWI pre) > 0.1 or MNDWI post > 0 (cell means over valid 10 m pixels), else 0; NaN if < 5 of 9 subpixels valid on either side. Pre/post Sentinel-2 dates per event: 2017 pre 18 Nov (fallback 13 Nov), post 23 Nov; 2022 pre 22 Nov, post 27 Nov. 2017: mostly NaN in the north of the city (post-event tiles cloudy there)", "Sentinel-2 L2A (Earth Search)", "10/20 m", "nearest/average to 10 m, mean per 30 m cell", "after event"),
    ("sampling_stratum", "label_meta", "-", "Stratum used for negative sampling: 'landcover|elev_bin' (stratified file) or 'elevq<k>[lo,hi)' positive-elevation quantile bin (elev_matched file)", "derived", "-", "-", "-"),
]

RAIN_CAUTION = ("USE WITH CAUTION — only 24 distinct values per event (~11 km IMERG cells) and a spurious "
                "spatial correlation with the label in BOTH events (IMERG rain is lower at the coastal positives)")
CAUTION = {c: RAIN_CAUTION for c in ("rain_antecedent_7d_mm", "rain_event_24h_mm", "rain_event_max30_mmhr",
                                     "rain_upstream_event_mm")}
CAUTION["s2_water_confirmed"] = "NOT A FEATURE — do not use for training (post-event label-quality metadata)"
for _c in ("label_source", "label_obs_utc", "flood_fraction", "sampling_stratum"):
    CAUTION[_c] = "NOT A FEATURE — label metadata"
CAUTION["rain_upstream_event_mm"] = (RAIN_CAUTION + "; nearly identical to rain_event_24h_mm "
                                     "(equal in ~96 % of 2022 rows and ~98 % of 2017 rows: most catchments lie inside "
                                     "one IMERG cell) — do not treat it as an independent variable")
CAUTION["sm_filled_from_neighbor"] = ("constant 0 in both events (no coastal SMAP fill was needed) — no information")
CAUTION["landcover_year_gap"] = ("metadata (constant per event) — do not use for training; WorldCover 2021 does not show "
                                 "2021-2023 demolitions, nor 2017 conditions")
CAUTION["imerg_run"] = "metadata (constant per event) — do not use for training"
CAUTION["gfm_version"] = "NOT A FEATURE — label metadata"
CAUTION["hand_ref_is_sea"] = ("metadata flag; strongly associated with the label (coastal proxy) — "
                              "do not use for training")

NAMES = [c[0] for c in COLUMNS]
ROLE = {c[0]: c[1] for c in COLUMNS}
MODEL_INPUTS = [c[0] for c in COLUMNS if c[1] in ("feature", "trigger")]
NEVER_FEATURES = [c[0] for c in COLUMNS if c[1] in ("label", "label_meta", "meta", "id")]


def write_dictionary(path) -> None:
    """data_dictionary.csv from COLUMNS (no table rebuild needed)."""
    import csv
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["column", "role", "use_as_model_input", "caution", "unit", "description", "source",
                    "native_resolution", "resampling_to_30m_grid", "time_reference"])
        for c in COLUMNS:
            use = "yes" if c[1] in ("feature", "trigger") else "NO"
            w.writerow([c[0], c[1], use, CAUTION.get(c[0], ""), c[2], c[3], c[4], c[5], c[6], c[7]])


if __name__ == "__main__":
    from pathlib import Path
    write_dictionary(Path(__file__).resolve().parents[1] / "data_dictionary.csv")
    print("data_dictionary.csv written")

WORLDCOVER = {10: "tree_cover", 20: "shrubland", 30: "grassland", 40: "cropland", 50: "built_up",
              60: "bare_sparse", 70: "snow_ice", 80: "permanent_water", 90: "herbaceous_wetland",
              95: "mangroves", 100: "moss_lichen"}
ELEV_BINS = [-float("inf"), 10, 25, 50, 100, float("inf")]
ELEV_LABELS = ["0-10", "10-25", "25-50", "50-100", ">100"]
