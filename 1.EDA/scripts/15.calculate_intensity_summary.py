#!/usr/bin/env python
# coding: utf-8

# Section 7: intensities across channels, compartments, patients.
#
# Summarizes ALL intensity stat types per channel x compartment, grouped by
# patient/treatment, using utils_analysis's feature-name parsers. Uses
# organoid-level tables for the whole_organoid compartment and sc-level tables
# for nucleus/cytoplasm/cell compartments, for both 2D (3 projections) and 3D
# (different feature-name ordering handled by the two separate parsers).

# In[ ]:


import warnings

import pandas as pd

warnings.filterwarnings("ignore")


# In[ ]:


from notebook_init_utils import init_notebook

root_dir, in_notebook = init_notebook()
from eda_helper_utils.utils_analysis import (
    COMPARTMENT_LABELS,
    PROJECTION_FILE_PREFIX,
    PROJECTIONS,
    harmonize_metadata,
    list_patient_dirs,
    parse_feature_2d,
    parse_feature_3d,
)

results_dir = root_dir / "1.EDA" / "results" / "intensity"
results_dir.mkdir(parents=True, exist_ok=True)

# Curated subset of intensity measurements to carry through to plotting (boxplots
# read the row-level output below, not the full ~15-18 stat types per channel).
SELECTED_INTENSITY_STATS = ["MeanIntensity", "IntegratedIntensity", "MedianIntensity"]


# In[ ]:


def summarize_intensity(df, parse_fn, group_cols):
    """Return (agg, raw_values) for SELECTED_INTENSITY_STATS only.

    agg: mean/median/sd/n per group x compartment x channel x stat.
    raw_values: one row per observation x compartment x channel x stat, for boxplots.
    """
    intensity_cols = []
    parsed_map = {}
    for c in df.columns:
        if c.startswith("Metadata_"):
            continue
        obj, category, channel, measurement = parse_fn(c)
        if category == "Intensity" and measurement in SELECTED_INTENSITY_STATS:
            intensity_cols.append(c)
            parsed_map[c] = (obj, channel, measurement)
    if not intensity_cols:
        return pd.DataFrame(), pd.DataFrame()
    long = df.melt(
        id_vars=group_cols,
        value_vars=intensity_cols,
        var_name="feature",
        value_name="value",
    )
    long["value"] = pd.to_numeric(long["value"], errors="coerce")
    parsed = long["feature"].map(parsed_map)
    long["compartment_raw"] = parsed.apply(lambda x: x[0])
    long["compartment"] = (
        long["compartment_raw"].map(COMPARTMENT_LABELS).fillna(long["compartment_raw"])
    )
    long["channel"] = parsed.apply(lambda x: x[1])
    long["stat"] = parsed.apply(lambda x: x[2])
    agg = (
        long.groupby(group_cols + ["compartment", "channel", "stat"])["value"]
        .agg(mean="mean", median="median", sd="std", n="count")
        .reset_index()
    )
    raw_values = long.drop(columns=["feature", "compartment_raw"])
    return agg, raw_values


all_2d = []
raw_2d = []
patients_2d = list_patient_dirs(root_dir / "data" / "profiles_2D")
for projection in PROJECTIONS:
    prefix = PROJECTION_FILE_PREFIX[projection]
    for kind in ["organoid", "sc"]:
        for patient in patients_2d:
            f = (
                root_dir
                / "data"
                / "profiles_2D"
                / patient
                / "5.normalized"
                / f"{prefix}_{kind}.parquet"
            )
            if not f.exists():
                continue
            df = pd.read_parquet(f)
            df = harmonize_metadata(df, "2D", patient)
            agg, raw = summarize_intensity(
                df, parse_feature_2d, ["Metadata_patient", "Metadata_treatment"]
            )
            if agg.empty:
                continue
            agg["projection"] = projection
            agg["object_kind"] = kind
            all_2d.append(agg)
            raw["projection"] = projection
            raw["object_kind"] = kind
            raw_2d.append(raw)
        print(f"2D {projection} {kind}: done")

intensity_2d = pd.concat(all_2d, ignore_index=True)
intensity_2d.to_parquet(results_dir / "intensity_summary_2D.parquet", index=False)
print(
    f"Wrote {(results_dir / 'intensity_summary_2D.parquet').relative_to(root_dir)} ({len(intensity_2d)} rows)"
)

intensity_values_2d = pd.concat(raw_2d, ignore_index=True)
intensity_values_2d.to_parquet(results_dir / "intensity_values_2D.parquet", index=False)
print(
    f"Wrote {(results_dir / 'intensity_values_2D.parquet').relative_to(root_dir)} ({len(intensity_values_2d)} rows)"
)

all_3d = []
raw_3d = []
patients_3d = list_patient_dirs(root_dir / "data" / "profiles_3D")
for kind, fname in [("organoid", "organoid_norm.parquet"), ("sc", "sc_norm.parquet")]:
    for patient in patients_3d:
        f = (
            root_dir
            / "data"
            / "profiles_3D"
            / patient
            / "5.normalized_profiles"
            / fname
        )
        if not f.exists():
            continue
        df = pd.read_parquet(f)
        df = harmonize_metadata(df, "3D", patient)
        agg, raw = summarize_intensity(
            df, parse_feature_3d, ["Metadata_patient", "Metadata_treatment"]
        )
        if agg.empty:
            continue
        agg["object_kind"] = kind
        all_3d.append(agg)
        raw["object_kind"] = kind
        raw_3d.append(raw)
    print(f"3D {kind}: done")

intensity_3d = pd.concat(all_3d, ignore_index=True)
intensity_3d.to_parquet(results_dir / "intensity_summary_3D.parquet", index=False)
print(
    f"Wrote {(results_dir / 'intensity_summary_3D.parquet').relative_to(root_dir)} ({len(intensity_3d)} rows)"
)

intensity_values_3d = pd.concat(raw_3d, ignore_index=True)
intensity_values_3d.to_parquet(results_dir / "intensity_values_3D.parquet", index=False)
print(
    f"Wrote {(results_dir / 'intensity_values_3D.parquet').relative_to(root_dir)} ({len(intensity_values_3d)} rows)"
)
