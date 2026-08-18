#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import numpy as np
import pandas as pd
import pyarrow as pq
from notebook_init_utils import bandicoot_check, init_notebook

root_dir, in_notebook = init_notebook()

if in_notebook:
    import tqdm.notebook as tqdm
else:
    import tqdm


# In[2]:


normalized_3D_profiles = pathlib.Path(
    f"{root_dir}/data/profiles_3D/all_patients/0.normalized_profiles/"
).resolve()
normalized_profiles = list(normalized_3D_profiles.glob("*.parquet"))
normalized_2D_profiles = pathlib.Path(
    f"{root_dir}/data/profiles_2D/all_patients/"
).resolve()
normalized_2D_profile_paths = [
    x
    for y in normalized_2D_profiles.glob("*")
    if y.is_dir()
    for x in y.glob("*")
    if not any(s in x.name for s in ["fs", "agg", "consensus"])
]
normalized_profiles = normalized_profiles + normalized_2D_profile_paths


# In[3]:


results_dir = pathlib.Path(f"{root_dir}/1.EDA/results/cell_counts/").resolve()
results_dir.mkdir(parents=True, exist_ok=True)


# In[4]:


# --- Pass 1: build counts_df for every profile type, split by 2D/3D ---
dict_of_count_dfs = {}
fov_lookup_frames = []

for path in tqdm.tqdm(normalized_profiles):
    profile_type = f"{path.stem}_{path.parent.parent.parent.name.split('_')[1]}"
    df = pd.read_parquet(path)
    metadata_cols = [x for x in df.columns if "Metadata" in x]
    df = df[metadata_cols]

    is_2d = "2D" in str(path)
    if is_2d:
        profile_type = f"{profile_type}_{path.parent.name}"
        df = df.rename(
            columns={
                "Metadata_patient_tumor": "Metadata_Biology_PatientTumor",
                "Metadata_treatment": "Metadata_Experiment_Treatment",
                "Metadata_dose": "Metadata_Experiment_Dose",
                "Metadata_Well": "Metadata_Experiment_Well",
            }
        )

    group_cols = [
        "Metadata_Biology_PatientTumor",
        "Metadata_Experiment_Treatment",
        "Metadata_Experiment_Dose",
    ]

    has_wellfov = "Metadata_Experiment_WellFOV" in df.columns

    agg_dict = {
        "Metadata_n_wells": pd.NamedAgg(
            column="Metadata_Experiment_Well", aggfunc="nunique"
        ),
    }
    if has_wellfov:
        agg_dict["Metadata_n_cells"] = pd.NamedAgg(
            column="Metadata_Experiment_WellFOV", aggfunc="count"
        )
        agg_dict["Metadata_n_fovs"] = pd.NamedAgg(
            column="Metadata_Experiment_WellFOV", aggfunc="nunique"
        )
    else:
        agg_dict["Metadata_n_cells"] = pd.NamedAgg(
            column="Metadata_Experiment_Well", aggfunc="count"
        )

    counts_df = df.groupby(group_cols).agg(**agg_dict).reset_index()

    if not has_wellfov:
        counts_df["Metadata_n_fovs"] = pd.NA

    counts_df["Metadata_profile_name"] = profile_type
    dict_of_count_dfs[profile_type] = {"df": counts_df, "is_2d": is_2d}

    # build the well-level FOV lookup from 3D data only
    if has_wellfov:
        well_fov_df = (
            df.groupby(group_cols + ["Metadata_Experiment_Well"])[
                "Metadata_Experiment_WellFOV"
            ]
            .nunique()
            .reset_index(name="Metadata_n_fovs_lookup")
        )
        fov_lookup_frames.append(well_fov_df)

# combine all 3D-derived fov lookups
fov_lookup_raw = pd.concat(fov_lookup_frames)

# dedupe: keep one row per patient/treatment/dose/well combo
# (multiple 3D profile types should all report the same FOV count per well;
# if they don't exactly agree, this keeps the first one seen — adjust if you'd rather average)
fov_lookup_by_well = fov_lookup_raw.drop_duplicates(
    subset=[
        "Metadata_Biology_PatientTumor",
        "Metadata_Experiment_Treatment",
        "Metadata_Experiment_Dose",
        "Metadata_Experiment_Well",
    ]
)

# now sum FOVs across wells within each patient/treatment/dose group, without double-counting
fov_lookup = (
    fov_lookup_by_well.groupby(
        [
            "Metadata_Biology_PatientTumor",
            "Metadata_Experiment_Treatment",
            "Metadata_Experiment_Dose",
        ]
    )["Metadata_n_fovs_lookup"]
    .sum()
    .reset_index()
)
# --- Pass 2: fill missing n_fovs in 2D profiles using the 3D-derived lookup ---
for profile_type, entry in dict_of_count_dfs.items():
    if entry["is_2d"]:
        df = (
            entry["df"]
            .drop(columns=["Metadata_n_fovs"])
            .merge(
                fov_lookup,
                on=[
                    "Metadata_Biology_PatientTumor",
                    "Metadata_Experiment_Treatment",
                    "Metadata_Experiment_Dose",
                ],
                how="left",
            )
            .rename(columns={"Metadata_n_fovs_lookup": "Metadata_n_fovs"})
        )
        dict_of_count_dfs[profile_type]["df"] = df

    dict_of_count_dfs[profile_type]["df"].to_parquet(
        results_dir / f"{profile_type}_counts.parquet", index=False
    )

# --- Pass 3: normalize the n_cells column by the number of FOVs, to get a per-FOV cell count ---
# annotate unique treatments
for profile_type, entry in dict_of_count_dfs.items():
    df = entry["df"]
    df["Metadata_n_cells_norm_by_well_fov"] = (
        df["Metadata_n_cells"] / df["Metadata_n_fovs"]
    )
    df["Metadata_unique_treatment"] = (
        df["Metadata_Experiment_Treatment"]
        + "_"
        + df["Metadata_Experiment_Dose"].astype(str)
    )
    dict_of_count_dfs[profile_type]["df"] = df
    dict_of_count_dfs[profile_type]["df"].to_parquet(
        results_dir / f"{profile_type}_counts.parquet", index=False
    )
