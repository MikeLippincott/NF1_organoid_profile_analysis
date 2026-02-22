#!/usr/bin/env python
# coding: utf-8

# In[1]:


import itertools
import pathlib
import time

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

# In[2]:


root_dir = pathlib.Path.cwd().parent.parent


# Output directory for correlation results
results_dir = root_dir / "1.EDA" / "results" / "correlation"
results_dir.mkdir(parents=True, exist_ok=True)
print(f"Results directory: {results_dir}")


# In[3]:


# Define data paths
# 2D profiles
path_2d_organoid = (
    root_dir
    / "data"
    / "2D_profiles"
    / "all_patient_profiles"
    / "max_projection"
    / "organoid_agg_profiles.parquet"
)
path_2d_sc = (
    root_dir
    / "data"
    / "2D_profiles"
    / "all_patient_profiles"
    / "max_projection"
    / "sc_agg_profiles.parquet"
)

# 3D profiles
path_3d_organoid = (
    root_dir
    / "data"
    / "3D_profiles"
    / "all_patient_profiles"
    / "organoid_agg_profiles.parquet"
)
path_3d_sc = (
    root_dir
    / "data"
    / "3D_profiles"
    / "all_patient_profiles"
    / "sc_agg_profiles.parquet"
)

# Load all data
df_2d_organoid = pd.read_parquet(path_2d_organoid)
df_2d_sc = pd.read_parquet(path_2d_sc)
df_3d_organoid = pd.read_parquet(path_3d_organoid)
df_3d_sc = pd.read_parquet(path_3d_sc)


# In[4]:


# Get patients present in both 2D and 3D organoid datasets
patients_organoid = sorted(
    set(df_2d_organoid["Metadata_patient"].unique())
    & set(df_3d_organoid["Metadata_patient"].unique())
)
patients_sc = sorted(
    set(df_2d_sc["Metadata_patient"].unique())
    & set(df_3d_sc["Metadata_patient"].unique())
)


# In[5]:


# Z-score standardize columns, set constant columns (std=0) by setting them to NaN
def standardize(arr):
    mean = np.nanmean(arr, axis=0)
    std = np.nanstd(arr, axis=0, ddof=1)
    std[std == 0] = np.nan
    return (arr - mean) / std


# ### Organoid Correlations per patient
#

# In[6]:


for patient_id in patients_organoid:
    # Match data in 2D and 3D ensuring wells are matched for a tumor
    tumors_2d = set(
        df_2d_organoid[df_2d_organoid["Metadata_patient"] == patient_id][
            "Metadata_tumor"
        ].unique()
    )
    tumors_3d = set(
        df_3d_organoid[df_3d_organoid["Metadata_patient"] == patient_id][
            "Metadata_tumor"
        ].unique()
    )
    shared_tumors = sorted(tumors_2d & tumors_3d)

    df_2d_patient = df_2d_organoid[
        (df_2d_organoid["Metadata_patient"] == patient_id)
        & (df_2d_organoid["Metadata_tumor"].isin(shared_tumors))
    ].copy()
    df_3d_patient = df_3d_organoid[
        (df_3d_organoid["Metadata_patient"] == patient_id)
        & (df_3d_organoid["Metadata_tumor"].isin(shared_tumors))
    ].copy()

    wells_2d = set(df_2d_patient["Metadata_Well"].unique())
    wells_3d = set(df_3d_patient["Metadata_Well"].unique())
    matched_wells = sorted(wells_2d & wells_3d)

    if len(matched_wells) == 0:
        raise ValueError(
            f"Patient {patient_id}: No matched wells ({len(matched_wells)})."
        )

    # Dataframe alignment
    df_2d_matched = (
        df_2d_patient[df_2d_patient["Metadata_Well"].isin(matched_wells)]
        .sort_values("Metadata_Well")
        .reset_index(drop=True)
    )
    df_3d_matched = (
        df_3d_patient[df_3d_patient["Metadata_Well"].isin(matched_wells)]
        .sort_values("Metadata_Well")
        .reset_index(drop=True)
    )

    assert list(df_2d_matched["Metadata_Well"]) == list(
        df_3d_matched["Metadata_Well"]
    ), "Wells are not aligned after sorting"

    # Extract feature columns
    features_2d = [
        col for col in df_2d_matched.columns if not col.startswith("Metadata_")
    ]
    features_3d = [
        col for col in df_3d_matched.columns if not col.startswith("Metadata_")
    ]

    # Z-score standardize feature arrays for Pearson correlation
    arr_2d = df_2d_matched[features_2d].values
    arr_3d = df_3d_matched[features_3d].values

    arr_2d_z = standardize(arr_2d)
    arr_3d_z = standardize(arr_3d)
    n_wells = arr_2d_z.shape[0]

    # Compute correlation via matrix multiplication
    corr_matrix = np.dot(arr_2d_z.T, arr_3d_z) / (n_wells - 1)

    # Convert correlation matrix to long-form DataFrame and save
    corr_df = pd.DataFrame(corr_matrix, index=features_2d, columns=features_3d)
    corr_long = (
        corr_df.reset_index()
        .melt(id_vars="index", var_name="feature_3d", value_name="pearson_r")
        .rename(columns={"index": "feature_2d"})
    )
    output_path = results_dir / f"{patient_id}_organoid_agg_correlation.parquet"
    corr_long.to_parquet(output_path, index=False)


# ### SC Correlations Per patient

# In[7]:


for patient_id in patients_sc:
    # Match data in 2D and 3D ensuring wells are matched for a tumor
    tumors_2d = set(
        df_2d_sc[df_2d_sc["Metadata_patient"] == patient_id]["Metadata_tumor"].unique()
    )
    tumors_3d = set(
        df_3d_sc[df_3d_sc["Metadata_patient"] == patient_id]["Metadata_tumor"].unique()
    )
    shared_tumors = sorted(tumors_2d & tumors_3d)

    df_2d_patient = df_2d_sc[
        (df_2d_sc["Metadata_patient"] == patient_id)
        & (df_2d_sc["Metadata_tumor"].isin(shared_tumors))
    ].copy()
    df_3d_patient = df_3d_sc[
        (df_3d_sc["Metadata_patient"] == patient_id)
        & (df_3d_sc["Metadata_tumor"].isin(shared_tumors))
    ].copy()

    wells_2d = set(df_2d_patient["Metadata_Well"].unique())
    wells_3d = set(df_3d_patient["Metadata_Well"].unique())
    matched_wells = sorted(wells_2d & wells_3d)

    if len(matched_wells) == 0:
        raise ValueError(
            f"Patient {patient_id}: No matched wells ({len(matched_wells)})."
        )

    # Dataframe alignment
    df_2d_matched = (
        df_2d_patient[df_2d_patient["Metadata_Well"].isin(matched_wells)]
        .sort_values("Metadata_Well")
        .reset_index(drop=True)
    )
    df_3d_matched = (
        df_3d_patient[df_3d_patient["Metadata_Well"].isin(matched_wells)]
        .sort_values("Metadata_Well")
        .reset_index(drop=True)
    )

    assert list(df_2d_matched["Metadata_Well"]) == list(
        df_3d_matched["Metadata_Well"]
    ), "Wells are not aligned after sorting"

    # Extract feature columns
    features_2d = [
        col for col in df_2d_matched.columns if not col.startswith("Metadata_")
    ]
    features_3d = [
        col for col in df_3d_matched.columns if not col.startswith("Metadata_")
    ]

    # Z-score standardize feature arrays for Pearson correlation
    arr_2d = df_2d_matched[features_2d].values
    arr_3d = df_3d_matched[features_3d].values

    arr_2d_z = standardize(arr_2d)
    arr_3d_z = standardize(arr_3d)
    n_wells = arr_2d_z.shape[0]

    # Compute correlation via matri multiplication
    corr_matrix = np.dot(arr_2d_z.T, arr_3d_z) / (n_wells - 1)

    # Convert correlation matrix to long-form DataFrame and save
    start_time = time.time()
    corr_df = pd.DataFrame(corr_matrix, index=features_2d, columns=features_3d)
    corr_long = (
        corr_df.reset_index()
        .melt(id_vars="index", var_name="feature_3d", value_name="pearson_r")
        .rename(columns={"index": "feature_2d"})
    )
    output_path = results_dir / f"{patient_id}_sc_agg_correlation.parquet"
    corr_long.to_parquet(output_path, index=False)
