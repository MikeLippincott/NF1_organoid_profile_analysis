#!/usr/bin/env python
# coding: utf-8

# In[1]:


import itertools
import pathlib

import numpy as np
import pandas as pd
from notebook_init_utils import init_notebook
from scipy.stats import pearsonr

root_dir, in_notebook = init_notebook()


# In[2]:


# Define input paths
# 2D profiles
path_2d_organoid = pathlib.Path(
    f"{root_dir}/data/2D_profiles/all_patient_profiles/max_projection/organoid_agg_profiles.parquet"
).resolve(strict=True)
path_2d_sc = pathlib.Path(
    f"{root_dir}/data/2D_profiles/all_patient_profiles/max_projection/sc_agg_profiles.parquet"
).resolve(strict=True)

# 3D profiles
path_3d_organoid = pathlib.Path(
    f"{root_dir}/data/3D_profiles/all_patient_profiles/organoid_agg_profiles.parquet"
).resolve(strict=True)
path_3d_sc = pathlib.Path(
    f"{root_dir}/data/3D_profiles/all_patient_profiles/sc_agg_profiles.parquet"
).resolve(strict=True)

# Define output paths
results_dir = pathlib.Path(
    f"{root_dir}/2.2d_vs_3d_analysis/results/correlation"
).resolve()
results_dir.mkdir(parents=True, exist_ok=True)

# Load all data
df_2d_organoid = pd.read_parquet(path_2d_organoid)
df_2d_sc = pd.read_parquet(path_2d_sc)
df_3d_organoid = pd.read_parquet(path_3d_organoid)
df_3d_sc = pd.read_parquet(path_3d_sc)


# In[3]:


# Get patients present in both 2D and 3D organoid datasets
patients_organoid = sorted(
    set(df_2d_organoid["Metadata_patient"].unique())
    & set(df_3d_organoid["Metadata_patient"].unique())
)
patients_sc = sorted(
    set(df_2d_sc["Metadata_patient"].unique())
    & set(df_3d_sc["Metadata_patient"].unique())
)


# ### Generate Correlation Files

# In[4]:


def standardize(arr: np.ndarray) -> np.ndarray:
    """
    Perform standard scalar normalization on the input array.
    This is NA aware!

    Parameters
    ----------
    arr : np.ndarray
        The input array to be normalized in the following dimensionality: (nxm)

    Returns
    ----------
    np.ndarray
        The standard scalar normalized array of the same dimensionality (nxm)
    """
    mean = np.nanmean(arr, axis=0)
    std = np.nanstd(arr, axis=0, ddof=0)
    std = np.where(std == 0, np.nan, std)
    return (arr - mean) / std


# In[5]:


def compute_and_save_correlation(
    df_2d: pd.DataFrame,
    df_3d: pd.DataFrame,
    patient_id: str,
    results_dir: pathlib.Path,
    profile_type: str,
) -> None:
    """
    Compute pairwise Pearson correlations between 2D and 3D features
    for a given patient, then save the results as a long-form parquet file.

    Parameters
    ----------
    df_2d : pd.DataFrame
        2D profile data containing Metadata_patient, Metadata_tumor,
        Metadata_Well, and feature columns.
    df_3d : pd.DataFrame
        3D profile data with the same metadata structure.
    patient_id : str
        The patient identifier to filter on.
    results_dir : pathlib.Path
        Directory to save the output parquet file.
    profile_type : str
        Label for the output file (e.g., "organoid_agg", "sc_agg").
    """
    # Match tumors present in both 2D and 3D
    tumors_2d = set(
        df_2d[df_2d["Metadata_patient"] == patient_id]["Metadata_tumor"].unique()
    )
    tumors_3d = set(
        df_3d[df_3d["Metadata_patient"] == patient_id]["Metadata_tumor"].unique()
    )
    shared_tumors = sorted(tumors_2d & tumors_3d)

    df_2d_patient = df_2d[
        (df_2d["Metadata_patient"] == patient_id)
        & (df_2d["Metadata_tumor"].isin(shared_tumors))
    ].copy()
    df_3d_patient = df_3d[
        (df_3d["Metadata_patient"] == patient_id)
        & (df_3d["Metadata_tumor"].isin(shared_tumors))
    ].copy()

    # Match wells present in both
    wells_2d = set(df_2d_patient["Metadata_Well"].unique())
    wells_3d = set(df_3d_patient["Metadata_Well"].unique())
    matched_wells = sorted(wells_2d & wells_3d)

    if len(matched_wells) == 0:
        raise ValueError(
            f"Patient {patient_id}: No matched wells ({len(matched_wells)})."
        )

    # Align dataframes by well
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

    # Extract and standardize feature arrays
    features_2d = [
        col for col in df_2d_matched.columns if not col.startswith("Metadata_")
    ]
    features_3d = [
        col for col in df_3d_matched.columns if not col.startswith("Metadata_")
    ]

    # Drop features that are all NaN
    df_2d_clean = df_2d_matched[features_2d].dropna(axis=1, how="all")
    df_3d_clean = df_3d_matched[features_3d].dropna(axis=1, how="all")

    features_2d = list(df_2d_clean.columns)
    features_3d = list(df_3d_clean.columns)

    arr_2d_z = standardize(df_2d_matched[features_2d].values)
    arr_3d_z = standardize(df_3d_matched[features_3d].values)
    n_wells = arr_2d_z.shape[0]

    # Compute correlation via matrix multiplication (NaN-aware)
    df_2d_std = pd.DataFrame(arr_2d_z, columns=features_2d)
    df_3d_std = pd.DataFrame(arr_3d_z, columns=features_3d)
    combined = pd.concat([df_2d_std, df_3d_std], axis=1)
    full_corr = combined.corr(method="pearson")
    corr_matrix = full_corr.loc[features_2d, features_3d].values

    # Convert to long-form and save
    corr_df = pd.DataFrame(corr_matrix, index=features_2d, columns=features_3d)
    corr_long = (
        corr_df.reset_index()
        .melt(id_vars="index", var_name="feature_3d", value_name="pearson_r")
        .rename(columns={"index": "feature_2d"})
    )
    output_path = results_dir / f"{patient_id}_{profile_type}_correlation.parquet"
    corr_long.to_parquet(output_path, index=False)


# ### Organoid Correlations per patient
#

# In[6]:


for patient_id in patients_organoid:
    compute_and_save_correlation(
        df_2d_organoid, df_3d_organoid, patient_id, results_dir, "organoid_agg"
    )


# ### SC Correlations Per patient

# In[7]:


for patient_id in patients_sc:
    compute_and_save_correlation(df_2d_sc, df_3d_sc, patient_id, results_dir, "sc_agg")


# ### All patient correlations

# In[8]:


def compute_and_save_correlation_all(
    df_2d: pd.DataFrame,
    df_3d: pd.DataFrame,
    results_dir: pathlib.Path,
    profile_type: str,
) -> None:
    """
    Compute pairwise Pearson correlations between 2D and 3D features
    across all patients combined, then save as a long-form parquet file.

    Parameters
    ----------
    df_2d : pd.DataFrame
        2D profile data containing Metadata_patient, Metadata_Well,
        and feature columns.
    df_3d : pd.DataFrame
        3D profile data with the same metadata structure.
    results_dir : pathlib.Path
        Directory to save the output parquet file.
    profile_type : str
        Label for the output file (e.g., "organoid_agg", "sc_agg").
    """
    merge_cols = ["Metadata_patient", "Metadata_tumor", "Metadata_Well"]

    # Find patient-tumor-well combinations present in both
    keys_2d = df_2d[merge_cols].drop_duplicates()
    keys_3d = df_3d[merge_cols].drop_duplicates()
    shared_keys = keys_2d.merge(keys_3d, on=merge_cols, how="inner")

    if len(shared_keys) == 0:
        raise ValueError("No matched patient-tumor-well combinations.")

    # Filter and align both dataframes to shared keys
    df_2d_matched = (
        df_2d.merge(shared_keys, on=merge_cols, how="inner")
        .sort_values(merge_cols)
        .reset_index(drop=True)
    )
    df_3d_matched = (
        df_3d.merge(shared_keys, on=merge_cols, how="inner")
        .sort_values(merge_cols)
        .reset_index(drop=True)
    )

    assert list(df_2d_matched[merge_cols].apply(tuple, axis=1)) == list(
        df_3d_matched[merge_cols].apply(tuple, axis=1)
    ), "Patient-tumor-well combinations are not aligned after sorting"

    # Extract and standardize feature arrays
    features_2d = [
        col for col in df_2d_matched.columns if not col.startswith("Metadata_")
    ]
    features_3d = [
        col for col in df_3d_matched.columns if not col.startswith("Metadata_")
    ]

    # Drop features that are all NaN
    df_2d_clean = df_2d_matched[features_2d].dropna(axis=1, how="all")
    df_3d_clean = df_3d_matched[features_3d].dropna(axis=1, how="all")

    features_2d = list(df_2d_clean.columns)
    features_3d = list(df_3d_clean.columns)

    arr_2d_z = standardize(df_2d_matched[features_2d].values)
    arr_3d_z = standardize(df_3d_matched[features_3d].values)
    n_wells = arr_2d_z.shape[0]

    # Compute correlation via matrix multiplication (NaN-aware)
    df_2d_std = pd.DataFrame(arr_2d_z, columns=features_2d)
    df_3d_std = pd.DataFrame(arr_3d_z, columns=features_3d)
    combined = pd.concat([df_2d_std, df_3d_std], axis=1)
    full_corr = combined.corr(method="pearson")
    corr_matrix = full_corr.loc[features_2d, features_3d].values

    # Convert to long-form and save
    corr_df = pd.DataFrame(corr_matrix, index=features_2d, columns=features_3d)
    corr_long = (
        corr_df.reset_index()
        .melt(id_vars="index", var_name="feature_3d", value_name="pearson_r")
        .rename(columns={"index": "feature_2d"})
    )
    output_path = results_dir / f"all_patients_{profile_type}_correlation.parquet"
    corr_long.to_parquet(output_path, index=False)


# In[9]:


# All patients combined
compute_and_save_correlation_all(
    df_2d_organoid, df_3d_organoid, results_dir, "organoid_agg"
)
compute_and_save_correlation_all(df_2d_sc, df_3d_sc, results_dir, "sc_agg")
