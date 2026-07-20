#!/usr/bin/env python
# coding: utf-8

# This notebooks looks at and finds the inter and intra patient distances for each patient profile.
# The goal of this is to find how similar or different profiles are within and across patients.
#

# In[1]:


import pathlib
from typing import Union

import numpy as np
import pandas as pd
import scipy
import scipy.spatial.distance
import sklearn.metrics.pairwise

# Get the current working directory
cwd = pathlib.Path.cwd()

if (cwd / ".git").is_dir():
    root_dir = cwd

else:
    root_dir = None
    for parent in cwd.parents:
        if (parent / ".git").is_dir():
            root_dir = parent
            break

# Check if a Git root directory was found
if root_dir is None:
    raise FileNotFoundError("No Git root directory found.")


# ## Load Data / Set Paths

# In[2]:


# Define input paths
path_2d_organoid = pathlib.Path(
    f"{root_dir}/data/2D_profiles/all_patient_profiles/max_projection/organoid_agg_profiles.parquet"
).resolve(strict=True)
path_2d_sc = pathlib.Path(
    f"{root_dir}/data/2D_profiles/all_patient_profiles/max_projection/sc_agg_profiles.parquet"
).resolve(strict=True)
path_3d_organoid = pathlib.Path(
    f"{root_dir}/data/3D_profiles/all_patient_profiles/organoid_agg_profiles.parquet"
).resolve(strict=True)
path_3d_sc = pathlib.Path(
    f"{root_dir}/data/3D_profiles/all_patient_profiles/sc_agg_profiles.parquet"
).resolve(strict=True)

# Define output paths
results_dir = pathlib.Path(f"{root_dir}/2.2d_vs_3d_analysis/results/mAP").resolve()
results_dir.mkdir(parents=True, exist_ok=True)
dist_results_dir = pathlib.Path(
    f"{root_dir}/2.2d_vs_3d_analysis/results/distance_metrics"
).resolve()
dist_results_dir.mkdir(parents=True, exist_ok=True)

# Load data
df_2d_organoid = pd.read_parquet(path_2d_organoid)
df_2d_sc = pd.read_parquet(path_2d_sc)
df_3d_organoid = pd.read_parquet(path_3d_organoid)
df_3d_sc = pd.read_parquet(path_3d_sc)

# Create treatment_dose column (combines treatment name and dose for dose-specific analysis)
for df in [df_2d_organoid, df_3d_organoid, df_2d_sc, df_3d_sc]:
    df["Metadata_treatment_dose"] = (
        df["Metadata_treatment"]
        + "_"
        + df["Metadata_dose"].fillna(0).astype(float).astype(int).astype(str)
    )

# Define metadata columns
metadata_cols_2d_organoid = [
    col for col in df_2d_organoid.columns if col.startswith("Metadata_")
]
metadata_cols_3d_organoid = [
    col for col in df_3d_organoid.columns if col.startswith("Metadata_")
]
metadata_cols_2d_sc = [col for col in df_2d_sc.columns if col.startswith("Metadata_")]
metadata_cols_3d_sc = [col for col in df_3d_sc.columns if col.startswith("Metadata_")]


# ## Define the functions

# In[3]:


def calculate_intra_patient_distance_metrics(
    df: pd.DataFrame,
    metadata_columns: list,
    col_for_reference: str = "Metadata_treatment",
    reference_group: str = "DMSO",
    output_path: pathlib.Path = None,
) -> Union[None, pd.DataFrame]:
    """
    Calculate intra-patient cosine distance metrics between each treatment
    and the reference group.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing profiles and metadata.
    metadata_columns : list
        List of metadata columns to exclude from distance calculations.
    col_for_reference : str
        Column name for treatment grouping.
    reference_group : str
        Reference treatment label.
    output_path : pathlib.Path, optional
        Path to save results as parquet.

    Returns
    ----------
    Union[None, pd.DataFrame]
        DataFrame with distance metrics per patient per treatment.
    """
    results = []
    for patient in df["Metadata_patient"].unique():
        patient_df = df[df["Metadata_patient"] == patient]
        for drug in patient_df[col_for_reference].unique():
            if drug == reference_group:
                continue

            drug_df = patient_df.loc[patient_df[col_for_reference] == drug].copy()
            dmso_df = patient_df.loc[
                patient_df[col_for_reference] == reference_group
            ].copy()

            drug_features = drug_df.drop(columns=metadata_columns)
            dmso_features = dmso_df.drop(columns=metadata_columns)

            valid_cols = drug_features.columns[
                ~drug_features.isna().all() & ~dmso_features.isna().all()
            ]
            drug_features = drug_features[valid_cols].fillna(0)
            dmso_features = dmso_features[valid_cols].fillna(0)

            if drug_features.shape[0] == 0 or dmso_features.shape[0] == 0:
                continue

            cosine_dist = sklearn.metrics.pairwise.cosine_distances(
                dmso_features.values, drug_features.values
            ).reshape(-1)

            results.append(
                {
                    "Metadata_patient": patient,
                    col_for_reference: drug,
                    "cosine_distance_mean": cosine_dist.mean(),
                    "cosine_distance_std": cosine_dist.std(),
                }
            )

    output_df = pd.DataFrame(results)
    if output_path is not None:
        output_df.to_parquet(output_path, index=False)
        return None
    return output_df


def calculate_inter_patient_distance_metrics(
    df: pd.DataFrame,
    metadata_columns: list,
    col_for_reference: str = "Metadata_treatment",
    reference_group: str = "DMSO",
    output_path: pathlib.Path = None,
) -> Union[None, pd.DataFrame]:
    """
    Calculate inter-patient cosine distance metrics between each treatment
    and the reference group across all patients.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing profiles and metadata.
    metadata_columns : list
        List of metadata columns to exclude from distance calculations.
    col_for_reference : str
        Column name for treatment grouping.
    reference_group : str
        Reference treatment label.
    output_path : pathlib.Path, optional
        Path to save results as parquet.

    Returns
    ----------
    Union[None, pd.DataFrame]
        DataFrame with distance metrics per treatment across all patients.
    """
    results = []
    for drug in df[col_for_reference].unique():
        if drug == reference_group:
            continue

        drug_df = df.loc[df[col_for_reference] == drug].copy()
        dmso_df = df.loc[df[col_for_reference] == reference_group].copy()

        drug_features = drug_df.drop(columns=metadata_columns)
        dmso_features = dmso_df.drop(columns=metadata_columns)

        valid_cols = drug_features.columns[
            ~drug_features.isna().all() & ~dmso_features.isna().all()
        ]
        drug_features = drug_features[valid_cols].fillna(0)
        dmso_features = dmso_features[valid_cols].fillna(0)

        if drug_features.shape[0] == 0 or dmso_features.shape[0] == 0:
            continue

        cosine_dist = sklearn.metrics.pairwise.cosine_distances(
            dmso_features.values, drug_features.values
        ).reshape(-1)

        results.append(
            {
                col_for_reference: drug,
                "cosine_distance_mean": cosine_dist.mean(),
                "cosine_distance_std": cosine_dist.std(),
            }
        )

    output_df = pd.DataFrame(results)
    if output_path is not None:
        output_df.to_parquet(output_path, index=False)
        return None
    return output_df


# ## Run the functions

# In[4]:


# Distance metrics output directory
dist_results_dir = pathlib.Path(
    f"{root_dir}/2.2d_vs_3d_analysis/results/distance_metrics"
).resolve()
dist_results_dir.mkdir(parents=True, exist_ok=True)

# Intra-patient distance metrics
dist_2d_organoid_intra = calculate_intra_patient_distance_metrics(
    df_2d_organoid,
    metadata_columns=metadata_cols_2d_organoid,
    col_for_reference="Metadata_treatment_dose",
    reference_group="DMSO_1",
    output_path=dist_results_dir / "2d_organoid_intra_patient_cosine_distance.parquet",
)

dist_3d_organoid_intra = calculate_intra_patient_distance_metrics(
    df_3d_organoid,
    metadata_columns=metadata_cols_3d_organoid,
    col_for_reference="Metadata_treatment_dose",
    reference_group="DMSO_1",
    output_path=dist_results_dir / "3d_organoid_intra_patient_cosine_distance.parquet",
)

dist_2d_sc_intra = calculate_intra_patient_distance_metrics(
    df_2d_sc,
    metadata_columns=metadata_cols_2d_sc,
    col_for_reference="Metadata_treatment_dose",
    reference_group="DMSO_1",
    output_path=dist_results_dir / "2d_sc_intra_patient_cosine_distance.parquet",
)

dist_3d_sc_intra = calculate_intra_patient_distance_metrics(
    df_3d_sc,
    metadata_columns=metadata_cols_3d_sc,
    col_for_reference="Metadata_treatment_dose",
    reference_group="DMSO_1",
    output_path=dist_results_dir / "3d_sc_intra_patient_cosine_distance.parquet",
)

# Inter-patient distance metrics
dist_2d_organoid_inter = calculate_inter_patient_distance_metrics(
    df_2d_organoid,
    metadata_columns=metadata_cols_2d_organoid,
    col_for_reference="Metadata_treatment_dose",
    reference_group="DMSO_1",
    output_path=dist_results_dir / "2d_organoid_inter_patient_cosine_distance.parquet",
)

dist_3d_organoid_inter = calculate_inter_patient_distance_metrics(
    df_3d_organoid,
    metadata_columns=metadata_cols_3d_organoid,
    col_for_reference="Metadata_treatment_dose",
    reference_group="DMSO_1",
    output_path=dist_results_dir / "3d_organoid_inter_patient_cosine_distance.parquet",
)

dist_2d_sc_inter = calculate_inter_patient_distance_metrics(
    df_2d_sc,
    metadata_columns=metadata_cols_2d_sc,
    col_for_reference="Metadata_treatment_dose",
    reference_group="DMSO_1",
    output_path=dist_results_dir / "2d_sc_inter_patient_cosine_distance.parquet",
)

dist_3d_sc_inter = calculate_inter_patient_distance_metrics(
    df_3d_sc,
    metadata_columns=metadata_cols_3d_sc,
    col_for_reference="Metadata_treatment_dose",
    reference_group="DMSO_1",
    output_path=dist_results_dir / "3d_sc_inter_patient_cosine_distance.parquet",
)
