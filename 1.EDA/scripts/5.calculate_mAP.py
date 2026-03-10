#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
from typing import Union

import numpy as np
import pandas as pd
from copairs import map
from copairs.matching import assign_reference_index
from notebook_init_utils.notebook_init_utils import init_notebook

root_dir, in_notebook = init_notebook()


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
results_dir = pathlib.Path(f"{root_dir}/1.EDA/results/mAP").resolve()
results_dir.mkdir(parents=True, exist_ok=True)

# Load data
df_2d_organoid = pd.read_parquet(path_2d_organoid)
df_2d_sc = pd.read_parquet(path_2d_sc)
df_3d_organoid = pd.read_parquet(path_3d_organoid)
df_3d_sc = pd.read_parquet(path_3d_sc)

# Create treatment_dose column
for df in [df_2d_organoid, df_3d_organoid, df_2d_sc, df_3d_sc]:
    df["Metadata_treatment_dose"] = (
        df["Metadata_treatment"] + "_" + df["Metadata_dose"].astype(str)
    )

print(
    df_2d_organoid[df_2d_organoid["Metadata_treatment"] == "DMSO"][
        "Metadata_treatment_dose"
    ].unique()
)


# In[3]:


metadata_cols_2d_organoid = [
    col for col in df_2d_organoid.columns if col.startswith("Metadata_")
]
metadata_cols_3d_organoid = [
    col for col in df_3d_organoid.columns if col.startswith("Metadata_")
]
metadata_cols_2d_sc = [col for col in df_2d_sc.columns if col.startswith("Metadata_")]
metadata_cols_3d_sc = [col for col in df_3d_sc.columns if col.startswith("Metadata_")]


# In[4]:


def calculate_mAP(
    df: pd.DataFrame,
    metadata_columns: list,
    col_for_reference: str = "Metadata_treatment",
    reference_group: str = "DMSO",
) -> Union[None, pd.DataFrame]:
    """
    Calculate mean Average Precision (mAP) for a given DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing the profiles and metadata.
    metadata_columns : list
        List of metadata columns to be used for grouping.
    col_for_reference : str
        Column name to be used for reference profiles.
    reference_group : str
        The value in col_for_reference that indicates reference profiles.

    Returns
    ----------
    Union[None, pd.DataFrame]
        DataFrame with mAP results, or None if calculation fails.
    """
    reference_col = "reference_index"

    df = assign_reference_index(
        df,
        f"{col_for_reference} == '{reference_group}'",
        reference_col=reference_col,
        default_value=-1,
    )
    feature_cols = [
        col
        for col in df.columns
        if not col.startswith("Metadata_") and col != reference_col
    ]
    df = df.drop(columns=[col for col in feature_cols if df[col].isna().any()])
    metadata_columns_with_ref = metadata_columns + [reference_col]
    metadata = df[metadata_columns_with_ref]
    profiles = df.drop(columns=metadata_columns_with_ref).values

    pos_sameby = [col_for_reference, reference_col]
    pos_diffby = []
    neg_sameby = []
    neg_diffby = [col_for_reference, reference_col]

    try:
        df_ap = map.average_precision(
            metadata, profiles, pos_sameby, pos_diffby, neg_sameby, neg_diffby
        )
        df_ap = df_ap.query(f"{col_for_reference} != '{reference_group}'")
        activity_map = map.mean_average_precision(
            df_ap, pos_sameby, null_size=1000000, threshold=0.05, seed=0
        )
        activity_map["-log10(p-value)"] = -activity_map["corrected_p_value"].apply(
            np.log10
        )
        return activity_map
    except Exception as e:
        print(f"Error calculating mAP: {e}")
        return None


def calculate_intra_patient_mAP(
    df: pd.DataFrame,
    metadata_columns: list,
    col_for_reference: str = "Metadata_treatment",
    reference_group: str = "DMSO",
    output_path: Union[None, pathlib.Path] = None,
) -> Union[None, pd.DataFrame]:
    """
    Calculate intra-patient mAP for each treatment per patient.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing the profiles and metadata.
    metadata_columns : list
        List of metadata columns to exclude from profiles.
    col_for_reference : str
        Column name for treatment grouping.
    reference_group : str
        Reference treatment label (e.g., "DMSO").
    output_path : Union[None, pathlib.Path]
        Path to save results as parquet.

    Returns
    ----------
    Union[None, pd.DataFrame]
        DataFrame with mAP results per patient per treatment.
    """
    list_of_dfs = []
    for patient in df["Metadata_patient"].unique():
        patient_df = df.loc[df["Metadata_patient"] == patient, :].copy()
        for drug in patient_df[col_for_reference].unique():
            if drug == reference_group:
                continue
            drug_df = patient_df.loc[patient_df[col_for_reference] == drug, :].copy()
            dmso_df = patient_df.loc[
                patient_df[col_for_reference] == reference_group, :
            ].copy()
            drug_df = pd.concat([drug_df, dmso_df], ignore_index=True)

            mAP = calculate_mAP(
                drug_df,
                metadata_columns=metadata_columns,
                col_for_reference=col_for_reference,
                reference_group=reference_group,
            )
            if mAP is not None:
                mAP["Metadata_patient"] = patient
                mAP["Metadata_treatment"] = drug
                list_of_dfs.append(mAP)

    if len(list_of_dfs) == 0:
        print("No mAP results computed.")
        return None

    output_df = pd.concat(list_of_dfs, ignore_index=True)
    if output_path is not None:
        output_df.to_parquet(output_path, index=False)
    return output_df


def calculate_inter_patient_mAP(
    df: pd.DataFrame,
    metadata_columns: list,
    col_for_reference: str = "Metadata_treatment",
    reference_group: str = "DMSO",
    output_path: Union[None, pathlib.Path] = None,
) -> Union[None, pd.DataFrame]:
    """
    Calculate inter-patient mAP across all patients for each treatment.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing the profiles and metadata.
    metadata_columns : list
        List of metadata columns to exclude from profiles.
    col_for_reference : str
        Column name for treatment grouping.
    reference_group : str
        Reference treatment label (e.g., "DMSO").
    output_path : Union[None, pathlib.Path]
        Path to save results as parquet.

    Returns
    ----------
    Union[None, pd.DataFrame]
        DataFrame with mAP results per treatment across all patients.
    """
    list_of_dfs = []
    for drug in df[col_for_reference].unique():
        if drug == reference_group:
            continue
        drug_df = df.loc[df[col_for_reference] == drug, :].copy()
        dmso_df = df.loc[df[col_for_reference] == reference_group, :].copy()
        drug_df = pd.concat([drug_df, dmso_df], ignore_index=True)

        mAP = calculate_mAP(
            drug_df,
            metadata_columns=metadata_columns,
            col_for_reference=col_for_reference,
            reference_group=reference_group,
        )
        if mAP is not None:
            mAP["Metadata_treatment"] = drug
            mAP["Metadata_patient"] = "all_patients"
            list_of_dfs.append(mAP)

    if len(list_of_dfs) == 0:
        print("No mAP results computed.")
        return None

    output_df = pd.concat(list_of_dfs, ignore_index=True)
    if output_path is not None:
        output_df.to_parquet(output_path, index=False)
    return output_df


# In[5]:


# Intra-patient mAP by treatment and dose
print("Computing 2D organoid intra-patient mAP by dose...")
mAP_2d_organoid_intra = calculate_intra_patient_mAP(
    df_2d_organoid,
    metadata_columns=metadata_cols_2d_organoid,
    col_for_reference="Metadata_treatment_dose",
    reference_group="DMSO_1",
    output_path=results_dir / "2d_organoid_intra_patient_mAP_by_dose.parquet",
)

print("Computing 3D organoid intra-patient mAP by dose...")
mAP_3d_organoid_intra = calculate_intra_patient_mAP(
    df_3d_organoid,
    metadata_columns=metadata_cols_3d_organoid,
    col_for_reference="Metadata_treatment_dose",
    reference_group="DMSO_1",
    output_path=results_dir / "3d_organoid_intra_patient_mAP_by_dose.parquet",
)

print("Computing 2D SC intra-patient mAP by dose...")
mAP_2d_sc_intra = calculate_intra_patient_mAP(
    df_2d_sc,
    metadata_columns=metadata_cols_2d_sc,
    col_for_reference="Metadata_treatment_dose",
    reference_group="DMSO_1",
    output_path=results_dir / "2d_sc_intra_patient_mAP_by_dose.parquet",
)

print("Computing 3D SC intra-patient mAP by dose...")
mAP_3d_sc_intra = calculate_intra_patient_mAP(
    df_3d_sc,
    metadata_columns=metadata_cols_3d_sc,
    col_for_reference="Metadata_treatment_dose",
    reference_group="DMSO_1",
    output_path=results_dir / "3d_sc_intra_patient_mAP_by_dose.parquet",
)


# In[6]:


# Inter-patient mAP by treatment and dose
print("Computing 2D organoid inter-patient mAP by dose...")
mAP_2d_organoid_inter = calculate_inter_patient_mAP(
    df_2d_organoid,
    metadata_columns=metadata_cols_2d_organoid,
    col_for_reference="Metadata_treatment_dose",
    reference_group="DMSO_1",
    output_path=results_dir / "2d_organoid_inter_patient_mAP_by_dose.parquet",
)

print("Computing 3D organoid inter-patient mAP by dose...")
mAP_3d_organoid_inter = calculate_inter_patient_mAP(
    df_3d_organoid,
    metadata_columns=metadata_cols_3d_organoid,
    col_for_reference="Metadata_treatment_dose",
    reference_group="DMSO_1",
    output_path=results_dir / "3d_organoid_inter_patient_mAP_by_dose.parquet",
)

print("Computing 2D SC inter-patient mAP by dose...")
mAP_2d_sc_inter = calculate_inter_patient_mAP(
    df_2d_sc,
    metadata_columns=metadata_cols_2d_sc,
    col_for_reference="Metadata_treatment_dose",
    reference_group="DMSO_1",
    output_path=results_dir / "2d_sc_inter_patient_mAP_by_dose.parquet",
)

print("Computing 3D SC inter-patient mAP by dose...")
mAP_3d_sc_inter = calculate_inter_patient_mAP(
    df_3d_sc,
    metadata_columns=metadata_cols_3d_sc,
    col_for_reference="Metadata_treatment_dose",
    reference_group="DMSO_1",
    output_path=results_dir / "3d_sc_inter_patient_mAP_by_dose.parquet",
)


# In[9]:


print(
    "3D organoid DMSO check:",
    (df_3d_organoid["Metadata_treatment_dose"] == "DMSO_1").sum(),
)
print("3D SC DMSO check:", (df_3d_sc["Metadata_treatment_dose"] == "DMSO_1").sum())
print()
print(
    "3D organoid unique treatments:",
    df_3d_organoid["Metadata_treatment_dose"].unique()[:5],
)


# In[11]:


import inspect

print(inspect.signature(calculate_intra_patient_mAP))
