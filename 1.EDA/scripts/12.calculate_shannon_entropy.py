#!/usr/bin/env python
# coding: utf-8

# # Calculate Shannon Entropy
#
# Compute per-feature Shannon entropy for 2D (max projection) and 3D organoid and single-cell aggregated profiles.
# Entropy is estimated by discretizing each feature into 50 histogram bins and computing the entropy of the resulting probability distribution.

# In[1]:


from typing import List

import numpy as np
import pandas as pd
from notebook_init_utils.notebook_init_utils import init_notebook
from scipy.stats import entropy

root_dir, in_notebook = init_notebook()
print(f"Root directory: {root_dir}")


# In[2]:


# Define input/output paths
_2d_dir = root_dir / "data" / "2D_profiles" / "all_patient_profiles" / "max_projection"
_3d_dir = root_dir / "data" / "3D_profiles" / "all_patient_profiles"
_entropy_dir = root_dir / "1.EDA" / "results" / "entropy"

data_dict = {
    # 2D organoid
    "2D_organoid": {
        "input": (_2d_dir / "organoid_profiles.parquet").resolve(strict=True),
        "output": (_entropy_dir / "2D_organoid_entropy.parquet").resolve(),
        "modality": "2D",
        "profile_type": "organoid",
    },
    "2D_organoid_fs": {
        "input": (_2d_dir / "organoid_fs_profiles.parquet").resolve(strict=True),
        "output": (_entropy_dir / "2D_organoid_fs_entropy.parquet").resolve(),
        "modality": "2D",
        "profile_type": "organoid",
    },
    "2D_organoid_agg": {
        "input": (_2d_dir / "organoid_agg_profiles.parquet").resolve(strict=True),
        "output": (_entropy_dir / "2D_organoid_agg_entropy.parquet").resolve(),
        "modality": "2D",
        "profile_type": "organoid",
    },
    # 2D single-cell
    "2D_sc": {
        "input": (_2d_dir / "sc_profiles.parquet").resolve(strict=True),
        "output": (_entropy_dir / "2D_sc_entropy.parquet").resolve(),
        "modality": "2D",
        "profile_type": "sc",
    },
    "2D_sc_fs": {
        "input": (_2d_dir / "sc_fs_profiles.parquet").resolve(strict=True),
        "output": (_entropy_dir / "2D_sc_fs_entropy.parquet").resolve(),
        "modality": "2D",
        "profile_type": "sc",
    },
    "2D_sc_agg": {
        "input": (_2d_dir / "sc_agg_profiles.parquet").resolve(strict=True),
        "output": (_entropy_dir / "2D_sc_agg_entropy.parquet").resolve(),
        "modality": "2D",
        "profile_type": "sc",
    },
    # 3D organoid
    "3D_organoid": {
        "input": (_3d_dir / "organoid_profiles.parquet").resolve(strict=True),
        "output": (_entropy_dir / "3D_organoid_entropy.parquet").resolve(),
        "modality": "3D",
        "profile_type": "organoid",
    },
    "3D_organoid_fs": {
        "input": (_3d_dir / "organoid_fs_profiles.parquet").resolve(strict=True),
        "output": (_entropy_dir / "3D_organoid_fs_entropy.parquet").resolve(),
        "modality": "3D",
        "profile_type": "organoid",
    },
    "3D_organoid_agg": {
        "input": (_3d_dir / "organoid_agg_profiles.parquet").resolve(strict=True),
        "output": (_entropy_dir / "3D_organoid_agg_entropy.parquet").resolve(),
        "modality": "3D",
        "profile_type": "organoid",
    },
    # 3D single-cell
    "3D_sc": {
        "input": (_3d_dir / "sc_profiles.parquet").resolve(strict=True),
        "output": (_entropy_dir / "3D_sc_entropy.parquet").resolve(),
        "modality": "3D",
        "profile_type": "sc",
    },
    "3D_sc_fs": {
        "input": (_3d_dir / "sc_fs_profiles.parquet").resolve(strict=True),
        "output": (_entropy_dir / "3D_sc_fs_entropy.parquet").resolve(),
        "modality": "3D",
        "profile_type": "sc",
    },
    "3D_sc_agg": {
        "input": (_3d_dir / "sc_agg_profiles.parquet").resolve(strict=True),
        "output": (_entropy_dir / "3D_sc_agg_entropy.parquet").resolve(),
        "modality": "3D",
        "profile_type": "sc",
    },
}


# In[3]:


def compute_feature_entropy(
    df: pd.DataFrame,
    n_bins: int = 50,
) -> pd.DataFrame:
    """Compute Shannon entropy for each feature column in the DataFrame.

    Feature columns are identified as those without the ``Metadata_`` prefix.
    Each feature's distribution is discretized into ``n_bins`` histogram bins,
    and Shannon entropy (base 2, in bits) is computed from the resulting
    probability distribution.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing both metadata and feature columns.
        Metadata columns are expected to carry the ``Metadata_`` prefix.
    n_bins : int, optional
        Number of histogram bins used for discretization. Default is 50.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ``feature_name`` (str) and ``entropy`` (float).
        Rows with no valid values are excluded.
    """
    # Separate feature columns from metadata columns
    feature_columns: List[str] = [
        col for col in df.columns if not col.startswith("Metadata_")
    ]

    results: Dict[str, List] = {"feature_name": [], "entropy": []}
    for feature in feature_columns:
        # Drop NaNs and skip empty features
        values = df[feature].dropna().values
        if len(values) == 0:
            continue

        # Discretize into histogram bins; entropy() normalizes counts internally
        counts, _ = np.histogram(values, bins=n_bins)
        feat_entropy: float = entropy(counts, base=2)

        results["feature_name"].append(feature)
        results["entropy"].append(feat_entropy)

    return pd.DataFrame(results)


# In[4]:


# Compute and save entropy for each profile type
for key, config in data_dict.items():
    df = pd.read_parquet(config["input"])

    entropy_df = compute_feature_entropy(df, n_bins=50)
    entropy_df["imaging_modality"] = config["modality"]
    entropy_df["profile_type"] = config["profile_type"]

    config["output"].parent.mkdir(parents=True, exist_ok=True)
    entropy_df.to_parquet(config["output"], index=False)
    print(config["output"])


# In[5]:


def compute_feature_entropy_per_patient(
    df: pd.DataFrame,
    patient_col: str = "Metadata_patient",
    n_bins: int = 50,
) -> pd.DataFrame:
    """Compute Shannon entropy per patient for each feature column.

    For each unique value in ``patient_col``, the feature values belonging
    to that patient are discretized into ``n_bins`` histogram bins and
    Shannon entropy (base 2) is computed from the resulting distribution.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with metadata and feature columns.
    patient_col : str, optional
        Name of the patient identifier column. Default is ``"Metadata_patient"``.
    n_bins : int, optional
        Number of histogram bins. Default is 50.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ``patient``, ``feature_name``, and ``entropy``.
    """
    feature_columns: List[str] = [
        col for col in df.columns if not col.startswith("Metadata_")
    ]

    records: List[dict] = []
    for patient, group in df.groupby(patient_col):
        for feature in feature_columns:
            values = group[feature].dropna().values
            if len(values) == 0:
                continue
            counts, _ = np.histogram(values, bins=n_bins)
            feat_entropy: float = entropy(counts, base=2)
            records.append(
                {"patient": patient, "feature_name": feature, "entropy": feat_entropy}
            )

    return pd.DataFrame(records)


# In[6]:


# Compute and save per-patient entropy for each profile type
for key, config in data_dict.items():
    df = pd.read_parquet(config["input"])

    # Some 2D normal profiles lack Metadata_patient; derive it from Metadata_patient_tumor
    if "Metadata_patient" not in df.columns and "Metadata_patient_tumor" in df.columns:
        df["Metadata_patient"] = (
            df["Metadata_patient_tumor"].str.rsplit("_", n=1).str[0]
        )

    per_patient_output = config["output"].parent / config["output"].name.replace(
        "_entropy.parquet", "_per_patient_entropy.parquet"
    )

    per_patient_df = compute_feature_entropy_per_patient(df, n_bins=50)
    per_patient_df["imaging_modality"] = config["modality"]
    per_patient_df["profile_type"] = config["profile_type"]

    per_patient_df.to_parquet(per_patient_output, index=False)
    print(per_patient_output)
