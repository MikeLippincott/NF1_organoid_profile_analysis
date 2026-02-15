#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import pandas as pd
from sklearn.decomposition import PCA

# In[2]:


# Configuration

root_dir = pathlib.Path.cwd().parent.parent

# PCA parameters
N_COMPONENTS = 2


# In[3]:


# Data Paths

# Create a comprehensive dictionary for all dimension and projection combinations
data_dict = {
    "2D_max_projection": {
        "sc": {
            "input": pathlib.Path(
                f"{root_dir}/data/2D_profiles/all_patient_profiles/max_projection/sc_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/pca/2D/max_projection/sc_pca.parquet"
            ).resolve(),
        },
        "sc_fs": {
            "input": pathlib.Path(
                f"{root_dir}/data/2D_profiles/all_patient_profiles/max_projection/sc_fs_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/pca/2D/max_projection/sc_fs_pca.parquet"
            ).resolve(),
        },
        "sc_agg": {
            "input": pathlib.Path(
                f"{root_dir}/data/2D_profiles/all_patient_profiles/max_projection/sc_agg_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/pca/2D/max_projection/sc_agg_pca.parquet"
            ).resolve(),
        },
        "sc_consensus": {
            "input": pathlib.Path(
                f"{root_dir}/data/2D_profiles/all_patient_profiles/max_projection/sc_consensus_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/pca/2D/max_projection/sc_consensus_pca.parquet"
            ).resolve(),
        },
    },
    "2D_middle_slice": {
        "sc": {
            "input": pathlib.Path(
                f"{root_dir}/data/2D_profiles/all_patient_profiles/middle_slice/sc_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/pca/2D/middle_slice/sc_pca.parquet"
            ).resolve(),
        },
        "sc_fs": {
            "input": pathlib.Path(
                f"{root_dir}/data/2D_profiles/all_patient_profiles/middle_slice/sc_fs_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/pca/2D/middle_slice/sc_fs_pca.parquet"
            ).resolve(),
        },
        "sc_agg": {
            "input": pathlib.Path(
                f"{root_dir}/data/2D_profiles/all_patient_profiles/middle_slice/sc_agg_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/pca/2D/middle_slice/sc_agg_pca.parquet"
            ).resolve(),
        },
        "sc_consensus": {
            "input": pathlib.Path(
                f"{root_dir}/data/2D_profiles/all_patient_profiles/middle_slice/sc_consensus_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/pca/2D/middle_slice/sc_consensus_pca.parquet"
            ).resolve(),
        },
    },
    "3D": {
        "sc": {
            "input": pathlib.Path(
                f"{root_dir}/data/3D_profiles/all_patient_profiles/sc_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/pca/3D/sc_pca.parquet"
            ).resolve(),
        },
        "sc_fs": {
            "input": pathlib.Path(
                f"{root_dir}/data/3D_profiles/all_patient_profiles/sc_fs_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/pca/3D/sc_fs_pca.parquet"
            ).resolve(),
        },
        "sc_agg": {
            "input": pathlib.Path(
                f"{root_dir}/data/3D_profiles/all_patient_profiles/sc_agg_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/pca/3D/sc_agg_pca.parquet"
            ).resolve(),
        },
        "sc_consensus": {
            "input": pathlib.Path(
                f"{root_dir}/data/3D_profiles/all_patient_profiles/sc_consensus_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/pca/3D/sc_consensus_pca.parquet"
            ).resolve(),
        },
    },
}

# Create output directories for all combinations
for projection_key in data_dict:
    data_dict[projection_key]["sc"]["output"].parent.mkdir(parents=True, exist_ok=True)


# In[ ]:


for projection_key in data_dict:
    for profile_type in data_dict[projection_key]:
        # Load the parquet file
        df = pd.read_parquet(data_dict[projection_key][profile_type]["input"])

        # Separate metadata columns from feature columns
        # Metadata columns start with "Metadata_"
        metadata_columns = [col for col in df.columns if "Metadata_" in col]

        # Keep a copy of metadata for later merging with PCA results
        metadata_df = df[metadata_columns].copy()

        # Extract only feature columns (drop all metadata)
        features_df = df.drop(columns=metadata_columns, errors="ignore")

        # Handle NaN values - drop rows with any NaN
        features_df = features_df.dropna(axis=0, how="any")

        # Update metadata to match cleaned features
        metadata_df = metadata_df.loc[features_df.index].reset_index(drop=True)
        features_df = features_df.reset_index(drop=True)

        # Initialize and fit PCA model
        pca_model = PCA(n_components=N_COMPONENTS)
        pca_model.fit(features_df)

        # Transform the data to PCA space
        pca_embeddings = pca_model.transform(features_df)

        # Create DataFrame with PCA coordinates
        pca_columns = [f"PC{i}" for i in range(N_COMPONENTS)]
        pca_df = pd.DataFrame(pca_embeddings, columns=pca_columns)

        # Combine metadata with PCA coordinates
        pca_df = pd.concat([metadata_df, pca_df], axis=1)

        # Save results to parquet file
        pca_df.to_parquet(
            data_dict[projection_key][profile_type]["output"], index=False
        )
