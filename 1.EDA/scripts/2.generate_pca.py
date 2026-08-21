#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from notebook_init_utils import init_notebook
from sklearn.decomposition import PCA

root_dir, in_notebook = init_notebook()

if in_notebook:
    import tqdm.notebook as tqdm
else:
    import tqdm


# In[2]:


# PCA parameters
N_COMPONENTS = 50
OVERWRITE_EXISTING = False


# In[3]:


all_patients_2D_dir = pathlib.Path(
    f"{root_dir}/data/profiles_2D/all_patients/"
).resolve(strict=True)
all_patients_3D_dir = pathlib.Path(
    f"{root_dir}/data/profiles_3D/all_patients/"
).resolve(strict=True)

all_patient_2D_dirs = [x for x in all_patients_2D_dir.iterdir() if x.is_dir()]
all_patient_3D_dirs = [x for x in all_patients_3D_dir.iterdir() if x.is_dir()]
all_patient_2D_dirs = [y for x in all_patient_2D_dirs for y in x.glob("*.parquet")]
all_patient_3D_dirs = [y for x in all_patient_3D_dirs for y in x.glob("*.parquet")]
all_patients_2D_and_3D_dirs = all_patient_2D_dirs + all_patient_3D_dirs
# filter for only fs, agg, consensus
all_patients_2D_and_3D_dirs = [
    x
    for x in all_patients_2D_and_3D_dirs
    if any(substr in x.name for substr in ["fs", "agg", "consensus"])
]


# In[4]:


# generate a dict for each path
data_dict = {
    "profile_type": [],
    "profile_strategy": [],
    "dimension": [],
    "input_path": [],
    "save_path": [],
}
for profile_path in all_patients_2D_and_3D_dirs:
    profile_type = profile_path.stem
    profile_strategy = profile_path.parent.name
    dimension = profile_path.parent.parent.parent.name.replace("profiles_", "")
    data_dict["save_path"].append(
        pathlib.Path(
            f"{root_dir}/1.EDA/results/pca/{dimension}_{profile_strategy}_{profile_type}_embeddings.parquet"
        ).resolve()
    )
    data_dict["profile_type"].append(profile_type)
    data_dict["profile_strategy"].append(profile_strategy)
    data_dict["dimension"].append(dimension)
    data_dict["input_path"].append(profile_path)

for i, (
    profile_type,
    profile_strategy,
    dimension,
    input_path,
    output_file_path,
) in tqdm.tqdm(
    enumerate(
        zip(
            data_dict["profile_type"],
            data_dict["profile_strategy"],
            data_dict["dimension"],
            data_dict["input_path"],
            data_dict["save_path"],
        )
    ),
    total=len(data_dict["input_path"]),
):
    explained_variance_output_file_path = (
        output_file_path.parent
        / f"{dimension}_{profile_strategy}_{profile_type}_explained_variance.parquet"
    )

    if (
        output_file_path.exists()
        and explained_variance_output_file_path.exists()
        and not OVERWRITE_EXISTING
    ):
        continue
    df = pd.read_parquet(input_path)
    # Separate metadata columns from feature columns
    # Metadata columns start with "Metadata_"
    metadata_columns = [col for col in df.columns if "Metadata_" in col]

    # Keep a copy of metadata for later merging with PCA results
    metadata_df = df[metadata_columns].copy()

    # Extract only feature columns (drop all metadata)
    features_df = df.drop(columns=metadata_columns, errors="ignore")

    # Coerce all feature columns to numeric before PCA
    features_df = features_df.apply(pd.to_numeric, errors="coerce")

    # Treat extreme sentinel values as missing before PCA
    # These are far beyond the scale of the normalized features in this dataset.
    max_abs_feature_value = 1e1
    features_df = features_df.mask(features_df.abs() > max_abs_feature_value, np.nan)

    # Replace infs with NaN, then drop rows with any missing values
    features_df = features_df.replace([np.inf, -np.inf], np.nan)
    features_df = features_df.dropna(axis=0)

    # Update metadata to match cleaned features
    metadata_df = metadata_df.loc[features_df.index].reset_index(drop=True)
    features_df = features_df.reset_index(drop=True)

    # Skip files with too few rows/features to fit the requested number of components
    n_components = min(N_COMPONENTS, features_df.shape[0], features_df.shape[1])
    if n_components < 1:
        print(f"Skipping {input_path}: not enough valid data after cleaning.")
        continue

    # Initialize and fit PCA model
    pca_model = PCA(n_components=n_components, svd_solver="full")
    features_array = features_df.to_numpy(dtype=np.float64, copy=False)
    pca_embeddings = pca_model.fit_transform(features_array)

    # Create DataFrame with PCA coordinates
    pca_columns = [f"PC{i}" for i in range(n_components)]
    pca_df = pd.DataFrame(pca_embeddings, columns=pca_columns)

    # Combine metadata with PCA coordinates
    pca_df = pd.concat([metadata_df, pca_df], axis=1)

    # Save results to parquet file
    output_file_path.parent.mkdir(parents=True, exist_ok=True)
    pca_df.to_parquet(
        output_file_path,
        index=False,
    )

    # define a new df for explained variance ratio
    explained_variance_df = pd.DataFrame(
        pca_model.explained_variance_ratio_.reshape(1, -1),
        columns=[f"PC{i}_explained_variance" for i in range(n_components)],
    )
    explained_variance_df.to_parquet(
        explained_variance_output_file_path,
        index=False,
    )
