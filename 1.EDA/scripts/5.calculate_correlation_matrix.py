#!/usr/bin/env python
# coding: utf-8

# This notebook generates the correlation matrices for each of the 2D and 3D patient profiles.
# The profiles cover single cell and organoid profiles:
# - 2D:
#     - maximum projection of the 3D volume
#     - middle slice of the 3D volume
#     - middle 3 slices of the 3D volume maxprojected
# - 3D:
#     - handcrafted 3D features
#     - masked 3D Deep Learning features
#     - nucleocentric 3D Deep Learning features
#     - nucleocentric 2D maxprojected Deep Learning features

# In[1]:


import pathlib

import pandas as pd
from notebook_init_utils import init_notebook

root_dir, in_notebook = init_notebook()

if in_notebook:
    import tqdm.notebook as tqdm
else:
    import tqdm


# In[2]:


# Discover every aggregated / consensus profile file (2D and 3D)
# Single-cell / organoid-level files (profiles, fs_profiles, norm_profile) are
# excluded here: with 20k-77k rows each, a row-by-row correlation matrix would
# blow up to tens of GB per file.
profile_dirs = [
    root_dir / "data/profiles_2D/all_patients",
    root_dir / "data/profiles_3D/all_patients",
]

input_files = []
for base_dir in profile_dirs:
    for file_path in sorted(base_dir.glob("**/*.parquet")):
        if file_path.stem.endswith("agg_profiles") or file_path.stem.endswith(
            "consensus_profiles"
        ):
            input_files.append(file_path)

output_dir = pathlib.Path(f"{root_dir}/1.EDA/results/correlation").resolve()
output_dir.mkdir(parents=True, exist_ok=True)

len(input_files)


# In[3]:


for file_path in tqdm.tqdm(input_files):
    df = pd.read_parquet(file_path)

    # row-by-row (sample-by-sample) correlation matrix over feature columns
    metadata_cols = [col for col in df.columns if col.startswith("Metadata_")]
    feature_cols = [col for col in df.columns if col not in metadata_cols]

    metadata_df = df[metadata_cols].reset_index(drop=True)
    features_df = df[feature_cols].reset_index(drop=True)

    corr_matrix = features_df.T.corr()
    corr_matrix.columns = [f"Sample_{i}" for i in corr_matrix.columns]
    corr_matrix = corr_matrix.reset_index(drop=True)

    # add the metadata back so each row can be identified/grouped when plotting in R
    corr_df = pd.concat([metadata_df, corr_matrix], axis=1)

    dimension = "2D" if "profiles_2D" in file_path.parts else "3D"
    output_name = f"{dimension}_{file_path.parent.name}_{file_path.stem}"
    corr_df.to_parquet(output_dir / output_name, index=False)

    print(f"Saved {output_name}: {corr_df.shape}")
