#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import sys

import pandas as pd
import umap

cwd = pathlib.Path.cwd()

if (cwd / ".git").is_dir():
    root_dir = cwd
else:
    root_dir = None
    for parent in cwd.parents:
        if (parent / ".git").is_dir():
            root_dir = parent
            break
sys.path.append(str(root_dir / "utils"))
from notebook_init_utils import bandicoot_check, init_notebook

root_dir, in_notebook = init_notebook()

profile_base_dir = bandicoot_check(
    pathlib.Path("/home/lippincm/mnt/bandicoot").resolve(), root_dir
)


# In[2]:


# Data Paths

# Create a comprehensive dictionary for all dimension and projection combinations
data_dict = {
    "2D_max_projection": {
        "sc": {
            "input": pathlib.Path(
                f"{root_dir}/data/2D_profiles/all_patient_profiles/max_projection/sc_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/umap/2D/max_projection/sc_umap.parquet"
            ).resolve(),
        },
        "sc_fs": {
            "input": pathlib.Path(
                f"{root_dir}/data/2D_profiles/all_patient_profiles/max_projection/sc_fs_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/umap/2D/max_projection/sc_fs_umap.parquet"
            ).resolve(),
        },
        "sc_agg": {
            "input": pathlib.Path(
                f"{root_dir}/data/2D_profiles/all_patient_profiles/max_projection/sc_agg_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/umap/2D/max_projection/sc_agg_umap.parquet"
            ).resolve(),
        },
        "sc_consensus": {
            "input": pathlib.Path(
                f"{root_dir}/data/2D_profiles/all_patient_profiles/max_projection/sc_consensus_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/umap/2D/max_projection/sc_consensus_umap.parquet"
            ).resolve(),
        },
    },
    "2D_middle_slice": {
        "sc": {
            "input": pathlib.Path(
                f"{root_dir}/data/2D_profiles/all_patient_profiles/middle_slice/sc_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/umap/2D/middle_slice/sc_umap.parquet"
            ).resolve(),
        },
        "sc_fs": {
            "input": pathlib.Path(
                f"{root_dir}/data/2D_profiles/all_patient_profiles/middle_slice/sc_fs_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/umap/2D/middle_slice/sc_fs_umap.parquet"
            ).resolve(),
        },
        "sc_agg": {
            "input": pathlib.Path(
                f"{root_dir}/data/2D_profiles/all_patient_profiles/middle_slice/sc_agg_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/umap/2D/middle_slice/sc_agg_umap.parquet"
            ).resolve(),
        },
        "sc_consensus": {
            "input": pathlib.Path(
                f"{root_dir}/data/2D_profiles/all_patient_profiles/middle_slice/sc_consensus_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/umap/2D/middle_slice/sc_consensus_umap.parquet"
            ).resolve(),
        },
    },
    "3D": {
        "sc": {
            "input": pathlib.Path(
                f"{root_dir}/data/3D_profiles/all_patient_profiles/sc_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/umap/3D/sc_umap.parquet"
            ).resolve(),
        },
        "sc_fs": {
            "input": pathlib.Path(
                f"{root_dir}/data/3D_profiles/all_patient_profiles/sc_fs_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/umap/3D/sc_fs_umap.parquet"
            ).resolve(),
        },
        "sc_agg": {
            "input": pathlib.Path(
                f"{root_dir}/data/3D_profiles/all_patient_profiles/sc_agg_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/umap/3D/sc_agg_umap.parquet"
            ).resolve(),
        },
        "sc_consensus": {
            "input": pathlib.Path(
                f"{root_dir}/data/3D_profiles/all_patient_profiles/sc_consensus_profiles.parquet"
            ).resolve(strict=True),
            "output": pathlib.Path(
                f"{root_dir}/1.EDA/results/umap/3D/sc_consensus_umap.parquet"
            ).resolve(),
        },
    },
}

# Create output directories for all combinations
for projection_key in data_dict:
    data_dict[projection_key]["sc"]["output"].parent.mkdir(parents=True, exist_ok=True)


# In[ ]:


umap_object = umap.UMAP(
    n_neighbors=15, min_dist=0.1, metric="euclidean", random_state=0
)

for projection_key in data_dict:
    for dataset, paths in data_dict[projection_key].items():
        # Load the data
        df = pd.read_parquet(paths["input"])
        metadata_columns = [x for x in df.columns if "Metadata_" in x]
        metadata_df = df.copy()
        metadata_df = df[metadata_columns]
        features_df = df.drop(columns=metadata_columns, errors="ignore")

        # Remove NaN values
        features_df = features_df.dropna(axis=0, how="any")

        # Extract features and apply UMAP
        umap_embedding = umap_object.fit_transform(features_df)

        # Create a DataFrame with UMAP results
        umap_df = pd.DataFrame(umap_embedding, columns=["UMAP1", "UMAP2"])
        umap_df = pd.concat([metadata_df.reset_index(drop=True), umap_df], axis=1)

        # Save the UMAP results
        umap_df.to_parquet(paths["output"], index=False)


# ## Individual umaps

# In[4]:


patients = pd.read_csv(
    pathlib.Path(f"{root_dir}/data/patient_IDs.txt").resolve(strict=True),
    header=None,
    names=["patient"],
)["patient"].to_list()


# In[ ]:


# Individual patients, 2D single cell aggregate and fs only
file_dict = {}

for patient in patients:
    file_dict[patient] = {
        "2D_max_projection": {
            "sc_fs": {
                "input": pathlib.Path(
                    f"{root_dir}/data/2D_profiles/{patient}/3.feature_selected/max_projected_sc.parquet"
                ).resolve(strict=True),
                "output": pathlib.Path(
                    f"{root_dir}/1.EDA/results/umap/patient_results/2D/max_projection/{patient}_sc_fs_umap.parquet"
                ).resolve(),
            },
        },
        "2D_middle_slice": {
            "sc_fs": {
                "input": pathlib.Path(
                    f"{root_dir}/data/2D_profiles/{patient}/3.feature_selected/middle_slice_sc.parquet"
                ).resolve(strict=True),
                "output": pathlib.Path(
                    f"{root_dir}/1.EDA/results/umap/patient_results/2D/middle_slice/{patient}_sc_fs_umap.parquet"
                ).resolve(),
            },
        },
    }

# Create output directories for all combinations
for patient in file_dict:
    for projection_key in file_dict[patient]:
        file_dict[patient][projection_key]["sc_fs"]["output"].parent.mkdir(
            parents=True, exist_ok=True
        )


# In[ ]:


umap_object = umap.UMAP(
    n_neighbors=15, min_dist=0.1, metric="euclidean", random_state=0
)

for patient in file_dict:
    for projection_key in file_dict[patient]:
        for dataset, paths in file_dict[patient][projection_key].items():
            # Load the data
            df = pd.read_parquet(paths["input"])
            metadata_columns = [x for x in df.columns if "Metadata_" in x]
            metadata_df = df.copy()
            metadata_df = df[metadata_columns]
            features_df = df.drop(columns=metadata_columns, errors="ignore")

            # Remove NaN values
            features_df = features_df.dropna(axis=0, how="any")

            # Extract features and apply UMAP
            umap_embedding = umap_object.fit_transform(features_df)

            # Create a DataFrame with UMAP results
            umap_df = pd.DataFrame(umap_embedding, columns=["UMAP1", "UMAP2"])
            umap_df = pd.concat([metadata_df.reset_index(drop=True), umap_df], axis=1)

            # Save the UMAP results
            umap_df.to_parquet(paths["output"], index=False)
