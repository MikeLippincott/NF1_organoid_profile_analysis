#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import pandas as pd
from notebook_init_utils import init_notebook
from pycytominer import aggregate
from pycytominer.cyto_utils import infer_cp_features

root_dir, in_notebook = init_notebook()

if in_notebook:
    from tqdm.notebook import tqdm
else:
    from tqdm import tqdm


# In[2]:


profile_dict = {
    "organoid_norm": {
        "input_profile_path": pathlib.Path(
            root_dir,
            "data/profiles_3D/all_patients/0.normalized_profiles/organoid_norm_norm_profile.parquet",
        ),
        "output_profile_path": pathlib.Path(
            root_dir, "4.linear_modeling/data/organoid_norm_aggregated_profile.parquet"
        ),
    },
    "single_cell_norm": {
        "input_profile_path": pathlib.Path(
            root_dir,
            "data/profiles_3D/all_patients/0.normalized_profiles/sc_norm_norm_profile.parquet",
        ),
        "output_profile_path": pathlib.Path(
            root_dir, "4.linear_modeling/data/sc_norm_aggregated_profile.parquet"
        ),
    },
}


# In[3]:


aggregate_strata = [
    "Metadata_Biology_PatientTumor",
    "Metadata_Experiment_Well",
    "Metadata_Experiment_Class",
    "Metadata_Experiment_Dose",
    "Metadata_Experiment_Target",
    "Metadata_Experiment_TherapeuticCategories",
    "Metadata_Experiment_Treatment",
    "Metadata_Experiment_Unit",
]


# In[4]:


for profile_name, profile_sub_dict in profile_dict.items():
    input_profile_path = profile_sub_dict["input_profile_path"]
    output_profile_path = profile_sub_dict["output_profile_path"]
    output_profile_path.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(input_profile_path)

    features_columns = [col for col in df.columns if not col.startswith("Metadata_")]
    if profile_name == "organoid_norm":
        features_columns += ["Metadata_Object_OrganoidSingleCellCount"]
        features_columns += ["Metadata_WellOrganoidCount"]

    agg_well_df = aggregate(
        population_df=df,
        strata=aggregate_strata,
        features=features_columns,
        operation="median",
        output_file=output_profile_path,
        output_type="parquet",
    )
