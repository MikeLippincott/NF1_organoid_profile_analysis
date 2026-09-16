#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pathlib
import warnings

import pandas as pd
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")  # Ignore all warnings
warnings.simplefilter("ignore")  # Additional suppression method

from notebook_init_utils.notebook_init_utils import init_notebook

root_dir, in_notebook = init_notebook()


# In[ ]:


profile_dict = {
    "organoid_fs": {
        "input_profile_path": pathlib.Path(
            root_dir, "4.linear_modeling/results/linear_modeling/organoid_norm.parquet"
        ).resolve(strict=True),
        "metadata_columns": [
            "patient",
            "object_id",
            "unit",
            "dose",
            "treatment",
            "Target",
            "Class",
            "image_set",
            "Well",
            "Therapeutic_Categories",
            "single_cell_count",
        ],
    },
    "single_cell_fs": {
        "input_profile_path": pathlib.Path(
            root_dir, "4.linear_modeling/results/linear_modeling/sc_norm.parquet"
        ).resolve(strict=True),
        "metadata_columns": [
            "patient",
            "object_id",
            "unit",
            "dose",
            "treatment",
            "Target",
            "Class",
            "image_set",
            "Well",
            "Therapeutic_Categories",
            "parent_organoid",
        ],
    },
}


# ## Filter significant features
# pvalue threshold is set to 0.05 - statistically significant features
# rsquared threshold is set to 0.5 - the explained variance is at least 50% of the total variance
# rsquared adjusted threshold is set to positive values - the model performs better than the mean
#

# ### Single Cell

# In[ ]:


df = pd.read_parquet(
    profile_dict["single_cell_fs"]["input_profile_path"],
)
print(df.shape)


# In[ ]:


pvalue_threshold_max = 0.05  # significance threshold for p-values
rsquared_threshold_min = 0.5  # 50% of variance explained by the model
rsquared_adj_threshold_min = 0  # the model performs better than the null model
coefficient_threshold_min = 1  # minimum effect size of 1


# In[ ]:


# filter significant features
df_filtered = df[
    (df["pvalue"] < pvalue_threshold_max)
    & (df["rsquared"] > rsquared_threshold_min)
    & (df["rsquared_adj"] > rsquared_adj_threshold_min)
    & (df["coefficient"].abs() > coefficient_threshold_min)
].copy()
print(df_filtered.shape)
df_filtered.head()


# In[ ]:


df_filtered["treatment"].unique()


# In[ ]:


df_filtered["patient"].unique()


# In[ ]:


df_filtered["feature"].unique()


# ### Organoid

# In[ ]:


df = pd.read_parquet(
    profile_dict["organoid_fs"]["input_profile_path"],
)
print(df.shape)


# In[ ]:


pvalue_threshold_max = 0.05
rsquared_threshold_min = 0.4
rsquared_adj_threshold_min = 0
coefficient_threshold_min = 1


# In[ ]:


# filter significant features
df_filtered = df[
    (df["pvalue"] < pvalue_threshold_max)
    & (df["rsquared"] > rsquared_threshold_min)
    & (df["rsquared_adj"] > rsquared_adj_threshold_min)
    & (df["coefficient"].abs() > coefficient_threshold_min)
].copy()
print(df_filtered.shape)
df_filtered.head()


# In[ ]:


df_filtered["treatment"].unique()


# In[ ]:


df_filtered["patient"].unique()


# In[ ]:


df_filtered["feature"].unique()


# In[ ]:


# In[ ]:
