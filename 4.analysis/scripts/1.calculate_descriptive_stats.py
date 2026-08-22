#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
"""Section 1: descriptive stats (mean/median/sd/n) per feature, grouped by
patient x treatment, for all organoid-level profile columns, 2D (3 projections)
and 3D. Uses raw organoid profiles (not sc) as the primary "profile" grain,
since all-column descriptive stats on ~2850-column sc tables would be
unmanageable in size; sc-level detail is covered separately in section 7
(intensity) and section 5 (neighbors).
"""


# In[2]:

import sys
import warnings

import pandas as pd

warnings.filterwarnings("ignore")


# In[3]:


from notebook_init_utils import init_notebook

root_dir, in_notebook = init_notebook()
sys.path.insert(0, str(root_dir / "4.analysis" / "scripts"))


# In[4]:


from utils_analysis import (
    PROJECTION_FILE_PREFIX,
    PROJECTIONS,
    harmonize_metadata,
    list_patient_dirs,
)

results_dir = root_dir / "4.analysis" / "results" / "descriptive_stats"
results_dir.mkdir(parents=True, exist_ok=True)


# In[5]:


def describe_group(df: pd.DataFrame, feature_cols, group_cols) -> pd.DataFrame:
    long = df.melt(
        id_vars=group_cols,
        value_vars=feature_cols,
        var_name="feature",
        value_name="value",
    )
    long["value"] = pd.to_numeric(long["value"], errors="coerce")
    agg = (
        long.groupby(group_cols + ["feature"])["value"]
        .agg(mean="mean", median="median", sd="std", n="count")
        .reset_index()
    )
    return agg


all_2d = []
patients_2d = list_patient_dirs(root_dir / "data" / "profiles_2D")
for projection in PROJECTIONS:
    prefix = PROJECTION_FILE_PREFIX[projection]
    for patient in patients_2d:
        f = (
            root_dir
            / "data"
            / "profiles_2D"
            / patient
            / "5.normalized"
            / f"{prefix}_organoid.parquet"
        )
        if not f.exists():
            continue
        df = pd.read_parquet(f)
        df = harmonize_metadata(df, "2D", patient)
        feature_cols = [c for c in df.columns if not c.startswith("Metadata_")]
        agg = describe_group(
            df, feature_cols, ["Metadata_patient", "Metadata_treatment"]
        )
        agg["projection"] = projection
        agg["modality"] = "2D"
        all_2d.append(agg)
        print(
            f"2D {projection} {patient}: {len(feature_cols)} features, {len(df)} organoids"
        )

desc_2d = pd.concat(all_2d, ignore_index=True)
desc_2d.to_parquet(results_dir / "descriptive_stats_2D.parquet", index=False)
print(f"Wrote {results_dir / 'descriptive_stats_2D.parquet'} ({len(desc_2d)} rows)")

all_3d = []
patients_3d = list_patient_dirs(root_dir / "data" / "profiles_3D")
for patient in patients_3d:
    f = (
        root_dir
        / "data"
        / "profiles_3D"
        / patient
        / "5.normalized_profiles"
        / "organoid_norm.parquet"
    )
    if not f.exists():
        continue
    df = pd.read_parquet(f)
    df = harmonize_metadata(df, "3D", patient)
    feature_cols = [c for c in df.columns if not c.startswith("Metadata_")]
    agg = describe_group(df, feature_cols, ["Metadata_patient", "Metadata_treatment"])
    agg["modality"] = "3D"
    all_3d.append(agg)
    print(f"3D {patient}: {len(feature_cols)} features, {len(df)} organoids")

desc_3d = pd.concat(all_3d, ignore_index=True)
desc_3d.to_parquet(results_dir / "descriptive_stats_3D.parquet", index=False)
print(f"Wrote {results_dir / 'descriptive_stats_3D.parquet'} ({len(desc_3d)} rows)")
