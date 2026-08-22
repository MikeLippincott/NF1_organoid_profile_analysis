#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
"""Section 1 (cont'd): Kruskal-Wallis test per feature across treatments,
within each patient, FDR-corrected (Benjamini-Hochberg). Uses the same
organoid-level tables as script 1. Saves a significant-feature table and a
per-patient/feature-category count summary for plotting in script 3.
"""


# In[2]:

import sys
import warnings

import pandas as pd
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests

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
    parse_feature_2d,
    parse_feature_3d,
)

results_dir = root_dir / "4.analysis" / "results" / "descriptive_stats"
results_dir.mkdir(parents=True, exist_ok=True)


# In[5]:


def test_patient(df: pd.DataFrame, feature_cols) -> pd.DataFrame:
    rows = []
    groups = df.groupby("Metadata_treatment")
    if len(groups) < 2:
        return pd.DataFrame(rows)
    for feat in feature_cols:
        vals_by_group = [
            pd.to_numeric(g[feat], errors="coerce").dropna().values
            for _, g in groups
            if len(g) > 0
        ]
        vals_by_group = [v for v in vals_by_group if len(v) >= 2]
        if len(vals_by_group) < 2 or all(len(set(v)) == 1 for v in vals_by_group):
            continue
        try:
            stat, p = kruskal(*vals_by_group)
        except ValueError:
            continue
        rows.append({"feature": feat, "statistic": stat, "p_value": p})
    return pd.DataFrame(rows)


# In[6]:


def run_modality(
    modality, patients, load_fn, parse_fn, extra_cols_fn=None, projection=None
):
    all_rows = []
    for patient in patients:
        df = load_fn(patient)
        if df is None:
            continue
        df = harmonize_metadata(df, modality, patient)
        feature_cols = [c for c in df.columns if not c.startswith("Metadata_")]
        res = test_patient(df, feature_cols)
        if res.empty:
            continue
        reject, p_fdr, _, _ = multipletests(res["p_value"], method="fdr_bh", alpha=0.05)
        res["p_fdr"] = p_fdr
        res["significant"] = reject
        res["Metadata_patient"] = patient
        parsed = res["feature"].apply(parse_fn)
        res["compartment"] = parsed.apply(lambda x: x[0])
        res["category"] = parsed.apply(lambda x: x[1])
        all_rows.append(res)
        n_sig = int(res["significant"].sum())
        print(
            f"{modality} {projection or ''} {patient}: {n_sig}/{len(res)} significant features (FDR<0.05)"
        )
    if not all_rows:
        return pd.DataFrame()
    return pd.concat(all_rows, ignore_index=True)


patients_2d = list_patient_dirs(root_dir / "data" / "profiles_2D")
all_2d = []
for projection in PROJECTIONS:
    prefix = PROJECTION_FILE_PREFIX[projection]

    def load_2d(patient, prefix=prefix):
        f = (
            root_dir
            / "data"
            / "profiles_2D"
            / patient
            / "5.normalized"
            / f"{prefix}_organoid.parquet"
        )
        return pd.read_parquet(f) if f.exists() else None

    res = run_modality(
        "2D", patients_2d, load_2d, parse_feature_2d, projection=projection
    )
    if not res.empty:
        res["projection"] = projection
        all_2d.append(res)

sig_2d = pd.concat(all_2d, ignore_index=True) if all_2d else pd.DataFrame()
sig_2d.to_parquet(results_dir / "feature_significance_2D.parquet", index=False)
print(f"Wrote {results_dir / 'feature_significance_2D.parquet'} ({len(sig_2d)} rows)")

patients_3d = list_patient_dirs(root_dir / "data" / "profiles_3D")


# In[7]:


def load_3d(patient):
    f = (
        root_dir
        / "data"
        / "profiles_3D"
        / patient
        / "5.normalized_profiles"
        / "organoid_norm.parquet"
    )
    return pd.read_parquet(f) if f.exists() else None


sig_3d = run_modality("3D", patients_3d, load_3d, parse_feature_3d)
sig_3d.to_parquet(results_dir / "feature_significance_3D.parquet", index=False)
print(f"Wrote {results_dir / 'feature_significance_3D.parquet'} ({len(sig_3d)} rows)")

# Summary counts per patient x category x modality, for the R plotting script
summary_rows = []
for name, df, modality in [("2D", sig_2d, "2D"), ("3D", sig_3d, "3D")]:
    if df.empty:
        continue
    group_cols = ["Metadata_patient", "category"] + (
        ["projection"] if "projection" in df.columns else []
    )
    summ = (
        df.groupby(group_cols)["significant"]
        .agg(n_significant="sum", n_tested="count")
        .reset_index()
    )
    summ["modality"] = modality
    summary_rows.append(summ)

summary = pd.concat(summary_rows, ignore_index=True) if summary_rows else pd.DataFrame()
summary.to_csv(results_dir / "feature_significance_summary.csv", index=False)
print(f"Wrote {results_dir / 'feature_significance_summary.csv'} ({len(summary)} rows)")
