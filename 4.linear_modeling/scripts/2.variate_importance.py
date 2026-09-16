#!/usr/bin/env python
# coding: utf-8

# # Variate importance: how much does each term contribute to each feature?
#
# `0.linear_modeling.ipynb` fits, per (treatment+DMSO combo, feature), the model
#
# ```
# feature ~ C(treatment) * C(patient) + cell_count + organoid_count + cell_per_organoid_count
# ```
#
# and saves per-term **coefficients and p-values** (`1.EDA/results/linear_modeling/*.parquet`).
# That table answers *"is term X significant, and in which direction"* but not
# *"how much of this feature's variance does term X actually account for"* —
# coefficients on different terms aren't on a comparable scale (a dummy
# contrast vs. a per-unit slope), so they can't be ranked against each other
# directly.
#
# This notebook answers that: for every (combo, feature) model, it decomposes
# the model's total sum of squares into the share attributable to each term
# via a **Type II ANOVA**, and expresses each term's contribution as a
# percentage of total variance explained. That is directly comparable across
# terms and features, and lets us rank `treatment`, `patient`,
# `treatment:patient`, `cell_count`, `organoid_count`, and
# `cell_per_organoid_count` by overall importance.
#
# Refits the same ~11k models as `0.linear_modeling.ipynb` (161 organoid +
# 278 single-cell features x 25 treatments), ~10 minutes total.
#

# In[1]:


import pathlib
import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from notebook_init_utils import init_notebook
from statsmodels.stats.anova import anova_lm

warnings.filterwarnings("ignore")
root_dir, in_notebook = init_notebook()

if in_notebook:
    from tqdm.notebook import tqdm
else:
    from tqdm import tqdm

results_path = pathlib.Path(root_dir, "4.linear_modeling/results/variate_importance")
results_path.mkdir(parents=True, exist_ok=True)

pd.set_option("display.max_columns", 50)


# In[2]:


profile_dict = {
    "organoid": pathlib.Path(
        root_dir,
        "data/profiles_3D/all_patients/1.feature_selected_profiles/organoid_norm_fs_profiles.parquet",
    ),
    "single_cell": pathlib.Path(
        root_dir,
        "data/profiles_3D/all_patients/1.feature_selected_profiles/sc_norm_fs_profiles.parquet",
    ),
}
organoid_counts_source = pathlib.Path(
    root_dir,
    "data/profiles_3D/all_patients/1.feature_selected_profiles/organoid_norm_fs_profiles.parquet",
)

# friendly names for the ANOVA table's term labels
TERM_RENAME = {
    "C(Metadata_treatment_full)": "treatment",
    "C(patient)": "patient",
    "C(Metadata_treatment_full):C(patient)": "treatment:patient",
    "cell_count": "cell_count",
    "organoid_count": "organoid_count",
    "cell_per_organoid_count": "cell_per_organoid_count",
    "Residual": "residual",
}


# ## Load + prep, exactly mirroring `0.linear_modeling.ipynb`
#
# Same patient exclusion, same combined treatment+dose label, same count
# covariates, same texture-feature drop, same outlier clipping — so these
# variance-partitioning results line up 1:1 with the existing coefficient
# table.
#

# In[3]:


def load_and_prep(profile_name):
    df = pd.read_parquet(profile_dict[profile_name])
    df = df.rename(
        columns={
            "Metadata_Biology_PatientTumor": "patient",
            "Metadata_Experiment_Treatment": "treatment",
        }
    )
    df = df.loc[df["patient"] != "NF0037_T1_CQ1"]
    df["Metadata_treatment_full"] = (
        df["treatment"].astype(str)
        + "_"
        + df["Metadata_Experiment_Dose"].astype(str)
        + df["Metadata_Experiment_Unit"].astype(str)
    )
    treatment_meta = (
        df[
            [
                "Metadata_treatment_full",
                "treatment",
                "Metadata_Experiment_TherapeuticCategories",
            ]
        ]
        .drop_duplicates()
        .set_index("Metadata_treatment_full")
    )

    if profile_name == "single_cell":
        organoid_counts_df = pd.read_parquet(
            organoid_counts_source,
            columns=[
                "Metadata_Biology_PatientTumor",
                "Metadata_Experiment_Well",
                "Metadata_WellOrganoidCount",
            ],
        ).drop_duplicates()
        organoid_counts_df = organoid_counts_df.rename(
            columns={"Metadata_Biology_PatientTumor": "patient"}
        )
        df = df.merge(
            organoid_counts_df, on=["patient", "Metadata_Experiment_Well"], how="left"
        )
        df = df.dropna(subset=["Metadata_WellOrganoidCount"])
        df["cell_count"] = df["Metadata_Object_WellSingleCellCount"]
        df["organoid_count"] = df["Metadata_WellOrganoidCount"]
    else:
        df["cell_count"] = df["Metadata_Object_OrganoidSingleCellCount"]
        df["organoid_count"] = df["Metadata_WellOrganoidCount"]
    df["cell_per_organoid_count"] = df["cell_count"] / df["organoid_count"]

    count_columns = ["cell_count", "organoid_count", "cell_per_organoid_count"]
    metadata_columns = (
        ["patient", "treatment"]
        + count_columns
        + [c for c in df.columns if c.startswith("Metadata_")]
    )
    df = df.drop(columns=[c for c in df.columns if "_Texture_" in c])
    feature_columns = [c for c in df.columns if c not in metadata_columns]
    df[feature_columns] = df[feature_columns].clip(lower=-1e1, upper=1e1)
    for col in df.columns:
        new_col = col.replace(".", "")
        df.rename(columns={col: new_col}, inplace=True)
    feature_columns = [c.replace(".", "") for c in feature_columns]

    dmso_label = df.loc[df["treatment"] == "DMSO", "Metadata_treatment_full"].unique()[
        0
    ]
    combo_list = [
        (dmso_label, t)
        for t in df["Metadata_treatment_full"].unique()
        if t != dmso_label
    ]
    return df, feature_columns, treatment_meta, combo_list


# In[4]:


def variance_partition(df, feature_columns, treatment_meta, combo_list, profile_name):
    rows = {
        "profile": [],
        "treatment": [],
        "drug": [],
        "therapeutic_category": [],
        "feature": [],
        "term": [],
        "sum_sq": [],
        "df": [],
        "fvalue": [],
        "pvalue": [],
        "pct_variance_explained": [],
        "rsquared": [],
    }
    for combo in tqdm(
        combo_list, desc=f"{profile_name}: treatment combos", unit="combo"
    ):
        drug_name = treatment_meta.loc[combo[1], "treatment"]
        therapeutic_category = treatment_meta.loc[
            combo[1], "Metadata_Experiment_TherapeuticCategories"
        ]

        df_trt = df.loc[df["Metadata_treatment_full"].isin(combo)].copy()
        df_trt["Metadata_treatment_full"] = pd.Categorical(
            df_trt["Metadata_treatment_full"], categories=list(combo)
        )
        patients_in_combo = sorted(df_trt["patient"].unique())
        df_trt["patient"] = pd.Categorical(
            df_trt["patient"], categories=patients_in_combo
        )

        for col in tqdm(feature_columns, desc="features", unit="feature", leave=False):
            formula = (
                f"Q('{col}') ~ C(Metadata_treatment_full) * C(patient)"
                " + cell_count + organoid_count + cell_per_organoid_count"
            )
            model = smf.ols(formula=formula, data=df_trt)
            results = model.fit()
            aov = anova_lm(results, typ=2)
            total_ss = aov["sum_sq"].sum()

            for term_label, term_row in aov.iterrows():
                term = TERM_RENAME.get(term_label, term_label)
                rows["profile"].append(profile_name)
                rows["treatment"].append(combo[1])
                rows["drug"].append(drug_name)
                rows["therapeutic_category"].append(therapeutic_category)
                rows["feature"].append(col)
                rows["term"].append(term)
                rows["sum_sq"].append(term_row["sum_sq"])
                rows["df"].append(term_row["df"])
                rows["fvalue"].append(term_row.get("F", np.nan))
                rows["pvalue"].append(term_row.get("PR(>F)", np.nan))
                rows["pct_variance_explained"].append(
                    term_row["sum_sq"] / total_ss * 100
                )
                rows["rsquared"].append(results.rsquared)
    return pd.DataFrame(rows)


# In[5]:


def split_feature_name(pdf):
    pdf = pdf.copy()
    pdf[["Compartment", "Channel", "Feature_type", "Measurement"]] = pdf[
        "feature"
    ].str.split("_", n=3, expand=True)
    pdf.loc[pdf["Feature_type"] == "AreaSizeShape", "Measurement"] = pdf["Channel"]
    pdf.loc[pdf["Feature_type"] == "AreaSizeShape", "Channel"] = None
    return pdf


partition_dfs = {}
for profile_name in profile_dict:
    out_path = results_path / f"{profile_name}_variance_partition.parquet"
    if out_path.exists():
        print(f"{profile_name}: loading cached fit from {out_path}")
        pdf = pd.read_parquet(out_path)
    else:
        df, feature_columns, treatment_meta, combo_list = load_and_prep(profile_name)
        pdf = variance_partition(
            df, feature_columns, treatment_meta, combo_list, profile_name
        )
        pdf.to_parquet(out_path, index=False)
        print(f"{profile_name}: {pdf.shape} -> {out_path}")
    partition_dfs[profile_name] = split_feature_name(pdf)

variance_df = pd.concat(partition_dfs.values(), ignore_index=True)
print(variance_df.head())


# ## Which variates matter most, overall?
#
# Average each term's share of variance explained across every (combo,
# feature) model fit. This is the headline ranking: on average, how much of a
# feature's variance does each term of the model account for?
#

# In[6]:


TERM_ORDER = [
    "treatment",
    "patient",
    "treatment:patient",
    "cell_count",
    "organoid_count",
    "cell_per_organoid_count",
    "residual",
]

term_summary_dfs = {}
for name, pdf in partition_dfs.items():
    overall = (
        pdf.groupby("term")["pct_variance_explained"]
        .agg(["mean", "median", "std"])
        .reindex(TERM_ORDER)
    )
    term_summary_dfs[name] = overall
    out_path = results_path / f"{name}_term_summary.parquet"
    overall.reset_index().rename(columns={"index": "term"}).to_parquet(
        out_path, index=False
    )
    print(
        f"\n==== {name}: mean % variance explained per term (across {pdf['feature'].nunique()} features x {pdf['treatment'].nunique()} treatments) ===="
    )
    print(overall)
    print(f"saved -> {out_path}")


# ## Per-feature breakdown: how does each variate contribute to *every* feature?
#
# Average each term's % variance explained across all 24 treatment combos,
# per feature — a `feature x term` importance matrix. Saved in full; shown
# here for the features where the model explains the most total variance
# (highest R^2), since a low-R^2 feature's term breakdown is mostly noise.
#

# In[7]:


def feature_term_matrix(pdf):
    matrix = pdf.pivot_table(
        index="feature", columns="term", values="pct_variance_explained", aggfunc="mean"
    )
    matrix = matrix.reindex(columns=TERM_ORDER)
    mean_rsq = pdf.groupby("feature")["rsquared"].mean()
    matrix["mean_rsquared"] = mean_rsq
    return matrix.sort_values("mean_rsquared", ascending=False)


feature_matrices = {}
for name, pdf in partition_dfs.items():
    matrix = feature_term_matrix(pdf)
    feature_matrices[name] = matrix
    out_path = results_path / f"{name}_feature_by_term_importance.parquet"
    matrix.reset_index().to_parquet(out_path, index=False)
    print(
        f"{name}: saved {matrix.shape[0]} features x {matrix.shape[1]} columns -> {out_path}"
    )
    print(matrix.head(15))


# ## Which channel / feature type is most driven by treatment vs. by nuisance covariates?
#
# Roll the per-(combo, feature) contributions up to `Channel` x `Feature_type`
# and compare the average share attributed to `treatment` (+ its patient
# interaction) against the share attributed to the three count covariates.
# Feature types dominated by counts are the ones a density-normalization
# would matter most for.
#

# In[8]:


def channel_feature_term_summary(pdf):
    pdf = pdf.copy()
    pdf["term_group"] = pdf["term"].map(
        {
            "treatment": "treatment (direct + response)",
            "treatment:patient": "treatment (direct + response)",
            "patient": "patient baseline",
            "cell_count": "count covariates",
            "organoid_count": "count covariates",
            "cell_per_organoid_count": "count covariates",
            "residual": "residual",
        }
    )
    return (
        pdf.groupby(["Channel", "Feature_type", "term_group"], dropna=False)[
            "pct_variance_explained"
        ]
        .mean()
        .reset_index()
        .pivot_table(
            index=["Channel", "Feature_type"],
            columns="term_group",
            values="pct_variance_explained",
        )
    )


channel_feature_summaries = {}
for name, pdf in partition_dfs.items():
    summary = channel_feature_term_summary(pdf)
    channel_feature_summaries[name] = summary
    out_path = results_path / f"{name}_channel_feature_term_summary.parquet"
    summary.reset_index().to_parquet(out_path, index=False)
    print(f"\n==== {name}: variance share by channel/feature-type ====")
    print(
        summary.sort_values("treatment (direct + response)", ascending=False).head(15)
    )
    print(f"saved -> {out_path}")


# ## Organoid vs. single-cell: does variate importance agree across scale?
#

# In[9]:


compare = pd.concat(
    [
        pdf.groupby("term")["pct_variance_explained"].mean().rename(name)
        for name, pdf in partition_dfs.items()
    ],
    axis=1,
).reindex(TERM_ORDER)
compare["difference (organoid - single_cell)"] = (
    compare["organoid"] - compare["single_cell"]
)
out_path = results_path / "term_importance_organoid_vs_sc.parquet"
compare.reset_index().rename(columns={"index": "term"}).to_parquet(
    out_path, index=False
)
print(compare)
print(f"saved -> {out_path}")


# ## Summary
#
# - **Ranking terms overall**: see the term-summary table above (also saved
#   to `{organoid,single_cell}_term_summary.parquet`) — plotted in
#   `plot_variate_importance_term_summary.ipynb`.
# - **Per-feature answer**: `4.linear_modeling/results/variate_importance/
#   {organoid,single_cell}_feature_by_term_importance.parquet` has, for every
#   feature, the average % of its variance attributable to `treatment`,
#   `patient`, `treatment:patient`, `cell_count`, `organoid_count`, and
#   `cell_per_organoid_count` — plotted in
#   `plot_variate_importance_feature_heatmap.ipynb`.
# - **Channel/feature-type roll-up** (saved to
#   `{organoid,single_cell}_channel_feature_term_summary.parquet`) flags which
#   measurement families are mostly density-driven vs. treatment-driven —
#   pairs with the density confound finding in `explore_leads.ipynb`.
# - **Organoid vs. single-cell comparison** saved to
#   `term_importance_organoid_vs_sc.parquet`, plotted in
#   `plot_variate_importance_organoid_vs_sc.ipynb`.
# - A large `treatment:patient` share relative to `treatment` on its own is
#   itself a finding — it means response heterogeneity across patients is a
#   bigger driver of that feature than the average drug effect, reinforcing
#   the patient-specific-responder leads from `explore_leads.ipynb`.
#
