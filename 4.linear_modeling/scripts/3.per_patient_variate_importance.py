#!/usr/bin/env python
# coding: utf-8

# # Per-patient variate importance
#
# `1.variate_importance.ipynb` fit one pooled model per (treatment+DMSO combo,
# feature) across all patients:
#
# ```
# feature ~ C(treatment) * C(patient) + cell_count + organoid_count + cell_per_organoid_count
# ```
#
# and found `patient` alone explains far more variance (~7.4%) than the direct
# `treatment` effect (~0.2%) — most of the "treatment" signal was showing up
# in the `treatment:patient` interaction instead, i.e. drug response is
# patient-specific.
#
# This notebook removes `patient` (and therefore `treatment:patient`, which is
# degenerate with a single patient) and fits **one model per patient**, using
# that patient's full treatment panel (all drugs + DMSO, not just one combo at
# a time — with `patient` fixed there's no need to pool patients pairwise
# against DMSO anymore):
#
# ```
# feature ~ C(treatment) + cell_count + organoid_count + cell_per_organoid_count
# ```
#
# Same variates as before, patient removed, as requested. This isolates each
# patient's own signal and checks whether `treatment`'s importance recovers
# once between-patient noise is no longer competing with it in the same
# model.
#

# In[1]:


import pathlib
import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from notebook_init_utils import init_notebook
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")
root_dir, in_notebook = init_notebook()

if in_notebook:
    from tqdm.notebook import tqdm
else:
    from tqdm import tqdm

results_path = pathlib.Path(
    root_dir, "4.linear_modeling/results/per_patient_variate_importance"
)
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

TERM_RENAME = {
    "C(Metadata_treatment_full)": "treatment",
    "cell_count": "cell_count",
    "organoid_count": "organoid_count",
    "cell_per_organoid_count": "cell_per_organoid_count",
    "Residual": "residual",
}
TERM_ORDER = [
    "treatment",
    "cell_count",
    "organoid_count",
    "cell_per_organoid_count",
    "residual",
]


# ## Load + prep
#
# Identical to `4.variate_importance.ipynb` up through building the count
# covariates and the combined treatment+dose label — the only change is
# downstream: we don't build DMSO-vs-single-drug combos anymore, since each
# patient's model uses their whole treatment panel at once.
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
    return df, feature_columns, treatment_meta, dmso_label


# ## Fit one model per (patient, feature)
#
# For each patient, `treatment` is releveled so DMSO is the reference; the
# model is `feature ~ C(treatment) + cell_count + organoid_count +
# cell_per_organoid_count`, fit only on that patient's rows. We collect both:
#
# - the **Type II ANOVA** decomposition (`treatment` here is one combined term
#   across every drug the patient received, vs. the three count covariates)
# - the **per-drug coefficient/p-value** (drug vs. DMSO, within that
#   patient), analogous to the `term == "treatment"` rows in
#   `1.EDA/results/linear_modeling/*.parquet`, but now a direct fit instead
#   of an interaction contrast.
#

# In[4]:


def per_patient_variate_importance(
    df, feature_columns, treatment_meta, dmso_label, profile_name
):
    anova_rows = {
        "profile": [],
        "patient": [],
        "feature": [],
        "term": [],
        "sum_sq": [],
        "df": [],
        "fvalue": [],
        "pvalue": [],
        "pct_variance_explained": [],
        "rsquared": [],
    }
    coef_rows = {
        "profile": [],
        "patient": [],
        "drug": [],
        "treatment": [],
        "therapeutic_category": [],
        "feature": [],
        "coefficient": [],
        "pvalue": [],
        "rsquared": [],
    }

    patients = sorted(df["patient"].unique())
    for patient in tqdm(patients, desc=f"{profile_name}: patients", unit="patient"):
        df_pat = df.loc[df["patient"] == patient].copy()
        treatments_present = [
            t for t in df_pat["Metadata_treatment_full"].unique() if t != dmso_label
        ]
        df_pat["Metadata_treatment_full"] = pd.Categorical(
            df_pat["Metadata_treatment_full"],
            categories=[dmso_label] + treatments_present,
        )

        for col in tqdm(feature_columns, desc="features", unit="feature", leave=False):
            formula = (
                f"Q('{col}') ~ C(Metadata_treatment_full)"
                " + cell_count + organoid_count + cell_per_organoid_count"
            )
            model = smf.ols(formula=formula, data=df_pat)
            results = model.fit()
            aov = anova_lm(results, typ=2)
            total_ss = aov["sum_sq"].sum()

            for term_label, term_row in aov.iterrows():
                term = TERM_RENAME.get(term_label, term_label)
                anova_rows["profile"].append(profile_name)
                anova_rows["patient"].append(patient)
                anova_rows["feature"].append(col)
                anova_rows["term"].append(term)
                anova_rows["sum_sq"].append(term_row["sum_sq"])
                anova_rows["df"].append(term_row["df"])
                anova_rows["fvalue"].append(term_row.get("F", np.nan))
                anova_rows["pvalue"].append(term_row.get("PR(>F)", np.nan))
                anova_rows["pct_variance_explained"].append(
                    term_row["sum_sq"] / total_ss * 100
                )
                anova_rows["rsquared"].append(results.rsquared)

            for trt in treatments_present:
                term_key = f"C(Metadata_treatment_full)[T.{trt}]"
                if term_key not in results.params:
                    continue
                coef_rows["profile"].append(profile_name)
                coef_rows["patient"].append(patient)
                coef_rows["drug"].append(treatment_meta.loc[trt, "treatment"])
                coef_rows["treatment"].append(trt)
                coef_rows["therapeutic_category"].append(
                    treatment_meta.loc[trt, "Metadata_Experiment_TherapeuticCategories"]
                )
                coef_rows["feature"].append(col)
                coef_rows["coefficient"].append(results.params[term_key].item())
                coef_rows["pvalue"].append(results.pvalues[term_key].item())
                coef_rows["rsquared"].append(results.rsquared)

    return pd.DataFrame(anova_rows), pd.DataFrame(coef_rows)


# In[5]:


def split_feature_name(pdf):
    pdf = pdf.copy()
    pdf[["Compartment", "Channel", "Feature_type", "Measurement"]] = pdf[
        "feature"
    ].str.split("_", n=3, expand=True)
    pdf.loc[pdf["Feature_type"] == "AreaSizeShape", "Measurement"] = pdf["Channel"]
    pdf.loc[pdf["Feature_type"] == "AreaSizeShape", "Channel"] = None
    return pdf


anova_dfs = {}
coef_dfs = {}
for profile_name in profile_dict:
    anova_out_path = results_path / f"{profile_name}_per_patient_anova.parquet"
    coef_out_path = results_path / f"{profile_name}_per_patient_coefficients.parquet"
    if anova_out_path.exists() and coef_out_path.exists():
        print(f"{profile_name}: loading cached fit")
        anova_pdf = pd.read_parquet(anova_out_path)
        coef_pdf = pd.read_parquet(coef_out_path)
    else:
        df, feature_columns, treatment_meta, dmso_label = load_and_prep(profile_name)
        anova_pdf, coef_pdf = per_patient_variate_importance(
            df, feature_columns, treatment_meta, dmso_label, profile_name
        )
        # FDR-correct the per-drug coefficient p-values, per patient (own testing family)
        coef_pdf["pvalue_fdr"] = float("nan")
        for patient, idx in coef_pdf.groupby("patient").groups.items():
            _, fdr, _, _ = multipletests(
                coef_pdf.loc[idx, "pvalue"].values, method="fdr_bh"
            )
            coef_pdf.loc[idx, "pvalue_fdr"] = fdr
        anova_pdf.to_parquet(anova_out_path, index=False)
        coef_pdf.to_parquet(coef_out_path, index=False)
        print(f"{profile_name}: anova {anova_pdf.shape} -> {anova_out_path}")
        print(f"{profile_name}: coefficients {coef_pdf.shape} -> {coef_out_path}")
    anova_dfs[profile_name] = split_feature_name(anova_pdf)
    coef_dfs[profile_name] = coef_pdf


# ## Which variates matter most, per patient, once `patient` is removed?
#
# Compare this directly against the pooled-model numbers from
# `4.variate_importance.ipynb` (`treatment` ~0.2%, `treatment:patient` ~1.5%,
# `patient` ~7.4%, mostly `residual`). If treatment's importance recovers
# once patient variance is no longer absorbing it, that confirms the pooled
# model's `patient`/`treatment:patient` split was masking a real per-patient
# treatment effect rather than the effect actually being negligible.
#

# In[6]:


TERM_ORDER = [
    "treatment",
    "cell_count",
    "organoid_count",
    "cell_per_organoid_count",
    "residual",
]

term_summary_dfs = {}
for name, pdf in anova_dfs.items():
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
        f"\n==== {name}: mean % variance explained per term, pooled across all {pdf['patient'].nunique()} patients ===="
    )
    print(overall)
    print(f"saved -> {out_path}")


# ## Per-(patient, feature) importance table
#
# Same shape as the per-feature CSV from `4.variate_importance.ipynb`, but
# now with `patient` as an extra axis: for every (patient, feature), the %
# variance explained by `treatment`, `cell_count`, `organoid_count`, and
# `cell_per_organoid_count`. This is the full answer to "how each variate
# contributes to every feature," now resolved per patient rather than
# averaged over the whole cohort.
#

# In[7]:


for name, pdf in anova_dfs.items():
    matrix = pdf.pivot_table(
        index=["patient", "feature"],
        columns="term",
        values="pct_variance_explained",
        aggfunc="mean",
    ).reindex(columns=TERM_ORDER)
    mean_rsq = pdf.groupby(["patient", "feature"])["rsquared"].mean()
    matrix["mean_rsquared"] = mean_rsq
    matrix = matrix.sort_values("mean_rsquared", ascending=False)
    out_path = results_path / f"{name}_patient_by_feature_by_term_importance.parquet"
    matrix.reset_index().to_parquet(out_path, index=False)
    print(f"{name}: saved {matrix.shape[0]} (patient, feature) rows -> {out_path}")
    print(matrix.head(15))


# ## Most important variate per (patient, feature) model
#
# Naively taking `argmax` over each model's four non-residual terms is biased:
# `treatment` is a single combined term spanning every drug the patient
# received (df = n_drugs - 1), so it mechanically accumulates more sum of
# squares than a single-df covariate like `cell_count`, even when it isn't
# truly more informative per parameter. We report both:
#
# - **raw winner** — which term has the single largest % variance explained
#   (what "most important variate" naively means)
# - **df-normalized winner** — % variance explained divided by the term's
#   degrees of freedom, i.e. average explanatory power *per parameter*, which
#   is the fairer comparison between a multi-level factor and a single slope
#
# Both are restricted to well-explained models (`mean_rsquared > 0.3`) as the
# headline number, since a "winner" from a model that explains 5% of a
# feature's variance overall isn't a meaningful lead — but the unrestricted
# counts are reported too for completeness.
#

# In[8]:


NON_RESIDUAL_TERMS = [t for t in TERM_ORDER if t != "residual"]


def safe_idxmax(frame):
    # frame.idxmax(axis=1) raises ValueError if any row is all-NaN (e.g. a
    # term was dropped from a given model, such as a collinear covariate);
    # skip those rows instead of failing the whole column.
    valid_rows = frame.notna().any(axis=1)
    result = pd.Series(np.nan, index=frame.index, dtype=object)
    result.loc[valid_rows] = frame.loc[valid_rows].idxmax(axis=1)
    return result


top_variate_dfs = {}
for name, pdf in anova_dfs.items():
    wide_pct = pdf.pivot_table(
        index=["patient", "feature"], columns="term", values="pct_variance_explained"
    )
    wide_df = pdf.pivot_table(index=["patient", "feature"], columns="term", values="df")
    wide_rsq = pdf.groupby(["patient", "feature"])["rsquared"].mean()

    per_df = wide_pct[NON_RESIDUAL_TERMS] / wide_df[NON_RESIDUAL_TERMS]

    out = pd.DataFrame(index=wide_pct.index)
    out["mean_rsquared"] = wide_rsq
    out["raw_top_variate"] = safe_idxmax(wide_pct[NON_RESIDUAL_TERMS])
    out["raw_top_variate_pct"] = wide_pct[NON_RESIDUAL_TERMS].max(axis=1)
    out["df_normalized_top_variate"] = safe_idxmax(per_df)
    out["df_normalized_top_variate_pct_per_df"] = per_df.max(axis=1)
    for term in NON_RESIDUAL_TERMS:
        out[f"{term}_pct_variance"] = wide_pct[term]
    out = out.reset_index().sort_values("mean_rsquared", ascending=False)

    out_path = results_path / f"{name}_top_variate_per_patient_feature.parquet"
    out.to_parquet(out_path, index=False)
    top_variate_dfs[name] = out
    print(f"{name}: saved {len(out)} (patient, feature) rows -> {out_path}")


# In[9]:


for name, out in top_variate_dfs.items():
    well = out[out["mean_rsquared"] > 0.3]
    print(
        f"\n==== {name}: n={len(out)} total, n={len(well)} well-explained (mean R^2 > 0.3) ===="
    )
    print("raw winner, well-explained models only:")
    print(
        (well["raw_top_variate"].value_counts(normalize=True) * 100)
        .round(1)
        .rename("% of models")
    )
    print("df-normalized winner, well-explained models only (the fairer comparison):")
    print(
        (well["df_normalized_top_variate"].value_counts(normalize=True) * 100)
        .round(1)
        .rename("% of models")
    )


# ## Strongest per-patient drug effects
#
# From the per-drug coefficient table: the largest-magnitude, FDR-significant
# drug-vs-DMSO effects fit directly within each patient (no interaction
# contrast involved).
#

# In[10]:


for name, coef_pdf in coef_dfs.items():
    print(f"\n==== {name}: strongest per-patient drug effects (FDR < 0.05) ====")
    sig = coef_pdf[coef_pdf["pvalue_fdr"] < 0.05].copy()
    sig["abs_coefficient"] = sig["coefficient"].abs()
    print(
        sig.sort_values("abs_coefficient", ascending=False)[
            [
                "patient",
                "drug",
                "therapeutic_category",
                "feature",
                "coefficient",
                "pvalue_fdr",
                "rsquared",
            ]
        ].head(15)
    )


# ## Summary
#
# - Removing `patient` and refitting **one model per patient** (`feature ~
#   C(treatment) + cell_count + organoid_count + cell_per_organoid_count`)
#   isolates each patient's own treatment signal instead of pooling it
#   through a treatment:patient interaction contrast.
# - Compare the "mean % variance explained per term" table above against
#   `1.variate_importance.ipynb`'s pooled numbers (treatment ~0.2%,
#   treatment:patient ~1.5%, patient ~7.4%) — plotted side by side in
#   `plot_per_patient_term_summary.ipynb`.
# - Per-patient treatment importance and the df-normalized winning-variate
#   breakdown are saved and plotted in `plot_per_patient_treatment_importance_by_patient.ipynb`
#   and `plot_per_patient_top_variate_by_patient.ipynb`.
# - Full per-(patient, feature) term breakdown:
#   `4.linear_modeling/results/per_patient_variate_importance/
#   {organoid,single_cell}_patient_by_feature_by_term_importance.parquet`.
# - Per-drug, per-patient significant effects (no pooling/interaction
#   assumptions) are saved in `{organoid,single_cell}_per_patient_
#   coefficients.parquet`, visualized in `plot_per_patient_coefficient_overview.ipynb`,
#   `plot_per_patient_volcano.ipynb`, and `plot_per_patient_treatment_heatmap.ipynb`.
# - See `significant_variates_per_patient_feature.ipynb` for which of these
#   terms are *statistically significant* (not just largest) per model, with
#   interpretation.
#
