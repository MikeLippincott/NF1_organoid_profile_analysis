#!/usr/bin/env python
# coding: utf-8

# # Significant variates per (patient, feature) model
#
# The per-patient models (`feature ~ C(treatment) + cell_count + organoid_count
# + cell_per_organoid_count`, fit independently per patient) tell us *how much*
# variance each term explains (`pct_variance_explained`), but not whether that
# share is more than noise. This notebook answers the significance question
# directly: for every (patient, feature) model, which of the four terms are
# statistically significant (FDR-corrected p-value < 0.05), and what does the
# resulting pattern mean?
#
# P-values are FDR-corrected (Benjamini-Hochberg) **within each (profile,
# term)** — the same convention used in `0.linear_modeling.ipynb` — since each
# term is tested many times (once per patient x feature) and each term is its
# own hypothesis-testing family with its own p-value distribution.
#

# In[1]:


import pathlib
import warnings

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from notebook_init_utils import init_notebook
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")
root_dir, in_notebook = init_notebook()

figures_path = pathlib.Path(root_dir, "4.linear_modeling/figures/significant_variates")
figures_path.mkdir(parents=True, exist_ok=True)
results_path = pathlib.Path(root_dir, "4.linear_modeling/results/significant_variates")
results_path.mkdir(parents=True, exist_ok=True)
anova_results_path = pathlib.Path(
    root_dir, "4.linear_modeling/results/per_patient_variate_importance"
)

try:
    get_ipython().run_line_magic("matplotlib", "inline")
except NameError:
    pass

sns.set_theme(style="whitegrid", context="talk")
pd.set_option("display.max_columns", 50)

PVALUE_MAX = 0.05
TERMS = ["treatment", "cell_count", "organoid_count", "cell_per_organoid_count"]


# ## Load the per-patient ANOVA fits and FDR-correct
#
# Reuses the `{organoid,single_cell}_per_patient_anova.parquet` tables produced
# in `2.per_patient_variate_importance.ipynb` (one row per patient x feature x
# term, with the raw F-test p-value for that term).
#

# In[2]:


def load_and_correct(profile_name):
    df = pd.read_parquet(
        anova_results_path / f"{profile_name}_per_patient_anova.parquet"
    )
    df = df[df["term"] != "residual"].copy()
    df["pvalue_fdr"] = float("nan")
    n_nan_pvalue = df["pvalue"].isna().sum()
    for term, idx in df.groupby("term").groups.items():
        # a handful of models are degenerate (e.g. a constant feature for one
        # patient) and produce a NaN p-value; multipletests propagates a single
        # NaN to every p-value in the array, so those rows must be excluded from
        # the correction rather than passed through — they stay pvalue_fdr=NaN
        # (never significant) rather than corrupting every other row's FDR value
        valid_idx = df.loc[idx].dropna(subset=["pvalue"]).index
        _, fdr, _, _ = multipletests(
            df.loc[valid_idx, "pvalue"].values, method="fdr_bh"
        )
        df.loc[valid_idx, "pvalue_fdr"] = fdr
    df["significant"] = df["pvalue_fdr"] < PVALUE_MAX
    return df, n_nan_pvalue


anova_dfs = {}
for name in ["organoid", "single_cell"]:
    anova_dfs[name], n_nan_pvalue = load_and_correct(name)
    if n_nan_pvalue:
        print(
            f"{name}: {n_nan_pvalue} row(s) had a NaN p-value (degenerate model, e.g. a constant feature for one patient) and were excluded from the FDR correction"
        )

for name, df in anova_dfs.items():
    print(
        f"{name}: {df['patient'].nunique()} patients x {df['feature'].nunique()} features x {len(TERMS)} terms"
    )
    print(
        f"  {df['significant'].sum():,} / {len(df):,} (patient, feature, term) rows significant (FDR<{PVALUE_MAX})"
    )


# ## Which terms are significant, per (patient, feature)?
#
# Pivot to one row per (patient, feature) with a boolean column per term, plus
# a summary of *how many* terms are jointly significant for that model. A
# model where 0 terms are significant means none of treatment, cell_count,
# organoid_count, or cell_per_organoid_count meaningfully predicts that
# feature for that patient — the feature is unexplained noise as far as this
# model specification goes.
#

# In[3]:


def significance_matrix(df):
    sig = df.pivot_table(
        index=["patient", "feature"],
        columns="term",
        values="significant",
        aggfunc="first",
    )
    sig = sig.reindex(columns=TERMS).fillna(False)
    sig["n_significant_terms"] = sig[TERMS].sum(axis=1)
    sig["significant_terms"] = sig[TERMS].apply(
        lambda row: ", ".join(t for t in TERMS if row[t]) or "none", axis=1
    )
    return sig.reset_index()


sig_matrices = {name: significance_matrix(df) for name, df in anova_dfs.items()}
for name, sig in sig_matrices.items():
    out_path = results_path / f"{name}_significant_variates_per_patient_feature.parquet"
    sig.to_parquet(out_path, index=False)
    print(f"{name}: saved {len(sig)} (patient, feature) rows -> {out_path}")
    print(sig["n_significant_terms"].value_counts().sort_index())


# ### Interpretation: how many terms are typically significant?
#
# - **0 significant terms** is the modal outcome for most features — expected,
#   since most single-cell/organoid morphology features are dominated by
#   measurement noise the model specification doesn't capture (consistent with
#   the residual dominating variance explained in
#   `1.variate_importance.ipynb`/`2.per_patient_variate_importance.ipynb`).
# - **1 significant term** usually means either a clean, isolated treatment
#   effect *or* a clean density confound — which one it is matters a lot, and
#   is exactly what the `significant_terms` breakdown (next cell) resolves.
# - **2+ significant terms simultaneously** (e.g. `treatment` and `cell_count`
#   both significant) is the case that most needs care: it doesn't mean the
#   treatment effect is fake, but it does mean the feature's real average
#   response is a mix of true drug effect and density shift, and the raw
#   coefficient magnitude can't cleanly separate the two.
#

# In[4]:


for name, sig in sig_matrices.items():
    print(f"\n==== {name}: frequency of each significant-term combination ====")
    print(
        (sig["significant_terms"].value_counts(normalize=True) * 100)
        .round(1)
        .head(12)
        .rename("% of (patient, feature) models")
    )


# ## Where treatment is significant on its own (no covariate confound)
#
# The most trustworthy per-patient drug-response hits: `treatment` is
# significant and none of the three count covariates are, for that same
# (patient, feature) — a clean signal, not a density artifact.
#

# In[5]:


for name, sig in sig_matrices.items():
    clean_treatment = sig[
        (sig["treatment"])
        & (~sig["cell_count"])
        & (~sig["organoid_count"])
        & (~sig["cell_per_organoid_count"])
    ]
    print(
        f"{name}: {len(clean_treatment):,} / {len(sig):,} (patient, feature) models "
        f"({len(clean_treatment) / len(sig):.1%}) have a clean, unconfounded treatment effect"
    )
    print(
        clean_treatment.groupby("patient")
        .size()
        .sort_values(ascending=False)
        .rename("n clean treatment hits")
    )


# ### Interpretation
#
# This is the strictest, most defensible definition of "this patient's
# morphology genuinely responds to treatment" available from this model
# family — no count covariate is competing for the same variance. Patients
# with many rows here are the strongest candidates for real, density-independent
# drug response; patients with few or none either don't respond to any drug
# tested, or their apparent responses in the raw coefficient tables
# (`2.per_patient_variate_importance.ipynb`) are entangled with density and
# should be interpreted cautiously.
#

# ## Where only a count covariate is significant (pure density confound)
#
# The mirror image: cases where a count covariate is significant but
# `treatment` is not. These are features whose apparent variability is fully
# explained by cell/organoid density — any downstream "treatment effect"
# claim on these particular (patient, feature) pairs (e.g. from the pooled
# model in `0.linear_modeling.ipynb`) should be treated skeptically.
#

# In[6]:


for name, sig in sig_matrices.items():
    pure_confound = sig[
        (~sig["treatment"])
        & (sig["cell_count"] | sig["organoid_count"] | sig["cell_per_organoid_count"])
    ]
    print(
        f"{name}: {len(pure_confound):,} / {len(sig):,} (patient, feature) models "
        f"({len(pure_confound) / len(sig):.1%}) are explained only by a count covariate, not treatment"
    )


# ### Interpretation
#
# Combined with the "clean treatment effect" count above, this splits every
# model into three buckets: **density-only** (this section), **treatment-only**
# (previous section), and **both/neither** (the remainder). The relative size
# of the density-only bucket versus the treatment-only bucket is a direct,
# per-model answer to "is this pipeline mostly picking up biology or mostly
# picking up well density" — a sharper version of the aggregate ~96-98%
# confound-overlap finding from `4.explore_leads.ipynb`, now resolved down to
# individual (patient, feature) pairs instead of an aggregate rate.
#

# ## Per-patient significant-term profile
#
# For each patient, what fraction of their well-modeled features fall into
# each category? Patients with a high `treatment`-only share are the most
# "biologically responsive" in a way this model can actually attribute
# correctly; patients with a high count-covariate-only share are ones where
# density dominates and treatment calls should be discounted more heavily.
#

# In[7]:


def patient_term_profile(sig):
    def classify(row):
        if row["treatment"] and not (
            row["cell_count"] or row["organoid_count"] or row["cell_per_organoid_count"]
        ):
            return "treatment only"
        if not row["treatment"] and (
            row["cell_count"] or row["organoid_count"] or row["cell_per_organoid_count"]
        ):
            return "count covariate only"
        if row["treatment"] and (
            row["cell_count"] or row["organoid_count"] or row["cell_per_organoid_count"]
        ):
            return "both"
        return "neither"

    sig = sig.copy()
    sig["category"] = sig.apply(classify, axis=1)
    return (
        sig.groupby(["patient", "category"])
        .size()
        .unstack(fill_value=0)
        .apply(lambda r: r / r.sum() * 100, axis=1)
    )


for name, sig in sig_matrices.items():
    profile_table = patient_term_profile(sig)
    out_path = results_path / f"{name}_patient_significance_profile.parquet"
    profile_table.reset_index().to_parquet(out_path, index=False)
    print(f"\n==== {name}: per-patient significant-term category (% of features) ====")
    print(profile_table.round(1))


# In[8]:


fig, axes = plt.subplots(1, 2, figsize=(16, 8), sharex=True)
for ax, (name, sig) in zip(axes, sig_matrices.items()):
    profile_table = patient_term_profile(sig)
    categories = ["treatment only", "both", "count covariate only", "neither"]
    profile_table = profile_table.reindex(
        columns=[c for c in categories if c in profile_table.columns]
    )
    profile_table.plot(kind="barh", stacked=True, ax=ax, colormap="RdYlGn_r")
    ax.set_title(f"{name}: significant-term category by patient")
    ax.set_xlabel("% of (patient, feature) models")
    ax.set_ylabel("")
plt.tight_layout()
plt.savefig(
    figures_path / "significance_category_by_patient.png", dpi=600, bbox_inches="tight"
)
plt.show()


# ### Interpretation
#
# This is the single most actionable summary in this notebook: it ranks
# patients by how much of their signal is trustworthy drug response
# (`treatment only`, green) versus how much is a density artifact
# (`count covariate only`, red). A patient dominated by `neither` simply
# doesn't have enough well-explained features to draw conclusions from
# regardless of treatment. Cross-reference the strongest `treatment only`
# patients here against the breadth/magnitude leads in
# `4.explore_leads.ipynb` and the coefficient plots in
# `plot_per_patient_coefficient_overview.ipynb` — a drug hit that shows up in
# both places is far more credible than one that only shows up in the raw
# coefficient table.
#

# ## Summary
#
# - For most (patient, feature) models, **no term is significant** — expected
#   given how noisy single-object morphology features are.
# - Splitting significant models into **treatment-only**, **count-covariate-only**,
#   and **both** gives a per-(patient, feature) resolution of the density-confound
#   problem flagged at the aggregate level in `4.explore_leads.ipynb`.
# - The **clean treatment effect** table is the most defensible list of
#   per-patient drug responses in this whole analysis: significant and not
#   simultaneously explained by density.
# - The **per-patient significant-term profile** ranks patients by how much of
#   their explainable signal is genuine biology vs. density — use it to weight
#   how much to trust a given patient's hits elsewhere in this pipeline.
# - Full per-(patient, feature) significance table saved to
#   `4.linear_modeling/results/significant_variates/
#   {organoid,single_cell}_significant_variates_per_patient_feature.parquet`.
#

#

#
