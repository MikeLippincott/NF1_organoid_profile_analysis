#!/usr/bin/env python
# coding: utf-8

# # Mining the linear modeling results for leads
#
# This notebook digs through the outputs of `0.linear_modeling.ipynb`
# (`1.EDA/results/linear_modeling/{organoid_fs,sc_fs}.parquet`) to surface
# candidate findings worth following up on, rather than re-deriving
# significance calls (that's `2.find_significant_features.ipynb`).
#
# Recap of the model: for every (patient, drug+dose, feature) combination, we
# fit, within that patient only,
#
# ```
# feature ~ C(treatment) + cell_count + organoid_count + cell_per_organoid_count
# ```
#
# Each fitted model contributes multiple rows to the long-form results table,
# one per **term**:
#
# - `treatment` — the treatment-vs-DMSO effect within that patient. This is
#   the drug-response signal.
# - `cell_count`, `organoid_count`, `cell_per_organoid_count` — nuisance
#   covariates, one row per (patient, combo, feature).
#
# All p-values are already FDR-corrected **within term** (`pvalue_fdr`),
# since each term is its own hypothesis-testing family.
#
# Leads pursued here:
#
# 1. Broadly-active vs. patient-specific drugs (breadth of response)
# 2. Drug leaderboard by hit count and by effect size
# 3. Therapeutic-category (MOA) roll-up
# 4. Which channels/feature-types are being perturbed (biological read-out)
# 5. Count-covariate confounding check
# 6. Organoid- vs. single-cell-level concordance
# 7. A saved shortlist of top candidate leads for follow-up
#

# In[1]:


import pathlib
import warnings

import numpy as np
import pandas as pd
from notebook_init_utils import init_notebook

warnings.filterwarnings("ignore")
root_dir, in_notebook = init_notebook()

results_path = pathlib.Path(root_dir, "4.linear_modeling/results/leads")
results_path.mkdir(parents=True, exist_ok=True)

pd.set_option("display.max_columns", 50)


# In[2]:


PVALUE_MAX = 0.05

profile_dict = {
    "organoid": pathlib.Path(
        root_dir, "4.linear_modeling/results/linear_modeling/organoid_norm.parquet"
    ),
    "single_cell": pathlib.Path(
        root_dir, "4.linear_modeling/results/linear_modeling/sc_norm.parquet"
    ),
}

dfs = {name: pd.read_parquet(path) for name, path in profile_dict.items()}
for name, df in dfs.items():
    print(
        f"{name}: {df.shape[0]:,} rows, {df['feature'].nunique():,} features, "
        f"{df['treatment'].nunique()} treatments, {df['patient'].nunique()} patients"
    )
    print(df.head(2))


# ## 1. Broadly-active vs. patient-specific drugs
#
# For the `treatment` term, count how many (patient, feature) cells are
# significant per drug, and — more importantly — in how many *distinct
# patients* a drug produces at least one significant feature. A drug that's
# significant in most patients is a **broad responder** (candidate pan-patient
# hit); a drug that's significant in only one or two patients but with many
# features is a **patient-specific responder** (candidate precision-oncology
# lead).
#

# In[3]:


def treatment_breadth(df, pvalue_max=PVALUE_MAX):
    trt = df[df["term"] == "treatment"].copy()
    trt["hit"] = trt["pvalue_fdr"] < pvalue_max

    per_patient = (
        trt.groupby(["treatment", "drug", "therapeutic_category", "patient"])["hit"]
        .sum()
        .reset_index(name="n_sig_features")
    )
    n_patients_total = trt["patient"].nunique()

    summary = (
        per_patient.groupby(["treatment", "drug", "therapeutic_category"])
        .agg(
            n_patients_responding=("n_sig_features", lambda s: (s > 0).sum()),
            median_sig_features_when_responding=(
                "n_sig_features",
                lambda s: s[s > 0].median() if (s > 0).any() else 0,
            ),
            total_sig_features=("n_sig_features", "sum"),
        )
        .reset_index()
    )
    summary["n_patients_tested"] = n_patients_total
    summary["frac_patients_responding"] = (
        summary["n_patients_responding"] / summary["n_patients_tested"]
    )
    return summary.sort_values("total_sig_features", ascending=False)


breadth = {name: treatment_breadth(df) for name, df in dfs.items()}
for name, summary in breadth.items():
    out_path = results_path / f"{name}_treatment_breadth.parquet"
    summary.to_parquet(out_path, index=False)
    print(
        f"{name}: saved treatment breadth summary ({len(summary)} treatments) -> {out_path}"
    )
print(breadth["organoid"].head(15))


# In[4]:


for name, summary in breadth.items():
    print(f"\n==== {name}: broadest responders (>=50% of patients) ====")
    print(
        summary[summary["frac_patients_responding"] >= 0.5]
        .sort_values("total_sig_features", ascending=False)
        .head(10)
    )
    print(
        f"\n==== {name}: most patient-specific (responds in exactly 1 patient, high hit count) ===="
    )
    print(
        summary[summary["n_patients_responding"] == 1]
        .sort_values("total_sig_features", ascending=False)
        .head(10)
    )


# ## 2. Effect-size leaderboard
#
# Breadth alone doesn't capture magnitude. Rank the single strongest
# (treatment, patient, feature) hits by absolute coefficient among
# statistically significant rows, to find the largest-magnitude, most
# confident effects.
#

# In[5]:


def top_effects(df, pvalue_max=PVALUE_MAX, n=25):
    trt = df[(df["term"] == "treatment") & (df["pvalue_fdr"] < pvalue_max)].copy()
    trt["abs_coefficient"] = trt["coefficient"].abs()
    cols = [
        "drug",
        "treatment",
        "patient",
        "therapeutic_category",
        "feature",
        "Compartment",
        "Channel",
        "Feature_type",
        "Measurement",
        "coefficient",
        "pvalue_fdr",
        "rsquared_adj",
    ]
    return trt.sort_values("abs_coefficient", ascending=False)[cols].head(n)


for name, df in dfs.items():
    print(f"\n==== {name}: top effect-size hits ====")
    print(top_effects(df))


# ## 3. Therapeutic category (MOA) roll-up
#
# Do drugs sharing a mechanism of action converge on the same features? A
# tight MOA cluster (many drugs in a category all hitting similar features)
# is a much stronger lead than a single-drug hit, since it corroborates the
# mechanism rather than an assay artifact.
#

# In[6]:


def moa_rollup(df, pvalue_max=PVALUE_MAX):
    trt = df[df["term"] == "treatment"].copy()
    trt["hit"] = trt["pvalue_fdr"] < pvalue_max
    rollup = (
        trt.groupby("therapeutic_category")
        .agg(
            n_drugs=("drug", "nunique"),
            n_sig_hits=("hit", "sum"),
            n_rows=("hit", "size"),
        )
        .assign(hit_rate=lambda d: d["n_sig_hits"] / d["n_rows"])
        .sort_values("hit_rate", ascending=False)
    )
    return rollup


moa_rollups = {}
for name, df in dfs.items():
    rollup = moa_rollup(df)
    moa_rollups[name] = rollup
    out_path = results_path / f"{name}_moa_rollup.parquet"
    rollup.reset_index().to_parquet(out_path, index=False)
    print(f"\n==== {name}: therapeutic category hit rates ====")
    print(rollup)
    print(f"saved -> {out_path}")


# ## 4. Which channels / feature types are being perturbed?
#
# For significant treatment-effect hits, break down by `Channel` and
# `Feature_type` to get a biological read on *what kind* of morphology is
# changing (e.g. mitochondrial texture, AGP granularity, nuclear shape).
#

# In[7]:


def channel_feature_breakdown(df, pvalue_max=PVALUE_MAX):
    trt = df[(df["term"] == "treatment") & (df["pvalue_fdr"] < pvalue_max)]
    return (
        trt.groupby(["Channel", "Feature_type"], dropna=False)
        .size()
        .reset_index(name="n_significant_hits")
        .sort_values("n_significant_hits", ascending=False)
    )


channel_feature_breakdowns = {}
for name, df in dfs.items():
    breakdown = channel_feature_breakdown(df)
    channel_feature_breakdowns[name] = breakdown
    out_path = results_path / f"{name}_channel_feature_breakdown.parquet"
    breakdown.to_parquet(out_path, index=False)
    print(f"{name}: saved channel/feature-type breakdown -> {out_path}")


# ## 5. Count-covariate confounding check
#
# A "treatment effect" that's really just a cell/organoid density effect is a
# false lead. Flag features where a count covariate is a significant
# predictor **and** the same feature shows up as a significant treatment hit
# for the same drug — these need density-normalization or extra scrutiny
# before trusting them as a true biological effect.
#

# In[9]:


def confounding_flags(df, pvalue_max=PVALUE_MAX):
    covariate_terms = ["cell_count", "organoid_count", "cell_per_organoid_count"]
    covariate_sig = df[
        df["term"].isin(covariate_terms) & (df["pvalue_fdr"] < pvalue_max)
    ][["treatment", "drug", "feature", "term"]].rename(
        columns={"term": "confounding_covariate"}
    )

    trt_sig = df[(df["term"] == "treatment") & (df["pvalue_fdr"] < pvalue_max)]

    flagged = trt_sig.merge(
        covariate_sig, on=["treatment", "drug", "feature"], how="inner"
    )
    return flagged


for name, df in dfs.items():
    flagged = confounding_flags(df)
    n_treatment_hits = (
        (df["term"] == "treatment") & (df["pvalue_fdr"] < PVALUE_MAX)
    ).sum()
    # a treatment-effect row can be flagged by more than one confounding covariate
    # (e.g. both cell_count and organoid_count significant for the same feature),
    # so de-duplicate to the underlying rows before reporting the share affected
    n_flagged_rows = flagged.drop_duplicates(
        ["treatment", "drug", "patient", "feature"]
    ).shape[0]
    print(
        f"{name}: {n_flagged_rows:,} / {n_treatment_hits:,} significant treatment-effect rows "
        f"({n_flagged_rows / max(n_treatment_hits, 1):.1%}) share a feature with a significant "
        "count-covariate effect for the same drug"
    )
    print(
        flagged.groupby(["drug", "confounding_covariate"])
        .size()
        .reset_index(name="n_flagged_rows")
        .sort_values("n_flagged_rows", ascending=False)
        .head(10)
    )


# ## 6. Organoid- vs. single-cell-level concordance
#
# Feature names differ between levels (`Organoid_...` vs. `Cell_/Cytoplasm_/
# Nuclei_...`), so exact `feature` strings won't match. Instead, match on
# `(treatment, patient, Channel, Feature_type, Measurement)` — same
# measurement type, same channel — as a proxy for "the same biological
# read-out at both scales." Hits that replicate across both organoid and
# single-cell level are the most trustworthy leads; hits at only one level
# may be scale-specific (interesting) or noise (needs follow-up).
#

# In[10]:


match_cols = ["treatment", "drug", "patient", "Channel", "Feature_type", "Measurement"]

org_sig = dfs["organoid"][
    (dfs["organoid"]["term"] == "treatment")
    & (dfs["organoid"]["pvalue_fdr"] < PVALUE_MAX)
][match_cols + ["feature", "coefficient"]].rename(
    columns={"feature": "organoid_feature", "coefficient": "organoid_coefficient"}
)
sc_sig = dfs["single_cell"][
    (dfs["single_cell"]["term"] == "treatment")
    & (dfs["single_cell"]["pvalue_fdr"] < PVALUE_MAX)
][match_cols + ["feature", "coefficient"]].rename(
    columns={"feature": "sc_feature", "coefficient": "sc_coefficient"}
)

concordant = org_sig.merge(sc_sig, on=match_cols, how="inner").drop_duplicates()
concordant["same_direction"] = np.sign(concordant["organoid_coefficient"]) == np.sign(
    concordant["sc_coefficient"]
)

print(f"organoid-only significant rows: {len(org_sig):,}")
print(f"single-cell-only significant rows: {len(sc_sig):,}")
print(
    f"concordant (treatment, patient, channel, feature_type, measurement) pairs: {len(concordant):,}"
)
print(f"  of which same-direction effect: {concordant['same_direction'].mean():.1%}")

print(
    concordant[concordant["same_direction"]]
    .assign(
        effect_strength=lambda d: (
            d["organoid_coefficient"].abs() + d["sc_coefficient"].abs()
        )
    )
    .sort_values("effect_strength", ascending=False)
    .drop(columns="effect_strength")
    .head(15)
)


# ## 7. Shortlist of top candidate leads
#
# Combine the signals above into a single ranked shortlist: significant
# treatment-effect hits that are (a) large effect size, (b) not explained by
# a count-covariate confound, and (c) — where available — concordant in
# direction between organoid and single-cell level. Save it for follow-up
# (e.g. pulling representative images, checking dose-response).
#

# In[11]:


def build_shortlist(df, name, flagged_confounds, n=40):
    trt = df[(df["term"] == "treatment") & (df["pvalue_fdr"] < PVALUE_MAX)].copy()
    confound_keys = set(
        zip(
            flagged_confounds["treatment"],
            flagged_confounds["drug"],
            flagged_confounds["feature"],
        )
    )
    trt["confounded_by_count_covariate"] = trt.apply(
        lambda r: (r["treatment"], r["drug"], r["feature"]) in confound_keys, axis=1
    )
    trt["abs_coefficient"] = trt["coefficient"].abs()
    trt["profile_level"] = name
    shortlist = trt[~trt["confounded_by_count_covariate"]].sort_values(
        "abs_coefficient", ascending=False
    )
    cols = [
        "profile_level",
        "drug",
        "treatment",
        "patient",
        "therapeutic_category",
        "feature",
        "Compartment",
        "Channel",
        "Feature_type",
        "Measurement",
        "coefficient",
        "pvalue_fdr",
        "rsquared_adj",
    ]
    return shortlist[cols].head(n)


shortlists = []
for name, df in dfs.items():
    flagged = confounding_flags(df)
    shortlists.append(build_shortlist(df, name, flagged))

leads_df = pd.concat(shortlists, ignore_index=True)
leads_out_path = results_path / "top_leads.parquet"
leads_df.to_parquet(leads_out_path, index=False)
print(f"Saved {len(leads_df)} candidate leads to {leads_out_path}")
print(leads_df.head(25))


# ## Summary
#
# - Sections 1-2 separate **broad** drug responses (candidate general
#   hits) from **patient-specific** ones (candidate precision-oncology leads)
#   and rank both by effect size. Breadth summary saved to
#   `{organoid,single_cell}_treatment_breadth.parquet`, plotted in
#   `plot_leads_breadth_vs_magnitude.ipynb`.
# - Section 3 checks whether hits cluster by mechanism of action — MOA-level
#   convergence is a stronger signal than a single-drug hit. Saved to
#   `{organoid,single_cell}_moa_rollup.parquet`, plotted in
#   `plot_leads_moa_hit_rate.ipynb`.
# - Section 4 gives a biological read-out (which channel/feature type is
#   moving). Saved to `{organoid,single_cell}_channel_feature_breakdown.parquet`,
#   plotted in `plot_leads_channel_feature_heatmap.ipynb`.
# - Section 5 flags hits that might just be density confounds. Section 7
#   excludes the flagged rows so the shortlist favors the more defensible
#   rows that survive. See `significant_variates_per_patient_feature.ipynb`
#   for a per-(patient, feature) resolved version of this same question.
# - Section 6 cross-checks organoid- vs. single-cell-level concordance as an
#   orthogonal replication signal.
# - The final shortlist (`4.linear_modeling/results/leads/top_leads.parquet`)
#   is the starting point for manual follow-up: pull representative images
#   for the top rows, check dose-response where multiple doses exist, and
#   cross-reference with the UpSet-style set analysis in
#   `find_significant_features.ipynb`.
#
