#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import warnings

import pandas as pd
import statsmodels.formula.api as smf
from notebook_init_utils import init_notebook
from statsmodels.stats.multitest import multipletests

root_dir, in_notebook = init_notebook()
warnings.filterwarnings("ignore")  # Ignore all warnings
warnings.simplefilter("ignore")  # Additional suppression method

if in_notebook:
    from tqdm.notebook import tqdm
else:
    from tqdm import tqdm


# In[2]:


profile_dict = {
    "organoid_fs": {
        "input_profile_path": pathlib.Path(
            root_dir,
            "data/profiles_3D/all_patients/0.normalized_profiles/organoid_norm_fs_profiles.parquet",
        ),
        "output_profile_path": pathlib.Path(
            root_dir, "4.linear_modeling/results/linear_modeling/organoid_fs.parquet"
        ),
    },
    "single_cell_fs": {
        "input_profile_path": pathlib.Path(
            root_dir,
            "data/profiles_3D/all_patients/0.normalized_profiles/sc_norm_fs_profiles.parquet",
        ),
        "output_profile_path": pathlib.Path(
            root_dir, "4.linear_modeling/results/linear_modeling/sc_fs.parquet"
        ),
    },
}


# ## Linear Modeling
#
# We want to predict each feature using information about the organoid per patient. We use linear regression for this.
#
# **General form:**
#
# $$y = \beta_0 + X_1\beta_1 + X_2\beta_2 + \dots + X_n\beta_n + \epsilon$$
#
# **Model specification:**
#
# $$y \sim \text{txt} + \text{patient tumor} + \text{txt} \times \text{patient tumor} + \text{cell count} + \text{organoid count} + \text{cell/organoid count}$$
#
# $$y = \beta_0 + X_1\beta_1 + X_2\beta_2 + X_3\beta_3 + X_4\beta_4 + X_5\beta_5 + (X_1 X_2)\beta_6 + \epsilon$$
#
# **Where:**
#
# | $X$ | $\beta$ | Description |
# |-----|---------|--------------|
# | — | $\beta_0$ | Intercept |
# | $X_1$ | $\beta_1$ | Treatment (e.g., control, drug + dosage) |
# | $X_2$ | $\beta_2$ | Patient tumor |
# | $X_3$ | $\beta_3$ | Cell count |
# | $X_4$ | $\beta_4$ | Organoid count |
# | $X_5$ | $\beta_5$ | Cell/organoid count |
# | $X_1 X_2$ | $\beta_6$ | Interaction between treatment and patient tumor (product of $X_1$ and $X_2$) |
#
# $y$ = feature to predict, $\epsilon$ = error term
#
# **For each model (feature), we compute the following statistics:**
#
# - **R-squared**: Proportion of variance explained by the model.
# - **p-value**: Significance of the model.
# - **F-statistic**: Overall significance of the model.
# - **Coefficients**: Effect size of each predictor.

# In[3]:


lm_equation_terms = (
    "C(Metadata_treatment_full) + C(patient) "
    "+ cell_count + organoid_count + cell_per_organoid_count "
    "+ C(Metadata_treatment_full):C(patient)"
)


# In[4]:


tumor_type_dict = {
    "NF0014_T1": "cNF",
    "NF0014_T2": "pNF",
    "NF0016_T1": "pNF",
    "NF0018_T6": "cNF",
    "NF0021_T1": "cNF",
    "NF0030_T1": "Other",
    "NF0035_T1": "cNF",
    "NF0037_T1": "cNF",
    "NF0040_T1": "Other",
    "NF0055_T1": "pNF",
    "SARCO219_T2": "MPNST",
    "SARCO361_T1": "MPNST",
}


# In[ ]:


for profile in tqdm(profile_dict.keys(), desc="Loading profiles"):
    # set the output dictionary for linear modeling results
    # per profile. Results are stored long-form: one row per
    # (combo, feature, term), where "term" identifies which piece of
    # the model specification the coefficient/pvalue belongs to
    # (treatment, patient, cell_count, organoid_count, cell_per_organoid_count).
    if profile_dict[profile]["output_profile_path"].exists():
        continue  # skip if the output already exists

    linear_modeling_results_dict = {
        "term": [],
        "patient": [],
        "treatment": [],
        "drug": [],
        "therapeutic_category": [],
        "feature": [],
        "rsquared": [],
        "rsquared_adj": [],
        "fvalue": [],
        "pvalue": [],
        "coefficient": [],
        "intercept": [],
    }
    df = pd.read_parquet(profile_dict[profile]["input_profile_path"])
    df = df.rename(
        columns={
            "Metadata_Biology_PatientTumor": "patient",
            "Metadata_Experiment_Treatment": "treatment",
        }
    )
    # drop the NF0037_T1_CQ1 patient
    df = df.loc[df["patient"] != "NF0037_T1_CQ1"]
    # map each patient to its tumor type via the manually defined lookup
    df["tumor_type"] = df["patient"].map(tumor_type_dict)
    # combine treatment, dose, and unit into a single column so that
    # different doses of the same treatment are modeled as distinct groups
    df["Metadata_treatment_full"] = (
        df["treatment"].astype(str)
        + "_"
        + df["Metadata_Experiment_Dose"].astype(str)
        + df["Metadata_Experiment_Unit"].astype(str)
    )
    # map each combined treatment label back to its raw drug name and
    # therapeutic category (MOA) so the modeling results can be grouped
    # by drug/MOA downstream without re-merging metadata
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

    # build the count covariates used in the model specification:
    # cell count, organoid count, and cell/organoid count
    if profile == "single_cell_fs":
        # the single-cell profile has no organoid count of its own,
        # so pull it in from the organoid-level profile via patient + well
        organoid_counts_df = pd.read_parquet(
            profile_dict["organoid_fs"]["input_profile_path"],
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
        # drop rows whose (patient, well) has no matching organoid-level data
        df = df.dropna(subset=["Metadata_WellOrganoidCount"])
        df["cell_count"] = df["Metadata_Object_WellSingleCellCount"]
        df["organoid_count"] = df["Metadata_WellOrganoidCount"]
    else:
        df["cell_count"] = df["Metadata_Object_OrganoidSingleCellCount"]
        df["organoid_count"] = df["Metadata_WellOrganoidCount"]
    df["cell_per_organoid_count"] = df["cell_count"] / df["organoid_count"]

    count_columns = ["cell_count", "organoid_count", "cell_per_organoid_count"]
    metadata_columns = (
        ["patient", "treatment", "tumor_type"]
        + count_columns
        + [col for col in df.columns if col.startswith("Metadata_")]
    )
    # TODO: temporarily drop texture features
    df = df.drop(columns=[col for col in df.columns if "_Texture_" in col])
    # clip feature values to reduce the influence of extreme outliers on the model fit
    feature_columns = [col for col in df.columns if col not in metadata_columns]
    df[feature_columns] = df[feature_columns].clip(lower=-1e1, upper=1e1)
    # rename feature columns as the "." dod not play nice with the formula
    # the linear model interprets the "." as an operator and not as part of the column name
    # track the sanitized -> original name mapping so the original feature
    # names can be recovered after the results are loaded in any other
    # environment/notebook (the "feature" column in the output only ever
    # contains the sanitized names)
    sanitized_to_original_col_map = {}
    for col in df.columns:
        new_col = col.replace(
            ".", ""
        )  # Replace . with empty string for compatibility in formula
        sanitized_to_original_col_map[new_col] = col
        df.rename(columns={col: new_col}, inplace=True)
    # redefine the feature columns after renaming
    feature_columns = [col for col in df.columns if col not in metadata_columns]

    # Filter for specific treatment/dose combinations
    # DMSO's combined label is consistent across all patients
    dmso_label = df.loc[df["treatment"] == "DMSO", "Metadata_treatment_full"].unique()[
        0
    ]
    combo_list = [
        (dmso_label, i)
        for i in df["Metadata_treatment_full"].unique()
        if i != dmso_label
    ]
    for combo in tqdm(
        combo_list,
        desc="Processing treatment combinations",
        unit="combo",
        leave=False,
    ):
        drug_name = treatment_meta.loc[combo[1], "treatment"]
        therapeutic_category = treatment_meta.loc[
            combo[1], "Metadata_Experiment_TherapeuticCategories"
        ]

        # pool across all patients that have this treatment (plus DMSO) so that
        # patient and the treatment x patient interaction have variance to fit on
        df_trt = df.loc[df["Metadata_treatment_full"].isin(combo)]
        # order the treatment column to ensure DMSO is first (reference level)
        df_trt = df_trt.copy()
        df_trt["Metadata_treatment_full"] = pd.Categorical(
            df_trt["Metadata_treatment_full"], categories=[combo[0], combo[1]]
        )
        patients_in_combo = sorted(df_trt["patient"].unique())
        reference_patient = patients_in_combo[0]
        df_trt["patient"] = pd.Categorical(
            df_trt["patient"], categories=[reference_patient] + patients_in_combo[1:]
        )

        for col in tqdm(
            feature_columns, desc="Processing features", unit="feature", leave=False
        ):
            # Prepare the formula for the linear model:
            # y ~ txt + patient + tumor_type + txt:patient + cell_count + organoid_count + cell/organoid_count
            formula = f"Q('{col}') ~ {lm_equation_terms}"
            model = smf.ols(formula=formula, data=df_trt)
            results = model.fit()

            def add_row(term, patient, coefficient, pvalue):
                linear_modeling_results_dict["term"].append(term)
                linear_modeling_results_dict["patient"].append(patient)
                linear_modeling_results_dict["treatment"].append(combo[1])
                linear_modeling_results_dict["drug"].append(drug_name)
                linear_modeling_results_dict["therapeutic_category"].append(
                    therapeutic_category
                )
                linear_modeling_results_dict["feature"].append(col)
                linear_modeling_results_dict["rsquared"].append(results.rsquared)
                linear_modeling_results_dict["rsquared_adj"].append(
                    results.rsquared_adj
                )
                linear_modeling_results_dict["fvalue"].append(results.fvalue)
                linear_modeling_results_dict["pvalue"].append(pvalue)
                linear_modeling_results_dict["coefficient"].append(coefficient)
                linear_modeling_results_dict["intercept"].append(
                    results.params["Intercept"].item()
                )

            # term: treatment effect, resolved per patient (main effect for the
            # reference patient, main effect + interaction contrast for the rest)
            treatment_term = f"C(Metadata_treatment_full)[T.{combo[1]}]"
            for patient in patients_in_combo:
                if patient == reference_patient:
                    coefficient = results.params[treatment_term].item()
                    pvalue = results.pvalues[treatment_term].item()
                else:
                    interaction_term = f"{treatment_term}:C(patient)[T.{patient}]"
                    contrast = f"{treatment_term} + {interaction_term} = 0"
                    contrast_test = results.t_test(contrast)
                    coefficient = contrast_test.effect.item()
                    pvalue = contrast_test.pvalue.item()
                add_row("treatment", patient, coefficient, pvalue)

            # term: patient main effect (baseline shift vs. the reference
            # patient, independent of treatment)
            for patient in patients_in_combo:
                if patient == reference_patient:
                    continue
                patient_term = f"C(patient)[T.{patient}]"
                add_row(
                    "patient",
                    patient,
                    results.params[patient_term].item(),
                    results.pvalues[patient_term].item(),
                )

            # terms: the count covariates (not patient-specific)
            for covariate in count_columns:
                add_row(
                    covariate,
                    None,
                    results.params[covariate].item(),
                    results.pvalues[covariate].item(),
                )
    linear_modeling_results_df = pd.DataFrame(linear_modeling_results_dict)
    # map the sanitized feature names back to their original (pre-".": removal)
    # names so downstream consumers loading this parquet in any other
    # env/notebook can recover the native column name without needing access
    # to the sanitization logic above
    linear_modeling_results_df["feature_original"] = linear_modeling_results_df[
        "feature"
    ].map(sanitized_to_original_col_map)
    # split the feature column into multiple columns
    # feature names follow the pattern: Compartment_Channel_Feature_type_Measurement
    linear_modeling_results_df[
        ["Compartment", "Channel", "Feature_type", "Measurement"]
    ] = linear_modeling_results_df["feature"].str.split("_", n=3, expand=True)

    # if feature type is area shape then make the measurement the channel and
    # set the channel to None
    # this because area size shape features are not channel specific
    linear_modeling_results_df.loc[
        linear_modeling_results_df["Feature_type"] == "AreaSizeShape", "Measurement"
    ] = linear_modeling_results_df["Channel"]
    linear_modeling_results_df.loc[
        linear_modeling_results_df["Feature_type"] == "AreaSizeShape", "Channel"
    ] = None
    # set compartment to None if is adjacent
    # this is because adjacent features are not compartment specific
    linear_modeling_results_df.loc[
        linear_modeling_results_df["Compartment"] == "adjacent", "Compartment"
    ] = None

    # run FDR on the p-values, separately per term since each term is its
    # own hypothesis-testing family (different scale/behavior of p-values)
    linear_modeling_results_df["pvalue_fdr"] = float("nan")
    for term, group_index in linear_modeling_results_df.groupby("term").groups.items():
        pvals = linear_modeling_results_df.loc[group_index, "pvalue"].values
        _, pvals_fdr, _, _ = multipletests(pvals, method="fdr_bh")
        linear_modeling_results_df.loc[group_index, "pvalue_fdr"] = pvals_fdr
    # Save the updated DataFrame with FDR p-values
    profile_dict[profile]["output_profile_path"].parent.mkdir(
        parents=True, exist_ok=True
    )
    linear_modeling_results_df.to_parquet(
        profile_dict[profile]["output_profile_path"], index=False
    )
    # persist the sanitized -> original feature name mapping on its own so it
    # can be loaded independently of the results file in any other env/notebook
    feature_name_mapping_df = pd.DataFrame(
        {
            "feature": list(sanitized_to_original_col_map.keys()),
            "feature_original": list(sanitized_to_original_col_map.values()),
        }
    )
    feature_name_mapping_path = (
        profile_dict[profile]["output_profile_path"].parent
        / f"{profile}_feature_name_mapping.parquet"
    )
    feature_name_mapping_df.to_parquet(feature_name_mapping_path, index=False)
