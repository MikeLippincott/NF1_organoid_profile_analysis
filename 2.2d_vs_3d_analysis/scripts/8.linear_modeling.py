#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import warnings

import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")  # Ignore all warnings
warnings.simplefilter("ignore")  # Additional suppression method

try:
    cfg = get_ipython().config
    in_notebook = True
except NameError:
    in_notebook = False
if in_notebook:
    from tqdm.notebook import tqdm
else:
    from tqdm import tqdm
# Get the current working directory
cwd = pathlib.Path.cwd()

if (cwd / ".git").is_dir():
    root_dir = cwd

else:
    root_dir = None
    for parent in cwd.parents:
        if (parent / ".git").is_dir():
            root_dir = parent
            break

# Check if a Git root directory was found
if root_dir is None:
    raise FileNotFoundError("No Git root directory found.")


# In[ ]:


profile_dict = {
    "organoid_fs": {
        "input_profile_path": pathlib.Path(
            root_dir,
            "data/profiles_3D/all_patients/1.feature_selected_profiles/organoid_norm_fs_profiles.parquet",
        ),
        "output_profile_path": pathlib.Path(
            root_dir, "2.2d_vs_3d_analysis/results/linear_modeling/organoid_fs.parquet"
        ),
    },
    "single_cell_fs": {
        "input_profile_path": pathlib.Path(
            root_dir,
            "data/profiles_3D/all_patients/1.feature_selected_profiles/sc_norm_fs_profiles.parquet",
        ),
        "output_profile_path": pathlib.Path(
            root_dir, "2.2d_vs_3d_analysis/results/linear_modeling/sc_fs.parquet"
        ),
    },
}


# ## Linear modeling
#
# We want to predict each feature given some information about the organoid per patient. We will use linear regression to do this.
# $y = X_0 * \beta_0 + X_1 * \beta_1 + ... + X_n * \beta_n + \epsilon$
#
# Where:
# $y$ = feature to predict
# $\beta_0$ = Intercept
# $X_1$ = The treatment (e.g. control, drug)
# $\beta_1$ = The coefficient for the treatment
# $\epsilon$ = The error term
#
# For each model(feature), we get the following statistics:
# - **R-squared**: Proportion of variance explained by the model.
# - **p-value**: Significance of the model.
# - **F-statistic**: Overall significance of the model.
# - **Coefficients**: Effect size of each predictor.

# In[ ]:


for profile in tqdm(profile_dict.keys(), desc="Loading profiles"):
    # set the output dictionary for linear modeling results
    # per profile
    linear_modeling_results_dict = {
        "patient": [],
        "treatment": [],
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
    # combine treatment, dose, and unit into a single column so that
    # different doses of the same treatment are modeled as distinct groups
    df["Metadata_treatment_full"] = (
        df["treatment"].astype(str)
        + "_"
        + df["Metadata_Experiment_Dose"].astype(str)
        + df["Metadata_Experiment_Unit"].astype(str)
    )
    metadata_columns = ["patient", "treatment"] + [
        col for col in df.columns if col.startswith("Metadata_")
    ]
    # rename feature columns as the "." dod not play nice with the formula
    for col in df.columns:
        new_col = col.replace(
            ".",
            "",  # we replace the "." with an empty string because it causes issues in the formula
            # the linear model interprets the "." as an operator and not as part of the column name
        )  # Replace . with empty string for compatibility in formula
        df.rename(columns={col: new_col}, inplace=True)

    for patient in tqdm(
        df["patient"].unique(), desc="Processing patients", unit="patient", leave=False
    ):
        df_patient = df.loc[df["patient"] == patient]

        # Filter for specific treatment/dose combinations
        dmso_label = df_patient.loc[
            df_patient["treatment"] == "DMSO", "Metadata_treatment_full"
        ].unique()[0]
        combo_list = [
            (dmso_label, i)
            for i in df_patient["Metadata_treatment_full"].unique()
            if i != dmso_label
        ]
        for combo in tqdm(
            combo_list,
            desc="Processing treatment combinations",
            unit="combo",
            leave=False,
        ):
            df_patient_trt = df_patient.loc[
                df_patient["Metadata_treatment_full"].isin(combo)
            ]
            # order the treatment column to ensure DMSO is first
            df_patient_trt["Metadata_treatment_full"] = pd.Categorical(
                df_patient_trt["Metadata_treatment_full"],
                categories=[dmso_label]
                + [
                    other_treatment
                    for other_treatment in df_patient[
                        "Metadata_treatment_full"
                    ].unique()
                    if other_treatment != dmso_label
                ],
            )
            for col in df_patient_trt.columns:
                if col not in metadata_columns:
                    # Prepare the formula for the linear model
                    formula = f"Q('{col}') ~ C(Metadata_treatment_full)"
                    # Import statsmodels and run the linear model
                    model = smf.ols(formula=formula, data=df_patient_trt)
                    results = model.fit()

                    linear_modeling_results_dict["patient"].append(patient)
                    linear_modeling_results_dict["treatment"].append(combo[1])
                    linear_modeling_results_dict["feature"].append(col)
                    linear_modeling_results_dict["rsquared"].append(results.rsquared)
                    linear_modeling_results_dict["rsquared_adj"].append(
                        results.rsquared_adj
                    )
                    linear_modeling_results_dict["fvalue"].append(results.fvalue)
                    linear_modeling_results_dict["pvalue"].append(
                        results.pvalues[f"C(Metadata_treatment_full)[T.{combo[1]}]"]
                    )
                    linear_modeling_results_dict["coefficient"].append(
                        results.params[
                            f"C(Metadata_treatment_full)[T.{combo[1]}]"
                        ].item()
                    )
                    linear_modeling_results_dict["intercept"].append(
                        results.params["Intercept"].item()
                    )
    linear_modeling_results_df = pd.DataFrame(linear_modeling_results_dict)
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

    # run FDR on the p-values
    pvals = linear_modeling_results_df["pvalue"].values
    _, pvals_fdr, _, _ = multipletests(pvals, method="fdr_bh")
    linear_modeling_results_df["pvalue_fdr"] = pvals_fdr
    # Save the updated DataFrame with FDR p-values
    profile_dict[profile]["output_profile_path"].parent.mkdir(
        parents=True, exist_ok=True
    )
    linear_modeling_results_df.to_parquet(
        profile_dict[profile]["output_profile_path"], index=False
    )
