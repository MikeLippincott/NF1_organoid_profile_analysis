#!/usr/bin/env python
# coding: utf-8

# ## Combine the profiles, viabilities, and platemap information
# this information is the new profiles will be already annotated in ad thus we will not need to do this step.

# In[1]:


import logging
import pathlib

import pandas as pd
from notebook_init_utils import bandicoot_check, init_notebook
from pycytominer import aggregate

root_dir, in_notebook = init_notebook()

if in_notebook:
    import tqdm.notebook as tqdm
else:
    import tqdm


# In[2]:


patient_ids = pd.read_csv(
    pathlib.Path(f"{root_dir}/data/patient_IDs.txt").resolve(strict=True),
    header=None,
    sep="\t",
    names=["patient_id"],
)["patient_id"].to_list()

viabilities_path = pathlib.Path(f"{root_dir}/data/viabilities/").resolve(strict=True)


# In[3]:


barcode_file_path = pathlib.Path(
    f"{root_dir}/config/platemaps/barcode_platemap.csv"
).resolve(strict=True)
platemap1_path = pathlib.Path(f"{root_dir}/config/platemaps/platemap1.csv").resolve(
    strict=True
)
platemap2_path = pathlib.Path(f"{root_dir}/config/platemaps/platemap2.csv").resolve(
    strict=True
)
barcode_df = pd.read_csv(barcode_file_path)
platemap1_df = pd.read_csv(platemap1_path)
platemap2_df = pd.read_csv(platemap2_path)
platemap_df_list = []

for patient in patient_ids:
    barcode = barcode_df.loc[
        barcode_df["patient_tumor_barcode"] == patient, "platemap_number"
    ].values[0]
    if barcode == "platemap1":
        tmp_df = platemap1_df.copy()
    elif barcode == "platemap2":
        tmp_df = platemap2_df.copy()
    else:
        raise ValueError(
            f"Unexpected barcode '{barcode}' for patient '{patient}' in barcode_platemap.csv"
        )
    tmp_df["patient_id"] = patient
    platemap_df_list.append(tmp_df)
platemap_df = pd.concat(platemap_df_list, axis=0)


# In[4]:


viabilities_file_path = pathlib.Path(
    f"{root_dir}/data/viabilities/raw_viabilities_combined.csv"
).resolve(strict=True)
viabilities_df = pd.read_csv(viabilities_file_path)

viabilities_df_list = []
for patient in patient_ids:
    # check if patient exists in viabilities_df
    if patient not in viabilities_df["Metadata_Biology_PatientTumor"].values:
        # skip this is expected, but log a warning
        logging.warning(f"Patient '{patient}' not found in viabilities_df, skipping.")
        continue
    # get the viabilities data for this patient
    tmp_df = viabilities_df.loc[
        viabilities_df["Metadata_Biology_PatientTumor"] == patient
    ].copy()
    # change DMSO dose to 1
    tmp_df.loc[tmp_df["Drug"] == "DMSO", "Concentration_uM"] = 1
    tmp_df.loc[tmp_df["Drug"] == "PD0325901", "Drug"] = "Mirdametinib"

    viabilities_df_list.append(tmp_df)

viabilities_df = pd.concat(viabilities_df_list, axis=0)
viabilities_df.rename(
    columns={"Viability_percentage": "Metadata_Viability_percentage"}, inplace=True
)


# In[5]:


# merge the viabilities with the platemap
platemap_viability_df = pd.merge(
    platemap_df,
    viabilities_df,
    how="left",
    left_on=["Treatment", "Dose", "patient_id"],
    right_on=["Drug", "Concentration_uM", "Metadata_Biology_PatientTumor"],
)

# Check (not just assume) that unmatched rows really are the expected
# empty/control wells, rather than a silent Treatment/Dose/patient_id
# naming mismatch that would disproportionately and invisibly drop real
# data. Done BEFORE dropping WellCol/WellPosition so the breakdown below
# can still see well position.
nan_rows = platemap_viability_df[platemap_viability_df.isna().any(axis=1)]
n_nan = len(nan_rows)
if n_nan > 0:
    breakdown = (
        nan_rows.groupby(["patient_id", "Treatment", "Dose"], dropna=False)
        .size()
        .sort_values(ascending=False)
    )
    logging.warning(
        f"{n_nan} of {len(platemap_viability_df)} row(s) "
        f"({n_nan / len(platemap_viability_df):.1%}) failed to merge platemap <-> "
        f"viability data and will be dropped. Breakdown by (patient_id, Treatment, Dose):\n"
        f"{breakdown.head(20)}"
    )
    if "WellPosition" in nan_rows.columns:
        non_b_well_nans = nan_rows[
            ~nan_rows["WellPosition"].astype(str).str.startswith("B")
        ]
        if len(non_b_well_nans) > 0:
            logging.warning(
                f"{len(non_b_well_nans)} of the {n_nan} dropped row(s) are NOT from a 'B' "
                f"well - this may indicate a real Treatment/Dose/patient_id naming "
                f"mismatch rather than the expected empty control wells. Inspect "
                f"`nan_rows` before trusting the drop below."
            )
        else:
            logging.info(
                "All dropped rows are from 'B' wells, matching the expected "
                "empty-control assumption."
            )

# now safe to drop these columns and the unmatched rows
platemap_viability_df = platemap_viability_df.drop(columns=["WellCol", "WellPosition"])
platemap_viability_df = platemap_viability_df.dropna().reset_index(drop=True)
platemap_viability_df.rename(
    columns={"Viability_percentage": "Metadata_Viability_percentage"}, inplace=True
)
# NOTE: this global, per-patient min-max column is kept here only for
# backward-compatibility with other notebooks that read
# combined_platemaps.parquet directly. The model-training pipeline further
# below does NOT use this column - it recomputes a fold-safe version from
# Metadata_Viability_percentage inside run_group_cv / run_random_split
# (see compute_fold_safe_target), since this global version is safe to use
# as-is for LOPO but leaks a held-out treatment's value into training
# targets for LOTO / random_split.
platemap_viability_df["min_max_viability"] = platemap_viability_df.groupby(
    "patient_id"
)["Metadata_Viability_percentage"].transform(
    lambda x: (x - x.min()) / (x.max() - x.min())
)
platemap_viability_df.drop(
    columns=["patient_id", "Drug", "Concentration_uM", "WellRow"], inplace=True
)
# save the combined platemaps
combined_platemaps_path = pathlib.Path(
    f"{root_dir}/data/viabilities/combined_platemaps.parquet"
).resolve()
# save this for use in later notebooks where viability is needed - no need to recombine the platemaps and viabilities each time
platemap_viability_df.to_parquet(combined_platemaps_path, index=False)


# ## Get all of the morphology profiles to work with

# In[6]:


# Consensus strata: one row per (patient, treatment) combination
consensus_strata_3D = [
    "Metadata_Biology_PatientTumor",
    "Metadata_Experiment_Treatment",
    "Metadata_Experiment_Dose",
    "Metadata_Experiment_Class",
    "Metadata_Experiment_Target",
    "Metadata_Experiment_TherapeuticCategories",
    "Metadata_Experiment_Unit",
]
consensus_strata_2D = [
    "Metadata_patient_tumor",
    "Metadata_treatment",
    "Metadata_dose",
    "Metadata_class",
    "Metadata_target",
    "Metadata_therapeutic_categories",
    "Metadata_dose_unit",
]


# In[7]:


consensus_profiles_3D_paths = pathlib.Path(
    f"{root_dir}/data/profiles_3D/all_patients/0.normalized_profiles/"
).resolve(strict=True)
consensus_profiles_3D_paths = list(consensus_profiles_3D_paths.glob("*.parquet"))
consensus_profiles_2D_paths = [
    pathlib.Path(
        f"{root_dir}/data/profiles_2D/all_patients/max_projection/organoid_profiles.parquet"
    ).resolve(strict=True),
    pathlib.Path(
        f"{root_dir}/data/profiles_2D/all_patients/max_projection/sc_profiles.parquet"
    ).resolve(strict=True),
]
paths_dict = {
    "3D": consensus_profiles_3D_paths,
    "2D": consensus_profiles_2D_paths,
}
paths_dict


# In[8]:


for dimension in paths_dict.keys():
    for profile_path in paths_dict[dimension]:
        print(f"Processing {profile_path.name} ({dimension})...")

        consensus_output_path = pathlib.Path(
            f"{root_dir}/3.viability_prediction_models/data/processed_profiles_{dimension}/{profile_path.stem.replace('_norm', '')}_consensus.parquet"
        ).resolve()
        consensus_output_path.parent.mkdir(parents=True, exist_ok=True)
        df = pd.read_parquet(profile_path)
        features_columns = [
            col for col in df.columns if not col.startswith("Metadata_")
        ]
        consensus_df = aggregate(
            population_df=df,
            strata=consensus_strata_3D if dimension == "3D" else consensus_strata_2D,
            features=features_columns,
            operation="median",
            output_file=consensus_output_path,
            output_type="parquet",
        )
        consensus_df = pd.read_parquet(consensus_output_path)
        # add viability to the consensus profiles
        if dimension == "3D":
            consensus_df = pd.merge(
                consensus_df,
                platemap_viability_df,
                how="left",
                left_on=[
                    "Metadata_Biology_PatientTumor",
                    "Metadata_Experiment_Treatment",
                    "Metadata_Experiment_Dose",
                    "Metadata_Experiment_Unit",
                ],
                right_on=["Metadata_Biology_PatientTumor", "Treatment", "Dose", "Unit"],
            ).drop(columns=["Treatment", "Dose", "Unit"])
        else:
            consensus_df = pd.merge(
                consensus_df,
                platemap_viability_df,
                how="left",
                left_on=[
                    "Metadata_patient_tumor",
                    "Metadata_treatment",
                    "Metadata_dose",
                    "Metadata_dose_unit",
                ],
                right_on=["Metadata_Biology_PatientTumor", "Treatment", "Dose", "Unit"],
            ).drop(
                columns=["Metadata_Biology_PatientTumor", "Treatment", "Dose", "Unit"]
            )
        consensus_df.to_parquet(consensus_output_path, index=False)
        print(consensus_df.shape)
