#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
"""Section 3 (cont'd): join derived cell counts (script 4) with viability
via combined_platemaps.parquet. Reports and logs dropped patients (those
in profiles but not in the platemap, and vice versa).
"""


# In[2]:

import warnings

import pandas as pd

warnings.filterwarnings("ignore")


# In[3]:


from notebook_init_utils import init_notebook

root_dir, in_notebook = init_notebook()

results_dir = root_dir / "4.analysis" / "results"
cell_counts = pd.read_parquet(results_dir / "cell_counts.parquet")
platemaps = pd.read_parquet(
    root_dir / "data" / "viabilities" / "combined_platemaps.parquet"
)

# Mean cell count per patient x treatment (across all organoids/wells)
mean_counts = (
    cell_counts.groupby(
        ["Metadata_patient_tumor", "Metadata_treatment", "modality", "projection"]
    )["cell_count"]
    .mean()
    .reset_index()
    .rename(columns={"cell_count": "mean_cell_count"})
)

joined = mean_counts.merge(
    platemaps,
    left_on=["Metadata_patient_tumor", "Metadata_treatment"],
    right_on=["patient_id", "Treatment"],
    how="inner",
)
joined.to_parquet(results_dir / "count_viability_joined.parquet", index=False)
print(f"Wrote {results_dir / 'count_viability_joined.parquet'} ({len(joined)} rows)")

profile_patients = set(mean_counts["Metadata_patient_tumor"].unique())
platemap_patients = set(platemaps["patient_id"].unique())
dropped_from_profiles = sorted(profile_patients - platemap_patients)
dropped_from_platemaps = sorted(platemap_patients - profile_patients)

with open(results_dir / "count_viability_dropped_patients.txt", "w") as fh:
    fh.write(
        "Patients in profiles but NOT in combined_platemaps.parquet (dropped from join):\n"
    )
    for p in dropped_from_profiles:
        fh.write(f"  {p}\n")
    fh.write(
        "\nPatients in combined_platemaps.parquet but NOT in profiles (unexpected, none expected):\n"
    )
    for p in dropped_from_platemaps:
        fh.write(f"  {p}\n")

print("Dropped from profiles (no viability/platemap coverage):", dropped_from_profiles)
print("In platemaps but not in profiles:", dropped_from_platemaps)
print(f"Wrote {results_dir / 'count_viability_dropped_patients.txt'}")
