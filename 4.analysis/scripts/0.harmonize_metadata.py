#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
"""Section 0: metadata harmonization report.

Verifies that 2D and 3D metadata columns can be harmonized to a common schema
and writes results/metadata_column_map.csv documenting the renames applied
per source and per-patient coverage. The actual harmonize_metadata() /
parse_feature_2d() / parse_feature_3d() utilities live in utils_analysis.py
and are imported by every downstream script in this module.
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
    column_rename_report,
    harmonize_metadata,
    list_patient_dirs,
)

results_dir = root_dir / "4.analysis" / "results"
results_dir.mkdir(parents=True, exist_ok=True)

report_rows = [column_rename_report("2D"), column_rename_report("3D")]

coverage_rows = []

# 2D: verify harmonization on one organoid + one sc file per patient per projection
patients_2d = list_patient_dirs(root_dir / "data" / "profiles_2D")
for patient in patients_2d:
    for projection in PROJECTIONS:
        prefix = PROJECTION_FILE_PREFIX[projection]
        for kind in ["organoid", "sc"]:
            f = (
                root_dir
                / "data"
                / "profiles_2D"
                / patient
                / "5.normalized"
                / f"{prefix}_{kind}.parquet"
            )
            if not f.exists():
                coverage_rows.append(
                    {
                        "modality": "2D",
                        "patient": patient,
                        "projection": projection,
                        "kind": kind,
                        "status": "MISSING_FILE",
                    }
                )
                continue
            df = pd.read_parquet(f)
            harmonized = harmonize_metadata(df, "2D", patient)
            common_cols_present = [
                c
                for c in [
                    "Metadata_patient",
                    "Metadata_treatment",
                    "Metadata_timepoint",
                    "Metadata_dose",
                    "Metadata_class",
                ]
                if c in harmonized.columns
            ]
            coverage_rows.append(
                {
                    "modality": "2D",
                    "patient": patient,
                    "projection": projection,
                    "kind": kind,
                    "status": "OK",
                    "n_rows": len(df),
                    "common_cols_present": len(common_cols_present),
                }
            )

# 3D: verify harmonization on organoid_norm + sc_norm per patient
patients_3d = list_patient_dirs(root_dir / "data" / "profiles_3D")
for patient in patients_3d:
    for kind, fname in [
        ("organoid", "organoid_norm.parquet"),
        ("sc", "sc_norm.parquet"),
    ]:
        f = (
            root_dir
            / "data"
            / "profiles_3D"
            / patient
            / "5.normalized_profiles"
            / fname
        )
        if not f.exists():
            coverage_rows.append(
                {
                    "modality": "3D",
                    "patient": patient,
                    "projection": "n/a",
                    "kind": kind,
                    "status": "MISSING_FILE",
                }
            )
            continue
        df = pd.read_parquet(f)
        harmonized = harmonize_metadata(df, "3D", patient)
        common_cols_present = [
            c
            for c in [
                "Metadata_patient",
                "Metadata_treatment",
                "Metadata_timepoint",
                "Metadata_dose",
                "Metadata_class",
            ]
            if c in harmonized.columns
        ]
        coverage_rows.append(
            {
                "modality": "3D",
                "patient": patient,
                "projection": "n/a",
                "kind": kind,
                "status": "OK",
                "n_rows": len(df),
                "common_cols_present": len(common_cols_present),
            }
        )

rename_report = pd.concat(report_rows, ignore_index=True)
rename_report.to_csv(results_dir / "metadata_column_map.csv", index=False)

coverage_report = pd.DataFrame(coverage_rows)
coverage_report.to_csv(results_dir / "metadata_harmonization_coverage.csv", index=False)

print("Rename map:")
print(rename_report.to_string(index=False))
print(f"\nWrote {results_dir / 'metadata_column_map.csv'}")
n_missing = (coverage_report["status"] == "MISSING_FILE").sum()
print(
    f"Coverage check: {len(coverage_report)} (patient, projection, kind) combos checked, "
    f"{n_missing} missing files."
)
print(f"Wrote {results_dir / 'metadata_harmonization_coverage.csv'}")
