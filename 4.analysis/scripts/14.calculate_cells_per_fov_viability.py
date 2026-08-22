#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
"""Viability vs. number of cells per field of view (FOV).

Cells per FOV = single-cell (nuclei) row count grouped by imaging field:
(Well, ImageNumber) for 2D, Well_FOV for 3D. This is a different granularity
than script 4's well-level average cell-per-organoid count - here every
imaged field is one observation, regardless of how many organoids it contains.
Mean cells-per-FOV per patient x treatment is then joined to viability via
combined_platemaps.parquet (see script 7 for the same join pattern).
"""

import sys
import warnings

import pandas as pd

# In[2]:


warnings.filterwarnings("ignore")


# In[3]:


from notebook_init_utils import init_notebook

root_dir, in_notebook = init_notebook()
sys.path.insert(0, str(root_dir / "4.analysis" / "scripts"))
from utils_analysis import (
    PROJECTION_FILE_PREFIX,
    PROJECTIONS,
    harmonize_metadata,
    list_patient_dirs,
)

results_dir = root_dir / "4.analysis" / "results" / "cells_per_fov"
results_dir.mkdir(parents=True, exist_ok=True)

# --- 2D: cells per (Well, ImageNumber) ---
rows_2d = []
patients_2d = list_patient_dirs(root_dir / "data" / "profiles_2D")
for projection in PROJECTIONS:
    prefix = PROJECTION_FILE_PREFIX[projection]
    for patient in patients_2d:
        f = (
            root_dir
            / "data"
            / "profiles_2D"
            / patient
            / "5.normalized"
            / f"{prefix}_sc.parquet"
        )
        if not f.exists():
            continue
        df = pd.read_parquet(
            f, columns=["Metadata_Well", "Metadata_ImageNumber", "Metadata_treatment"]
        )
        df = harmonize_metadata(df, "2D", patient)
        per_fov = (
            df.groupby(["Metadata_well", "Metadata_ImageNumber", "Metadata_treatment"])
            .size()
            .reset_index(name="cells_per_fov")
        )
        g = (
            per_fov.groupby("Metadata_treatment")["cells_per_fov"]
            .mean()
            .reset_index()
            .rename(columns={"cells_per_fov": "mean_cells_per_fov"})
        )
        g["Metadata_patient_tumor"] = patient
        g["modality"] = "2D"
        g["projection"] = projection
        rows_2d.append(g)

cells_per_fov_2d = pd.concat(rows_2d, ignore_index=True)

# --- 3D: cells per Well_FOV ---
rows_3d = []
patients_3d = list_patient_dirs(root_dir / "data" / "profiles_3D")
for patient in patients_3d:
    f = (
        root_dir
        / "data"
        / "profiles_3D"
        / patient
        / "5.normalized_profiles"
        / "sc_norm.parquet"
    )
    if not f.exists():
        continue
    df = pd.read_parquet(
        f, columns=["Metadata_Experiment_WellFOV", "Metadata_Experiment_Treatment"]
    )
    df = harmonize_metadata(df, "3D", patient)
    per_fov = (
        df.groupby(["Metadata_Experiment_WellFOV", "Metadata_treatment"])
        .size()
        .reset_index(name="cells_per_fov")
    )
    g = (
        per_fov.groupby("Metadata_treatment")["cells_per_fov"]
        .mean()
        .reset_index()
        .rename(columns={"cells_per_fov": "mean_cells_per_fov"})
    )
    g["Metadata_patient_tumor"] = patient
    g["modality"] = "3D"
    g["projection"] = None
    rows_3d.append(g)

cells_per_fov_3d = pd.concat(rows_3d, ignore_index=True)

cells_per_fov = pd.concat([cells_per_fov_2d, cells_per_fov_3d], ignore_index=True)
cells_per_fov.to_parquet(results_dir / "cells_per_fov.parquet", index=False)
print(f"Wrote {results_dir / 'cells_per_fov.parquet'} ({len(cells_per_fov)} rows)")

# --- join with viability ---
platemaps = pd.read_parquet(
    root_dir / "data" / "viabilities" / "combined_platemaps.parquet"
)
joined = cells_per_fov.merge(
    platemaps,
    left_on=["Metadata_patient_tumor", "Metadata_treatment"],
    right_on=["patient_id", "Treatment"],
    how="inner",
)
joined.to_parquet(results_dir / "cells_per_fov_viability_joined.parquet", index=False)
print(
    f"Wrote {results_dir / 'cells_per_fov_viability_joined.parquet'} ({len(joined)} rows)"
)

profile_patients = set(cells_per_fov["Metadata_patient_tumor"].unique())
platemap_patients = set(platemaps["patient_id"].unique())
dropped_from_profiles = sorted(profile_patients - platemap_patients)
with open(results_dir / "cells_per_fov_viability_dropped_patients.txt", "w") as fh:
    fh.write(
        "Patients with cells-per-FOV data but NOT in combined_platemaps.parquet (dropped from join):\n"
    )
    for p in dropped_from_profiles:
        fh.write(f"  {p}\n")
print("Dropped (no viability/platemap coverage):", dropped_from_profiles)
