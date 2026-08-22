#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
"""Area (2D, all 3 projections) vs. volume (3D), aggregated per patient x treatment.

Aggregating to patient x treatment (rather than joining on organoid object ID)
avoids the well-scoped-object-ID collision that affects a raw per-organoid join:
2D and 3D use separate pipelines/segmentations, so organoid IDs are not directly
comparable across modalities anyway. Patient x treatment is the only key both
modalities share.
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

results_dir = root_dir / "4.analysis" / "results" / "area_vs_volume"
results_dir.mkdir(parents=True, exist_ok=True)

AREA_COL = "Organoid_AreaShape_Area"
VOLUME_COL = "Organoid_NoChannel_AreaSizeShape_Volume"

# --- 2D area, per projection, per patient x treatment ---
area_rows = []
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
            / f"{prefix}_organoid.parquet"
        )
        if not f.exists():
            continue
        df = pd.read_parquet(f)
        if AREA_COL not in df.columns:
            continue
        df = harmonize_metadata(df, "2D", patient)
        g = (
            df.groupby(["Metadata_patient", "Metadata_treatment"])[AREA_COL]
            .agg(mean_area="mean", n_organoids_2d="count")
            .reset_index()
        )
        g["projection"] = projection
        area_rows.append(g)

area_by_patient_treatment = pd.concat(area_rows, ignore_index=True)

# --- 3D volume, per patient x treatment ---
volume_rows = []
patients_3d = list_patient_dirs(root_dir / "data" / "profiles_3D")
for patient in patients_3d:
    f = (
        root_dir
        / "data"
        / "profiles_3D"
        / patient
        / "5.normalized_profiles"
        / "organoid_norm.parquet"
    )
    if not f.exists():
        continue
    df = pd.read_parquet(f)
    if VOLUME_COL not in df.columns:
        continue
    df = harmonize_metadata(df, "3D", patient)
    g = (
        df.groupby(["Metadata_patient", "Metadata_treatment"])[VOLUME_COL]
        .agg(mean_volume="mean", n_organoids_3d="count")
        .reset_index()
    )
    volume_rows.append(g)

volume_by_patient_treatment = pd.concat(volume_rows, ignore_index=True)

# --- join (inner: only patient x treatment combos present in both modalities) ---
merged = area_by_patient_treatment.merge(
    volume_by_patient_treatment,
    on=["Metadata_patient", "Metadata_treatment"],
    how="inner",
)
merged.to_parquet(
    results_dir / "area_vs_volume_by_patient_treatment.parquet", index=False
)
print(
    f"Wrote {results_dir / 'area_vs_volume_by_patient_treatment.parquet'} ({len(merged)} rows)"
)

dropped_2d_patients = sorted(
    set(area_by_patient_treatment["Metadata_patient"]) - set(merged["Metadata_patient"])
)
dropped_3d_patients = sorted(
    set(volume_by_patient_treatment["Metadata_patient"])
    - set(merged["Metadata_patient"])
)
with open(results_dir / "area_vs_volume_dropped_patients.txt", "w") as fh:
    fh.write(
        f"Patients with 2D area but no matching 3D volume (any treatment overlap): {dropped_2d_patients}\n"
    )
    fh.write(
        f"Patients with 3D volume but no matching 2D area (any treatment overlap): {dropped_3d_patients}\n"
    )
print(f"2D-only patients: {dropped_2d_patients}")
print(f"3D-only patients: {dropped_3d_patients}")
