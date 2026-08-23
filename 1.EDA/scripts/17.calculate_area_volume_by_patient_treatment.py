#!/usr/bin/env python
# coding: utf-8

# Area (2D, all 3 projections) vs. volume (3D), at both organoid and single-cell
# level. Writes raw per-record area/volume tables (patient_tumor x treatment x
# dose x projection, no aggregation) - script 14 does the pairing (a full
# cross-join within each patient_tumor x treatment x dose group, since 2D and 3D
# use separate pipelines/segmentations and organoid/cell IDs are not directly
# comparable across modalities) and plots every point.

# In[1]:


import warnings

import pandas as pd

warnings.filterwarnings("ignore")


# In[2]:


from notebook_init_utils import init_notebook

root_dir, in_notebook = init_notebook()
from eda_helper_utils.utils_analysis import (
    PROJECTION_FILE_PREFIX,
    PROJECTIONS,
    harmonize_metadata,
    list_patient_dirs,
)

results_dir = root_dir / "1.EDA" / "results" / "area_vs_volume"
results_dir.mkdir(parents=True, exist_ok=True)

KEY_COLS = ["Metadata_patient_tumor", "Metadata_treatment", "Metadata_dose"]
COL_2D = {"organoid": "Organoid_AreaShape_Area", "cell": "Cells_AreaShape_Area"}
COL_3D = {
    "organoid": "Organoid_NoChannel_AreaSizeShape_Volume",
    "cell": "Cell_NoChannel_AreaSizeShape_Volume",
}

patients_2d = list_patient_dirs(root_dir / "data" / "profiles_2D")
patients_3d = list_patient_dirs(root_dir / "data" / "profiles_3D")

# --- 2D area, per projection, raw per-record (organoid and single-cell) ---
for kind, area_col in COL_2D.items():
    rows = []
    suffix = "organoid" if kind == "organoid" else "sc"
    for projection in PROJECTIONS:
        prefix = PROJECTION_FILE_PREFIX[projection]
        for patient in patients_2d:
            f = (
                root_dir
                / "data"
                / "profiles_2D"
                / patient
                / "5.normalized"
                / f"{prefix}_{suffix}.parquet"
            )
            if not f.exists():
                continue
            df = pd.read_parquet(f)
            if area_col not in df.columns:
                continue
            df = harmonize_metadata(df, "2D", patient)
            g = df[KEY_COLS + [area_col]].rename(columns={area_col: "area"})
            g["projection"] = projection
            rows.append(g)
    area_df = pd.concat(rows, ignore_index=True)
    area_df.to_parquet(results_dir / f"area_2D_{kind}.parquet", index=False)
    print(f"Wrote {results_dir / f'area_2D_{kind}.parquet'} ({len(area_df)} rows)")

# --- 3D volume, raw per-record (organoid and single-cell) ---
for kind, volume_col in COL_3D.items():
    rows = []
    fname = "organoid_norm.parquet" if kind == "organoid" else "sc_norm.parquet"
    for patient in patients_3d:
        f = (
            root_dir
            / "data"
            / "profiles_3D"
            / patient
            / "5.normalized_profiles"
            / fname
        )
        if not f.exists():
            continue
        df = pd.read_parquet(f)
        if volume_col not in df.columns:
            continue
        df = harmonize_metadata(df, "3D", patient)
        g = df[KEY_COLS + [volume_col]].rename(columns={volume_col: "volume"})
        rows.append(g)
    volume_df = pd.concat(rows, ignore_index=True)
    volume_df.to_parquet(results_dir / f"volume_3D_{kind}.parquet", index=False)
    print(f"Wrote {results_dir / f'volume_3D_{kind}.parquet'} ({len(volume_df)} rows)")

# --- coverage diagnostics (organoid-level patient_tumor keys only in one modality) ---
area_organoid = pd.read_parquet(results_dir / "area_2D_organoid.parquet")
volume_organoid = pd.read_parquet(results_dir / "volume_3D_organoid.parquet")
dropped_2d_only = sorted(
    set(area_organoid["Metadata_patient_tumor"])
    - set(volume_organoid["Metadata_patient_tumor"])
)
dropped_3d_only = sorted(
    set(volume_organoid["Metadata_patient_tumor"])
    - set(area_organoid["Metadata_patient_tumor"])
)
dropped_patients = pd.DataFrame(
    [{"patient": p, "present_in": "2D_only"} for p in dropped_2d_only]
    + [{"patient": p, "present_in": "3D_only"} for p in dropped_3d_only],
    columns=["patient", "present_in"],
)
dropped_patients.to_parquet(
    results_dir / "area_vs_volume_dropped_patients.parquet", index=False
)
print(f"2D-only patients: {dropped_2d_only}")
print(f"3D-only patients: {dropped_3d_only}")
