from typing import Dict, List

import numpy as np
import pandas as pd
from notebook_init_utils.notebook_init_utils import init_notebook
from scipy.stats import entropy

root_dir, in_notebook = init_notebook()
print(f"Root directory: {root_dir}")

_2d_dir = root_dir / "data" / "profiles_2D" / "all_patients" / "max_projection"
_3d_dir = root_dir / "data" / "profiles_3D" / "all_patients"
_3d_normal_dir = _3d_dir / "0.normalized_profiles"
_3d_fs_dir = _3d_dir / "1.feature_selected_profiles"
_3d_agg_dir = _3d_dir / "2.aggregated_profiles"
_entropy_dir = root_dir / "2.2d_vs_3d_analysis" / "results" / "entropy"

# (group_label, granularity, filename) for each 2D type x stage
_2d_entries = [
    ("2D", "organoid", "normal", _2d_dir / "organoid_profiles.parquet"),
    ("2D", "organoid", "fs", _2d_dir / "organoid_fs_profiles.parquet"),
    ("2D", "organoid", "agg", _2d_dir / "organoid_agg_profiles.parquet"),
    ("2D", "sc", "normal", _2d_dir / "sc_profiles.parquet"),
    ("2D", "sc", "fs", _2d_dir / "sc_fs_profiles.parquet"),
    ("2D", "sc", "agg", _2d_dir / "sc_agg_profiles.parquet"),
]

# (slug, group_label, granularity, file_prefix) for each 3D type
_3d_types = [
    ("organoid_handcrafted", "Handcrafted", "organoid", "organoid"),
    ("organoid_sammed", "DL (SAM-Med3D)", "organoid", "sammed_organoid"),
    ("sc_handcrafted", "Handcrafted", "sc", "sc"),
    ("sc_sammed", "DL (SAM-Med3D)", "sc", "sammed_sc"),
    ("sc_sammed_nucleocentric", "Nucleocentric DL", "sc", "sammed_nucleocentric"),
    (
        "sc_nucleocentric_morphem",
        "Nucleocentric MorphEM",
        "sc",
        "nucleocentric_morphem",
    ),
]

data_dict = {}

for group_label, granularity, stage, path in _2d_entries:
    key = f"2D_{granularity}" + ("" if stage == "normal" else f"_{stage}")
    suffix = "" if stage == "normal" else f"_{stage}"
    data_dict[key] = {
        "input": path,
        "output": (_entropy_dir / f"2D_{granularity}{suffix}_entropy.parquet"),
        "modality": "2D",
        "group_label": group_label,
        "granularity": granularity,
        "stage": stage,
    }

for slug, group_label, granularity, file_prefix in _3d_types:
    stage_paths = {
        "normal": _3d_normal_dir / f"{file_prefix}_norm_norm_profile.parquet",
        "fs": _3d_fs_dir / f"{file_prefix}_norm_fs_profiles.parquet",
        "agg": _3d_agg_dir / f"{file_prefix}_norm_sc_agg_profiles.parquet",
    }
    for stage, path in stage_paths.items():
        suffix = "" if stage == "normal" else f"_{stage}"
        key = f"3D_{slug}{suffix}"
        data_dict[key] = {
            "input": path,
            "output": (_entropy_dir / f"3D_{slug}{suffix}_entropy.parquet"),
            "modality": "3D",
            "group_label": group_label,
            "granularity": granularity,
            "stage": stage,
        }

print(f"{len(data_dict)} profile-type x stage combinations defined.")


def compute_feature_entropy(
    df: pd.DataFrame,
    n_bins: int = 50,
) -> pd.DataFrame:
    """Compute Shannon entropy for each feature column in the DataFrame.

    Feature columns are identified as those without the ``Metadata_`` prefix.
    Each feature's distribution is discretized into ``n_bins`` histogram bins,
    and Shannon entropy (base 2, in bits) is computed from the resulting
    probability distribution.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing both metadata and feature columns.
        Metadata columns are expected to carry the ``Metadata_`` prefix.
    n_bins : int, optional
        Number of histogram bins used for discretization. Default is 50.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ``feature_name`` (str) and ``entropy`` (float).
        Rows with no valid values are excluded.
    """
    # Separate feature columns from metadata columns
    feature_columns: List[str] = [
        col for col in df.columns if not col.startswith("Metadata_")
    ]

    results: Dict[str, List] = {"feature_name": [], "entropy": []}
    for feature in feature_columns:
        # Drop NaNs/infs and skip empty features
        values = df[feature].dropna().values
        values = values[np.isfinite(values)]
        if len(values) == 0:
            continue

        # Discretize into histogram bins; entropy() normalizes counts internally
        # Some upstream texture features contain near-float64-max outliers that
        # overflow numpy's bin-edge calculation even though they're finite -
        # skip those features rather than crash the whole run.
        try:
            counts, _ = np.histogram(values, bins=n_bins)
        except ValueError:
            continue
        feat_entropy: float = entropy(counts, base=2)

        results["feature_name"].append(feature)
        results["entropy"].append(feat_entropy)

    return pd.DataFrame(results)


_entropy_dir.mkdir(parents=True, exist_ok=True)

for key, config in data_dict.items():
    df = pd.read_parquet(config["input"])

    entropy_df = compute_feature_entropy(df, n_bins=50)
    entropy_df["imaging_modality"] = config["modality"]
    entropy_df["group_label"] = config["group_label"]
    entropy_df["granularity"] = config["granularity"]
    entropy_df["stage"] = config["stage"]

    entropy_df.to_parquet(config["output"], index=False)
    print(config["output"])


def derive_patient_column(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure a plain ``Metadata_patient`` column exists.

    2D and 3D profiles use different combined patient/tumor column names
    (``Metadata_patient_tumor`` vs ``Metadata_Biology_PatientTumor``). This
    derives ``Metadata_patient`` from whichever is present by dropping the
    tumor/replicate suffix (e.g. ``"NF0037_T1_CQ1"`` -> ``"NF0037"``),
    matching the same logic used in the sparse CCA refactor.
    """
    if "Metadata_patient" in df.columns:
        return df

    source_col = None
    for candidate in ("Metadata_patient_tumor", "Metadata_Biology_PatientTumor"):
        if candidate in df.columns:
            source_col = candidate
            break

    if source_col is None:
        raise ValueError(
            "No patient/tumor metadata column found to derive Metadata_patient from."
        )

    df["Metadata_patient"] = df[source_col].str.replace(r"_T\d+.*$", "", regex=True)
    return df


def compute_feature_entropy_per_patient(
    df: pd.DataFrame,
    patient_col: str = "Metadata_patient",
    n_bins: int = 50,
) -> pd.DataFrame:
    """Compute Shannon entropy per patient for each feature column.

    For each unique value in ``patient_col``, the feature values belonging
    to that patient are discretized into ``n_bins`` histogram bins and
    Shannon entropy (base 2) is computed from the resulting distribution.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with metadata and feature columns.
    patient_col : str, optional
        Name of the patient identifier column. Default is ``"Metadata_patient"``.
    n_bins : int, optional
        Number of histogram bins. Default is 50.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ``patient``, ``feature_name``, and ``entropy``.
    """
    feature_columns: List[str] = [
        col for col in df.columns if not col.startswith("Metadata_")
    ]

    records: List[dict] = []
    for patient, group in df.groupby(patient_col):
        for feature in feature_columns:
            values = group[feature].dropna().values
            values = values[np.isfinite(values)]
            if len(values) == 0:
                continue
            try:
                counts, _ = np.histogram(values, bins=n_bins)
            except ValueError:
                continue
            feat_entropy: float = entropy(counts, base=2)
            records.append(
                {"patient": patient, "feature_name": feature, "entropy": feat_entropy}
            )

    return pd.DataFrame(records)


for key, config in data_dict.items():
    df = pd.read_parquet(config["input"])
    df = derive_patient_column(df)

    per_patient_output = config["output"].parent / config["output"].name.replace(
        "_entropy.parquet", "_per_patient_entropy.parquet"
    )

    per_patient_df = compute_feature_entropy_per_patient(df, n_bins=50)
    per_patient_df["imaging_modality"] = config["modality"]
    per_patient_df["group_label"] = config["group_label"]
    per_patient_df["granularity"] = config["granularity"]
    per_patient_df["stage"] = config["stage"]

    per_patient_df.to_parquet(per_patient_output, index=False)
    print(per_patient_output)
