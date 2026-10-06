"""One render function per analysis type, grouped into one function per module."""

import re

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from background import subpanel_background
from data_io import (
    DIFFERENTIAL_PRODUCED_BY,
    DIFFERENTIAL_RESULTS,
    EDA_RESULTS,
    EXCLUDED_PATIENTS,
    Dataset,
    load_barcode_platemap,
    load_dataset,
    load_platemaps,
    parquet_columns,
    registry,
)
from palettes import (
    BIOLOGICAL_TERMS,
    TECHNICAL_TERMS,
    TERM_ORDER,
    TREATMENT_CLASS_DEFAULT,
    TREATMENT_CLASS_MAP,
    TREATMENT_MOA_MAP,
    TUMOR_TYPE_LOOKUP,
    TUMOR_TYPE_PALETTE,
    humanize_label,
    palette_for,
)
from plotly.subplots import make_subplots
from plots import (
    COL_TRACK_NAMES,
    NONE,
    ROW_TRACK_NAMES,
    annotated_significance_heatmap,
    apply_global_filters,
    correlation_heatmap,
    explorer,
    local_subset,
    missing_notice,
    platemap_heatmap,
    png_download,
    regroup_combos_by_identity,
    upset_plot,
)

Filters = dict[str, list[str]]


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------
def dataset_section(
    datasets: dict[str, Dataset],
    key: str,
    filters: Filters,
    produced_by: str,
    directory,
    defaults_for=None,
    kinds=None,
    columns: tuple[str, ...] | None = None,
    prepare=None,
) -> None:
    """Pick one of several result files and explore it."""
    if not datasets:
        missing_notice(key.replace("_", " "), produced_by, directory)
        return
    label = st.selectbox(
        "Dataset", list(datasets), key=f"{key}_dataset", format_func=humanize_label
    )
    ds = datasets[label]
    df = load_dataset(str(ds.path), columns)
    if prepare:
        df = prepare(df, label)
    defaults = defaults_for(label, df) if defaults_for else {}
    explorer(df, f"{key}_{label}", filters, kinds=kinds, defaults=defaults, title=label)


def _first(df: pd.DataFrame, *names: str) -> str | None:
    return next((n for n in names if n in df.columns), None)


# ---------------------------------------------------------------------------
# 0.Overview -- static experiment design (not computed, just config/platemaps/*)
# ---------------------------------------------------------------------------
def platemap_section(filters: Filters) -> None:
    plates = load_platemaps()
    if not plates:
        missing_notice("platemaps", "(none)", "config/platemaps")
        return
    barcodes = load_barcode_platemap()
    name = st.selectbox(
        "Platemap", list(plates), key="overview_platemap", format_func=humanize_label
    )
    plate = plates[name]
    patients = sorted(
        barcodes.loc[barcodes["platemap_number"] == name, "patient_tumor"]
    )
    if patients:
        st.caption(f"Run on {len(patients)} patient(s): {', '.join(patients)}")
    palette = palette_for("treatment", sorted(plate["Treatment"].unique())) or {}
    fig = platemap_heatmap(plate, palette, title=humanize_label(name))
    st.plotly_chart(fig, width="stretch", key="overview_platemap_chart")
    png_download(fig, "overview_platemap", f"{name}_layout")


def drugs_section(filters: Filters) -> None:
    plates = load_platemaps()
    if not plates:
        missing_notice("platemaps", "(none)", "config/platemaps")
        return
    combined = pd.concat(plates.values(), ignore_index=True)
    drugs = (
        combined[combined["Treatment"] != "DMSO"]
        .drop_duplicates(["Treatment", "Dose", "Unit"])
        .assign(
            moa=lambda d: d["Treatment"].map(TREATMENT_MOA_MAP),
            drug_class=lambda d: (
                d["Treatment"].map(TREATMENT_CLASS_MAP).fillna(TREATMENT_CLASS_DEFAULT)
            ),
        )
        .sort_values(["moa", "Treatment", "Dose"])
        .rename(
            columns={
                "Treatment": "treatment",
                "Dose": "dose",
                "Unit": "unit",
                "drug_class": "class",
            }
        )
    )
    drugs = local_subset(
        drugs[["treatment", "class", "moa", "dose", "unit"]], "overview_drugs"
    )
    st.caption(f"{drugs['treatment'].nunique()} drugs across {len(plates)} platemap(s)")
    st.dataframe(drugs, width="stretch", hide_index=True)


def patients_section(filters: Filters) -> None:
    barcodes = load_barcode_platemap()
    if barcodes.empty:
        missing_notice("patients", "(none)", "config/platemaps/barcode_platemap.csv")
        return
    st.caption(f"{len(barcodes)} patient tumor samples")
    st.dataframe(
        barcodes.rename(
            columns={"platemap_number": "platemap", "tumor_type": "tumor manifestation"}
        ).sort_values("patient_tumor"),
        width="stretch",
        hide_index=True,
    )


OVERVIEW_SECTIONS = {
    "Platemap": platemap_section,
    "Drugs": drugs_section,
    "Patients & tumor manifestations": patients_section,
}


# ---------------------------------------------------------------------------
# 1.EDA
# ---------------------------------------------------------------------------
_DIMENSION_RE = re.compile(r"(?:^|_)(2D|3D)(?:_|$)")


def _dimension_of(key: str) -> str | None:
    """ "2D"/"3D" token in a dataset key (e.g. ``patient_specific_2D_maxproj_scfs_umap``
    -> "2D"), or None if the key doesn't carry one."""
    m = _DIMENSION_RE.search(key)
    return m.group(1) if m else None


def _default_index(options: list, preferred) -> int:
    """Position of ``preferred`` in ``options``, or 0 when it isn't there."""
    return options.index(preferred) if preferred in options else 0


# default profile for the UMAP / PCA / correlation tabs: pooled 3D single-cell
# feature-selected profiles ("combined", i.e. not per patient)
DEFAULT_UMAP = "3D_scfs_umap"
DEFAULT_PCA = "3D_1.feature_selected_profiles_sc_norm_fs_profiles_embeddings"
DEFAULT_PAIRS = "3D_sc_correlation_pairs"


def umap_section(filters: Filters) -> None:
    datasets = registry()["umap"]
    dims = sorted({d for d in (_dimension_of(k) for k in datasets) if d})
    if dims:
        dim = st.radio(
            "Dimensionality",
            dims,
            index=_default_index(dims, "3D"),
            horizontal=True,
            key="umap_dim",
        )
        datasets = {k: v for k, v in datasets.items() if _dimension_of(k) == dim}
    if not datasets:
        missing_notice("umap", "1.EDA/scripts/0.generate_umap.py", EDA_RESULTS / "umap")
        return

    label = st.selectbox(
        "Dataset",
        list(datasets),
        index=_default_index(list(datasets), DEFAULT_UMAP),
        key="umap_dataset",
        format_func=humanize_label,
    )
    df = load_dataset(str(datasets[label].path))

    is_per_patient = label.startswith("patient_specific")
    individual = True
    if is_per_patient:
        individual = st.checkbox(
            "Individual patients (one plot per patient, each with its own axes "
            "-- combining them is misleading since each patient has its own "
            "independent UMAP fit)",
            value=True,
            key="umap_pp_individual",
        )

    facet_by_patient = is_per_patient and individual
    defaults = {
        "kind": "scatter",
        "x": "UMAP1",
        "y": "UMAP2",
        "color": (
            _first(df, "treatment")
            if facet_by_patient
            else _first(df, "treatment", "patient_tumor")
        ),
        "facet": "patient_tumor" if facet_by_patient else None,
    }
    explorer(
        df,
        f"umap_{label}_{'ind' if facet_by_patient else 'comb'}",
        filters,
        kinds=["scatter", "histogram", "heatmap"],
        defaults=defaults,
        title=label,
        free_facet_axes=facet_by_patient,
    )


def pca_section(filters: Filters) -> None:
    datasets = registry()["pca"]
    dims = sorted({d for d in (_dimension_of(k) for k in datasets) if d})
    if dims:
        dim = st.radio(
            "Dimensionality",
            dims,
            index=_default_index(dims, "3D"),
            horizontal=True,
            key="pca_dim",
        )
        datasets = {k: v for k, v in datasets.items() if _dimension_of(k) == dim}
    if not datasets:
        missing_notice("PCA", "1.EDA/scripts/2.generate_pca.py", EDA_RESULTS / "pca")
        return
    label = st.selectbox(
        "Dataset",
        list(datasets),
        index=_default_index(list(datasets), DEFAULT_PCA),
        key="pca_dataset",
        format_func=humanize_label,
    )
    path = datasets[label].path
    df = load_dataset(str(path))
    variance_path = path.with_name(
        path.name.replace("_embeddings", "_explained_variance")
    )
    variance = None
    if variance_path.exists():
        variance = pd.read_parquet(variance_path).iloc[0]
    x_col, y_col = st.columns(2)
    pcs = [c for c in df.columns if c.startswith("PC") and c[2:].isdigit()]
    x_pc = x_col.selectbox("X component", pcs, index=0, key="pca_x_pc")
    y_pc = y_col.selectbox(
        "Y component", pcs, index=min(1, len(pcs) - 1), key="pca_y_pc"
    )
    explorer(
        df,
        f"pca_{label}",
        filters,
        kinds=["scatter", "histogram", "heatmap"],
        defaults={
            "kind": "scatter",
            "x": x_pc,
            "y": y_pc,
            "color": _first(df, "treatment", "patient_tumor"),
        },
        title=label,
    )
    if variance is not None:
        explained = variance.rename(lambda n: n.replace("_explained_variance", ""))
        st.caption(
            f"Explained variance: {x_pc} {explained.get(x_pc, float('nan')):.1%}, "
            f"{y_pc} {explained.get(y_pc, float('nan')):.1%}"
        )
        scree = px.bar(
            x=explained.index,
            y=explained.values,
            labels={"x": "Component", "y": "Explained variance ratio"},
            title="Scree plot",
        )
        subpanel_background("pca_scree")
        st.plotly_chart(scree, width="stretch", key=f"pca_scree_{label}")


def correlation_section(filters: Filters) -> None:
    corr_dir = EDA_RESULTS / "correlation"
    pair_files = (
        sorted(corr_dir.glob("*_correlation_pairs.parquet"))
        if corr_dir.exists()
        else []
    )
    matrix_files = (
        sorted(corr_dir.glob("*per_patient_correlation_matrices.parquet"))
        if corr_dir.exists()
        else []
    )
    if not pair_files and not matrix_files:
        missing_notice(
            "correlation matrices",
            "1.EDA/scripts/4.calculate_correlation_matrix.py",
            corr_dir,
        )
        return

    mode = st.radio(
        "Matrix type",
        ["Replicate/treatment-level (pairs)", "Per-patient single-cell"],
        horizontal=True,
        key="corr_mode",
    )
    subpanel_background(
        "corr_pairs" if mode.startswith("Replicate") else "corr_per_patient"
    )
    if mode.startswith("Replicate"):
        _pairs_heatmap(pair_files, filters)
    else:
        _per_patient_heatmap(matrix_files)


def _pairs_heatmap(pair_files, filters: Filters, preset: dict | None = None) -> None:
    """Sample-by-sample heatmap. ``preset`` maps group columns to the value
    selected by default (e.g. ``{"profile_type": "consensus"}``)."""
    if not pair_files:
        st.info("No `*_correlation_pairs.parquet` files found.")
        return
    pairs_path = st.selectbox(
        "Pairs file",
        pair_files,
        index=_default_index([p.stem for p in pair_files], DEFAULT_PAIRS),
        format_func=lambda p: humanize_label(p.stem),
        key="corr_pairs_file",
    )
    samples_path = pairs_path.with_name(pairs_path.name.replace("_pairs", "_samples"))
    if not samples_path.exists():
        st.info(f"Missing sample metadata file `{samples_path.name}`.")
        return
    pairs = pd.read_parquet(pairs_path)
    samples = load_dataset(str(samples_path))
    group_cols = [
        c for c in pairs.columns if c not in ("sample_i", "sample_j", "correlation")
    ]
    selection = {}
    cols = st.columns(len(group_cols))
    for col, name in zip(cols, group_cols):
        options = sorted(pairs[name].unique())
        selection[name] = col.selectbox(
            name,
            options,
            index=_default_index(options, (preset or {}).get(name)),
            key=f"corr_{name}",
        )
    pair_mask = np.ones(len(pairs), dtype=bool)
    sample_mask = np.ones(len(samples), dtype=bool)
    for name, value in selection.items():
        pair_mask &= (pairs[name] == value).to_numpy()
        sample_mask &= (samples[name] == value).to_numpy()
    pairs = pairs[pair_mask]
    samples = samples[sample_mask].reset_index(drop=True)

    n = int(samples["sample_index"].max()) + 1
    matrix = np.full((n, n), np.nan)
    matrix[pairs["sample_i"].to_numpy(), pairs["sample_j"].to_numpy()] = pairs[
        "correlation"
    ].to_numpy()
    matrix[pairs["sample_j"].to_numpy(), pairs["sample_i"].to_numpy()] = pairs[
        "correlation"
    ].to_numpy()

    filtered = local_subset(apply_global_filters(samples, filters), "corr")
    order_options = [
        c
        for c in filtered.columns
        if 1 < filtered[c].nunique() <= 100 and c not in group_cols
    ]
    order_by = st.selectbox(
        "Order/label samples by",
        order_options,
        index=0 if order_options else None,
        key="corr_order",
    )
    if order_by:
        filtered = filtered.sort_values(order_by, kind="stable")
    idx = filtered["sample_index"].to_numpy()
    if len(idx) == 0:
        st.warning("No samples left after subsetting.")
        return
    sub = matrix[np.ix_(idx, idx)]
    labels = (
        filtered[order_by].astype(str).to_numpy()
        if order_by
        else np.arange(len(idx)).astype(str)
    )
    labels = [f"{i}: {v}" for i, v in enumerate(labels)]

    c1, c2 = st.columns(2)
    cluster = c1.checkbox(
        "Hierarchical clustering (rows & cols)", key="corr_pairs_cluster"
    )
    track_options = [c for c in order_options if filtered[c].nunique() <= 60]
    tracks = c2.multiselect(
        "Color bars (rows & cols)", track_options, key="corr_pairs_tracks"
    )
    meta = (
        filtered[tracks].reset_index(drop=True) if tracks else pd.DataFrame(index=idx)
    )
    fig = correlation_heatmap(
        sub,
        labels,
        meta,
        [(t, t, t) for t in tracks],
        cluster,
        title=" | ".join(f"{k}={v}" for k, v in selection.items()),
    )
    st.plotly_chart(fig, width="content", key="corr_pairs_chart")
    png_download(fig, "corr_pairs", "correlation_heatmap")


def _per_patient_index(path: str) -> pd.DataFrame:
    """Small (variant, patient, n_samples) table; skips the big matrix columns."""
    df = pd.read_parquet(path, columns=["variant", "patient", "n_samples"])
    return df[~df["patient"].astype(str).isin(EXCLUDED_PATIENTS)].reset_index(drop=True)


def _per_patient_row(path: str, variant: str, patient: str) -> dict:
    """Read only the one matrix the user picked, not the whole file."""
    row = pd.read_parquet(
        path,
        columns=["n_samples", "correlation", "treatment"],
        filters=[("variant", "==", variant), ("patient", "==", patient)],
    ).iloc[0]
    return {
        "n": int(row["n_samples"]),
        "correlation": np.asarray(row["correlation"]),
        "treatment": np.asarray(row["treatment"]),
    }


def _per_patient_heatmap(matrix_files) -> None:
    if not matrix_files:
        st.info("No per-patient correlation matrix file found.")
        return
    path = str(matrix_files[0])
    index = _per_patient_index(path)
    c1, c2 = st.columns(2)
    variant = c1.selectbox(
        "Normalization variant",
        sorted(index["variant"].unique()),
        key="corr_variant",
        format_func=humanize_label,
    )
    patient = c2.selectbox(
        "Patient",
        sorted(index.loc[index["variant"] == variant, "patient"].unique()),
        key="corr_patient",
    )
    row = _per_patient_row(path, variant, patient)
    n = row["n"]
    matrix = row["correlation"].reshape(n, n)
    treatments = row["treatment"]
    order_treatments = st.multiselect(
        "Keep treatments", sorted(set(treatments)), key="corr_pp_treat"
    )
    keep = (
        np.isin(treatments, order_treatments)
        if order_treatments
        else np.ones(n, dtype=bool)
    )
    idx = np.where(keep)[0]
    idx = idx[np.argsort(treatments[idx], kind="stable")]

    c3, c4 = st.columns(2)
    cluster = c3.checkbox(
        "Hierarchical clustering (rows & cols)", key="corr_pp_cluster"
    )
    show_track = c4.checkbox("Color bar: treatment", value=True, key="corr_pp_track")
    meta = pd.DataFrame({"treatment": treatments[idx]})
    tracks = [("Treatment", "treatment", "treatment")] if show_track else []
    order_note = "clustered" if cluster else "sorted by treatment"
    fig = correlation_heatmap(
        matrix[np.ix_(idx, idx)],
        [str(i) for i in idx],
        meta,
        tracks,
        cluster,
        title=f"{patient} - {humanize_label(variant)} ({len(idx)} cells, {order_note})",
    )
    st.plotly_chart(fig, width="content", key="corr_pp_chart")
    png_download(fig, "corr_pp", "correlation_heatmap_per_patient")


def cell_counts_section(filters: Filters) -> None:
    def defaults(label, df):
        y = _first(df, "n_cells", "total_cells", "mean_cells_per_organoid")
        return {
            "kind": "box",
            "x": _first(df, "treatment"),
            "y": y,
            "color": _first(df, "patient_tumor"),
        }

    dataset_section(
        registry()["cell_counts"],
        "cell counts",
        filters,
        "1.EDA/scripts/7.generate_cell_counts.py",
        EDA_RESULTS / "cell_counts",
        defaults,
    )


def area_volume_section(filters: Filters) -> None:
    def defaults(label, df):
        return {
            "kind": "violin",
            "x": "treatment",
            "y": _first(df, "volume", "area"),
            "color": _first(df, "treatment"),
            "facet": _first(df, "patient_tumor"),
        }

    dataset_section(
        registry()["area_vs_volume"],
        "area and volume",
        filters,
        "1.EDA/scripts/17.calculate_area_volume_by_patient_treatment.py",
        EDA_RESULTS / "area_vs_volume",
        defaults,
    )


def neighbors_section(filters: Filters) -> None:
    def defaults(label, df):
        measure = next(
            (
                c
                for c in df.columns
                if c.startswith(("Organoid_", "Nuclei_"))
                and pd.api.types.is_numeric_dtype(df[c])
            ),
            None,
        )
        return {
            "kind": "box",
            "x": "treatment",
            "y": measure,
            "color": _first(df, "patient_tumor"),
        }

    dataset_section(
        registry()["neighbors"],
        "neighbors",
        filters,
        "1.EDA/scripts/13.calculate_neighbor_features.py",
        EDA_RESULTS / "neighbors",
        defaults,
    )


def intensity_section(filters: Filters) -> None:
    def defaults(label, df):
        return {
            "kind": "box",
            "x": "treatment",
            "y": "value",
            "color": _first(df, "patient", "patient_tumor"),
            "facet": _first(df, "channel"),
        }

    dataset_section(
        registry()["intensity"],
        "intensity",
        filters,
        "1.EDA/scripts/15.calculate_intensity_values.py",
        EDA_RESULTS / "intensity",
        defaults,
    )


def count_viability_section(filters: Filters) -> None:
    def defaults(label, df):
        numeric = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
        x = _first(df, "mean_cell_count", "mean_cells_per_organoid") or (
            numeric[0] if numeric else None
        )
        y = next((c for c in df.columns if "iab" in c.lower() and c in numeric), None)
        return {
            "kind": "scatter",
            "x": x,
            "y": y,
            "color": _first(df, "patient_tumor", "treatment"),
        }

    dataset_section(
        registry()["count_viability"],
        "count vs viability",
        filters,
        "1.EDA/scripts/10.calculate_count_viability_join.py",
        EDA_RESULTS / "count_viability",
        defaults,
    )


def consensus_heatmaps_section(filters: Filters) -> None:
    corr_dir = EDA_RESULTS / "correlation"
    pair_files = (
        sorted(corr_dir.glob("*_correlation_pairs.parquet"))
        if corr_dir.exists()
        else []
    )
    if not pair_files:
        missing_notice(
            "consensus heatmaps",
            "1.EDA/scripts/4.calculate_correlation_matrix.py",
            corr_dir,
        )
        return
    subpanel_background("corr_pairs")
    _pairs_heatmap(pair_files, filters, preset={"profile_type": "consensus"})


def correlation_viability_section(filters: Filters) -> None:
    def defaults(label, df):
        return {
            "kind": "scatter",
            "x": _first(df, "correlation"),
            "y": _first(df, "group1_group2_viability_diff"),
            "color": _first(df, "quadrant"),
        }

    dataset_section(
        registry()["correlation_viability"],
        "correlation vs viability",
        filters,
        "1.EDA/scripts/5a.find_correlation_pairs_for_montages.py",
        EDA_RESULTS / "correlation",
        defaults,
    )


AREA_VOLUME_PAIRS_PER_GROUP = 200


def _area_volume_pairs(area_path: str, volume_path: str) -> pd.DataFrame:
    """Random area/volume pairs within each patient x treatment.

    Area and volume come from separate pipelines with no shared organoid ID, so
    each pair is two independent draws from the same group (the same assumption
    as 18.plot_area_vs_volume.r). It shows the joint range of the two
    distributions, not a per-organoid relationship.
    """
    keys = ["patient_tumor", "treatment"]
    area = load_dataset(area_path)
    volume = load_dataset(volume_path)
    area = area[np.isfinite(area["area"])]
    volume = volume[np.isfinite(volume["volume"])]
    volume_groups = dict(tuple(volume.groupby(keys)))
    parts = []
    for key, area_group in area.groupby(keys):
        volume_group = volume_groups.get(key)
        if volume_group is None:
            continue
        n = min(len(area_group), len(volume_group), AREA_VOLUME_PAIRS_PER_GROUP)
        parts.append(
            pd.DataFrame(
                {
                    "patient_tumor": key[0],
                    "treatment": key[1],
                    "area": area_group["area"].sample(n, random_state=0).to_numpy(),
                    "volume": volume_group["volume"]
                    .sample(n, random_state=0)
                    .to_numpy(),
                }
            )
        )
    if not parts:
        return pd.DataFrame(columns=[*keys, "area", "volume"])
    return pd.concat(parts, ignore_index=True)


def area_vs_volume_section(filters: Filters) -> None:
    raw = registry()["area_vs_volume"]
    area, volume = raw.get("area_2D_organoid_raw"), raw.get("volume_3D_organoid_raw")
    if not area or not volume:
        missing_notice(
            "area vs volume",
            "1.EDA/scripts/17.calculate_area_volume_by_patient_treatment.py",
            EDA_RESULTS / "area_vs_volume",
        )
        return
    pairs = _area_volume_pairs(str(area.path), str(volume.path))
    if pairs.empty:
        st.warning("No patient x treatment group has both area and volume.")
        return
    st.caption(
        "Area and volume have no shared organoid ID, so each point pairs a random "
        "area and a random volume from the same patient x treatment (up to "
        f"{AREA_VOLUME_PAIRS_PER_GROUP} per group). It shows the joint range of the "
        "two distributions, not a per-organoid relationship."
    )
    explorer(
        pairs,
        "area_vs_volume",
        filters,
        kinds=["scatter", "histogram", "heatmap"],
        defaults={
            "kind": "scatter",
            "x": "area",
            "y": "volume",
            "color": "treatment",
            "facet": "patient_tumor",
        },
        title="Area (2D) vs volume (3D), randomly paired",
        free_facet_axes=True,
    )


def _count_measure_by_condition(
    counts_path: str, measure_path: str, measure: str
) -> pd.DataFrame:
    """One row per patient x treatment x dose: mean cells per organoid (over
    wells) against the mean organoid ``measure``. The tables share no organoid
    ID, so this is condition-level rather than per organoid."""
    keys = ["patient_tumor", "treatment", "dose"]
    counts = load_dataset(counts_path)
    cells = counts.groupby(keys, as_index=False).agg(
        cells_per_organoid=("mean_cells_per_organoid", "mean"),
        n_wells=("well", "nunique"),
    )
    values = load_dataset(measure_path)
    values = values[np.isfinite(values[measure])]
    values = values.groupby(keys, as_index=False).agg(
        **{f"mean_{measure}": (measure, "mean"), "n_organoids": (measure, "size")}
    )
    # doses are int in one table and float/Int64 in the other; compare as text
    for frame in (cells, values):
        frame["dose"] = frame["dose"].astype(str)
    return cells.merge(values, on=keys, how="inner")


def volume_area_vs_count_section(filters: Filters) -> None:
    counts = registry()["cell_counts"].get("organoid_cell_counts")
    raw = registry()["area_vs_volume"]
    measure = st.radio(
        "Measure",
        ["Volume (3D)", "Area (2D)"],
        horizontal=True,
        key="vol_area_count_measure",
    )
    col = "volume" if measure.startswith("Volume") else "area"
    source = raw.get(
        "volume_3D_organoid_raw" if col == "volume" else "area_2D_organoid_raw"
    )
    if not counts or not source:
        missing_notice(
            "volume/area vs count",
            "1.EDA/scripts/7.generate_cell_counts.py and "
            "1.EDA/scripts/17.calculate_area_volume_by_patient_treatment.py",
            EDA_RESULTS,
        )
        return
    conditions = _count_measure_by_condition(str(counts.path), str(source.path), col)
    if conditions.empty:
        st.warning("No patient x treatment x dose has both counts and this measure.")
        return
    st.caption(
        "One point per patient x treatment x dose. Cells per organoid is the mean "
        "over wells; the measure is the mean over organoids. The two tables share "
        "no organoid ID, so this is not a per-organoid relationship."
    )
    explorer(
        conditions,
        f"vol_area_count_{col}",
        filters,
        kinds=["scatter", "histogram", "heatmap"],
        defaults={
            "kind": "scatter",
            "x": "cells_per_organoid",
            "y": f"mean_{col}",
            "color": "treatment",
            "facet": "patient_tumor",
        },
        title=f"Mean {col} vs mean cells per organoid",
        free_facet_axes=True,
    )


EDA_SECTIONS = {
    "UMAP": umap_section,
    "PCA": pca_section,
    "Correlation heatmaps": correlation_section,
    "Consensus heatmaps": consensus_heatmaps_section,
    "Correlation vs viability": correlation_viability_section,
    "Cell counts": cell_counts_section,
    "Area & volume": area_volume_section,
    "Area vs volume": area_vs_volume_section,
    "Neighbors": neighbors_section,
    "Intensity": intensity_section,
    "Count vs viability": count_viability_section,
    "Volume & area vs count": volume_area_vs_count_section,
}


# ---------------------------------------------------------------------------
# 3.viability_prediction_models
# ---------------------------------------------------------------------------
def _viability_dataset(name: str) -> Dataset | None:
    return registry()["viability_models"].get(name)


def model_performance_section(filters: Filters) -> None:
    datasets = {
        k: v
        for k, v in registry()["viability_models"].items()
        if k in ("combined_fold_metrics", "combined_summary_metrics")
    }

    def defaults(label, df):
        return {
            "kind": "box" if "fold" in label else "bar",
            "x": "profile_type",
            "y": "R2",
            "color": _first(df, "shuffle_status"),
            "facet": _first(df, "split_method"),
        }

    dataset_section(
        datasets,
        "model performance",
        filters,
        "3.viability_prediction_models/scripts/1.viability_prediction.py",
        "3.viability_prediction_models/model_results",
        defaults,
    )


def predicted_vs_actual_section(filters: Filters) -> None:
    ds = _viability_dataset("combined_predicted_viabilities")
    if ds is None:
        missing_notice(
            "predicted viabilities",
            "3.viability_prediction_models/scripts/1.viability_prediction.py",
            "3.viability_prediction_models/model_results",
        )
        return
    # the file has ~23k feature columns; read only the metadata and target columns
    keep = tuple(
        c
        for c in parquet_columns(str(ds.path))
        if c.startswith("Metadata_")
        or c in ("Actual_Viability", "Predicted_Viability", "min_max_viability")
    )
    df = load_dataset(str(ds.path), keep)
    explorer(
        df,
        "vpred",
        filters,
        kinds=["scatter", "box", "violin", "histogram", "heatmap"],
        defaults={
            "kind": "scatter",
            "x": "Actual_Viability",
            "y": "Predicted_Viability",
            "color": "patient_tumor",
            "facet": _first(df, "split_method"),
        },
        title="Predicted vs actual viability",
    )


def feature_importance_section(filters: Filters) -> None:
    ds = _viability_dataset("combined_feature_importances")
    if ds is None:
        missing_notice(
            "feature importances",
            "3.viability_prediction_models/scripts/1.viability_prediction.py",
            "3.viability_prediction_models/model_results",
        )
        return
    df = load_dataset(str(ds.path))
    df = apply_global_filters(df, filters)
    c1, c2, c3, c4 = st.columns(4)
    selectors = {}
    for col, name in zip(
        (c1, c2, c3, c4),
        ("profile_type", "split_method", "shuffle_status", "image_mode"),
    ):
        if name in df.columns:
            options = sorted(df[name].astype(str).unique())
            selectors[name] = col.multiselect(
                name,
                options,
                default=options[:1] if name == "profile_type" else options,
                key=f"fi_{name}",
                format_func=humanize_label,
            )
    for name, values in selectors.items():
        if values:
            df = df[df[name].astype(str).isin(values)]
    df = local_subset(df, "fi")
    if df.empty:
        st.warning("No rows left after subsetting.")
        return
    c5, c6 = st.columns(2)
    top_n = c5.slider("Top N features", 5, 100, 25, key="fi_topn")
    color = c6.selectbox(
        "Color by",
        [NONE, "split_method", "shuffle_status", "profile_type", "held_out_group"],
        key="fi_color",
        format_func=lambda v: v if v == NONE else humanize_label(v),
    )
    color = None if color == NONE or color not in df.columns else color
    top = df.groupby("feature")["importance"].mean().nlargest(top_n).index
    grouped = df[df["feature"].isin(top)]
    group = ["feature"] + ([color] if color else [])
    grouped = grouped.groupby(group, observed=True)["importance"].mean().reset_index()
    fig = px.bar(
        grouped,
        x="importance",
        y="feature",
        color=color,
        barmode="group",
        orientation="h",
        category_orders={"feature": list(top)},
        title=f"Top {top_n} features by mean importance",
    )
    fig.update_layout(height=max(450, 22 * top_n), template="plotly_white")
    st.plotly_chart(fig, width="stretch", key="fi_chart")
    png_download(fig, "fi", "feature_importances")


VIABILITY_SECTIONS = {
    "Model performance": model_performance_section,
    "Predicted vs actual": predicted_vs_actual_section,
    "Feature importances": feature_importance_section,
}


# ---------------------------------------------------------------------------
# 4.linear_modeling
# ---------------------------------------------------------------------------
def _normalize_lm_df(df: pd.DataFrame) -> pd.DataFrame:
    """Column/value fixups shared by every linear-modeling results table.

    The "_technical_model" files keep their raw Metadata_-prefixed schema, so
    ``canonicalize_columns`` maps their patient column to "patient_tumor" (not
    "patient") and their treatment-term rows still read the raw column name.
    """
    if "patient" not in df.columns and "patient_tumor" in df.columns:
        df = df.rename(columns={"patient_tumor": "patient"})
    if "term" in df.columns:
        df["term"] = df["term"].replace({"Metadata_Experiment_Treatment": "treatment"})
    if "tumor_type" not in df.columns and "patient" in df.columns:
        df["tumor_type"] = df["patient"].map(TUMOR_TYPE_LOOKUP)
    return df


def _linear_model_data(key: str, extra_defaults=None):
    datasets = registry()["linear_modeling"]
    if not datasets:
        missing_notice(
            "linear-model results",
            "4.linear_modeling/scripts/0.linear_modeling.py",
            "4.linear_modeling/results/linear_modeling",
        )
        return None, None
    label = st.selectbox(
        "Model results",
        list(datasets),
        key=f"{key}_dataset",
        format_func=humanize_label,
    )
    df = _normalize_lm_df(load_dataset(str(datasets[label].path)))
    if "pvalue_fdr" in df.columns:
        df["neg_log10_pvalue_fdr"] = -np.log10(df["pvalue_fdr"].clip(lower=1e-300))
        df["significant_fdr_0.05"] = np.where(
            df["pvalue_fdr"] < 0.05, "significant", "n.s."
        )
    if "coefficient" in df.columns:
        df["abs_coefficient"] = df["coefficient"].abs()
    return label, df


def lm_volcano_section(filters: Filters) -> None:
    label, df = _linear_model_data("lmvol")
    if df is None:
        return
    explorer(
        df,
        f"lmvol_{label}",
        filters,
        kinds=["scatter", "heatmap", "histogram"],
        defaults={
            "kind": "scatter",
            "x": "coefficient",
            "y": "neg_log10_pvalue_fdr",
            "color": "significant_fdr_0.05",
            "facet": _first(df, "term"),
        },
        title=f"{label}: volcano",
    )


def lm_effects_section(filters: Filters) -> None:
    label, df = _linear_model_data("lmeff")
    if df is None:
        return
    explorer(
        df,
        f"lmeff_{label}",
        filters,
        kinds=["box", "violin", "bar", "histogram", "scatter"],
        defaults={
            "kind": "box",
            "x": _first(df, "therapeutic_category", "treatment"),
            "y": "coefficient",
            "color": _first(df, "term"),
            "facet": _first(df, "patient"),
        },
        title=f"{label}: effect sizes",
    )


def lm_fit_section(filters: Filters) -> None:
    label, df = _linear_model_data("lmfit")
    if df is None:
        return
    explorer(
        df,
        f"lmfit_{label}",
        filters,
        kinds=["histogram", "box", "violin", "bar", "scatter"],
        defaults={
            "kind": "histogram",
            "x": "rsquared_adj",
            "color": _first(df, "Feature_type"),
            "facet": _first(df, "Compartment"),
        },
        title=f"{label}: model fit",
    )


# ported from 6.plot_variate_importance.r's `scopes` list: scope name -> its group columns
UPSET_SCOPES = {
    "all_models": [],
    "per_patient_treatment": ["patient", "treatment"],
    "per_patient": ["patient"],
    "per_treatment": ["treatment"],
    "per_treatment_tumor_type": ["treatment", "tumor_type"],
    "per_tumor_type": ["tumor_type"],
}

# profile x model_set -> the linear_modeling dataset that holds it
UPSET_PROFILE_DATASETS = {
    ("organoid", "original"): "organoid_norm",
    ("organoid", "technical"): "organoid_norm_technical_model",
    ("sc", "original"): "sc_norm",
    ("sc", "technical"): "sc_norm_technical_model",
}


def _build_membership(
    hits: pd.DataFrame, terms: list[str], index_cols: list[str]
) -> pd.DataFrame:
    """One row per (index_cols..., feature), one boolean column per term: True
    when any model at that index is a hit for that term. Ported from
    ``5.calculate_variate_importance.py``'s ``build_membership()``."""
    return (
        hits.loc[hits["term"].isin(terms), index_cols + ["term", "hit"]]
        .pivot_table(index=index_cols, columns="term", values="hit", aggfunc="max")
        .reindex(columns=terms)
        .fillna(False)
        .astype(bool)
    )


def _upset_combinations(membership: pd.DataFrame, terms: list[str]) -> pd.DataFrame:
    """Features per exact term combination, largest first. Ported from
    ``5.calculate_variate_importance.py``'s ``upset_combinations()``."""
    combos = (
        membership.groupby(terms, observed=True)
        .size()
        .rename("n_features")
        .reset_index()
    )
    combos = (
        combos.loc[combos[terms].any(axis=1)]
        .sort_values("n_features", ascending=False)
        .reset_index(drop=True)
    )
    if combos.empty:
        return combos.assign(combination=[], treatment_specific=[])
    combos["combination"] = combos[terms].apply(
        lambda r: " + ".join(t for t in terms if r[t]), axis=1
    )
    other = [t for t in terms if t != "treatment"]
    combos["treatment_specific"] = (
        combos["treatment"] & ~combos[other].any(axis=1)
        if "treatment" in terms and other
        else pd.Series(False, index=combos.index)
    )
    return combos.rename(columns={t: f"in_{t}" for t in terms})


def _set_sizes(membership: pd.DataFrame, terms: list[str]) -> pd.DataFrame:
    return pd.DataFrame({"term": terms, "set_size": membership[terms].sum().to_numpy()})


def lm_upset_section(filters: Filters) -> None:
    """UpSet plots of which variates (model terms) a feature is a 'hit' for,
    computed on the fly (ported from ``5.calculate_variate_importance.py`` /
    ``6.plot_variate_importance.r``'s ``plot_upset()``): a feature is a hit for
    a term when a model has ``pvalue_fdr < threshold`` and ``coefficient >
    min`` for it. Below the plot, an overlap matrix shows how many features
    a chosen combination shares between every pair of groups in the scope
    (e.g. MPNST vs pNF for the tumor-type scope).
    """
    c1, c2 = st.columns(2)
    profile = c1.radio(
        "Profile", ["organoid", "sc"], horizontal=True, key="upset_profile"
    )
    model_set = c2.radio(
        "Model set", ["original", "technical"], horizontal=True, key="upset_model_set"
    )
    subpanel_background("lm_upset_profile")
    dataset_key = UPSET_PROFILE_DATASETS[(profile, model_set)]
    ds = registry()["linear_modeling"].get(dataset_key)
    if ds is None:
        missing_notice(
            f"{profile} ({model_set}) linear-model results",
            "4.linear_modeling/scripts/0.linear_modeling.py",
            "4.linear_modeling/results/linear_modeling",
        )
        return
    df = _normalize_lm_df(load_dataset(str(ds.path)))
    df = apply_global_filters(df, filters)
    if df.empty:
        st.warning("No rows left after subsetting.")
        return
    feature_meta = df.drop_duplicates("feature")[
        ["feature", "Feature_type", "Channel", "Compartment"]
    ]

    c3, c4 = st.columns(2)
    fdr_max = c3.slider("FDR threshold", 0.001, 0.25, 0.05, key="upset_fdr")
    coef_min = c4.slider(
        "Min coefficient (increases only, matches 5.calculate_variate_importance.py)",
        0.0,
        1.0,
        0.1,
        key="upset_coef",
    )
    hits = df.assign(hit=(df["pvalue_fdr"] < fdr_max) & (df["coefficient"] > coef_min))
    terms = [t for t in TERM_ORDER if t in hits["term"].unique()]

    scope = st.selectbox(
        "Scope", list(UPSET_SCOPES), key="upset_scope", format_func=humanize_label
    )
    group_cols = UPSET_SCOPES[scope]
    index_cols = group_cols + ["feature"]
    membership_all = _build_membership(hits, terms, index_cols)
    if membership_all.empty:
        st.warning("No hits at the current thresholds.")
        return

    group_label = ""
    if group_cols:
        cols = st.columns(len(group_cols))
        selection = {}
        for col_widget, name in zip(cols, group_cols):
            values = sorted(membership_all.index.get_level_values(name).unique())
            selection[name] = col_widget.selectbox(name, values, key=f"upset_{name}")
        mask = pd.Series(True, index=membership_all.index)
        for name, value in selection.items():
            mask &= membership_all.index.get_level_values(name) == value
        membership = membership_all[mask]
        membership.index = membership.index.get_level_values("feature")
        group_label = " | ".join(f"{k}={v}" for k, v in selection.items())
    else:
        membership = membership_all

    combos = _upset_combinations(membership, terms)
    sizes = _set_sizes(membership, terms)
    if combos.empty:
        st.warning("No term combinations for this selection.")
        return

    c5, c6 = st.columns([2, 1])
    top_n = c5.slider("Top N combinations", 5, 40, 25, key="upset_topn")
    term_group = c6.radio(
        "Variate group",
        ["Individual variates", "Grouped (Biological vs Technical)"],
        horizontal=True,
        key="upset_term_group",
    )
    if term_group == "Grouped (Biological vs Technical)":
        groups = {
            "Biological": [t for t in BIOLOGICAL_TERMS if t in terms],
            "Technical": [t for t in TECHNICAL_TERMS if t in terms],
        }
        groups = {name: g_terms for name, g_terms in groups.items() if g_terms}
        if len(groups) < 2:
            st.warning(
                "This model set doesn't have both biological and technical terms."
            )
            return
        combos, sizes = regroup_combos_by_identity(combos, groups)
        term_order = list(groups)
    else:
        term_order = terms

    n_total = len(membership)
    n_none = int((~membership[terms].any(axis=1)).sum())
    st.caption(f"{n_total:,} features ({n_none:,} not a hit for any term)")

    title = f"{profile} {model_set} model, {scope}" + (
        f" ({group_label})" if group_label else ""
    )
    fig = upset_plot(sizes, combos, term_order, top_n=top_n, title=title)
    if fig is None:
        st.warning("No term combinations to show.")
        return
    st.plotly_chart(fig, width="stretch", key="upset_chart")
    png_download(fig, "upset", f"upset_{profile}_{model_set}_{scope}")

    if not group_cols:
        return
    st.divider()
    st.caption(f"Feature overlap across {' x '.join(group_cols)} groups in this scope")
    group_keys = membership_all.index.to_frame(index=False)[
        group_cols
    ].drop_duplicates()
    group_keys["_label"] = (
        group_keys[group_cols[0]].astype(str)
        if len(group_cols) == 1
        else group_keys[group_cols].astype(str).agg(" | ".join, axis=1)
    )
    n_groups = len(group_keys)
    if n_groups < 2:
        st.info("Only one group in this scope; nothing to compare.")
        return
    if n_groups > 60:
        st.caption(
            f"{n_groups} groups -- rendering the full pairwise matrix, cell counts "
            "hidden for legibility (hover to read a value)."
        )

    combo_options = combos["combination"].tolist()
    chosen_combo = st.selectbox(
        "Combination to compare across groups", combo_options, key="upset_overlap_combo"
    )
    combo_terms = term_order if term_group == "Individual variates" else list(groups)
    pattern_row = combos.loc[combos["combination"] == chosen_combo].iloc[0]
    match_mask = pd.Series(True, index=membership_all.index)
    for t in combo_terms:
        want = bool(pattern_row[f"in_{t}"])
        col = (
            membership_all[t]
            if term_group == "Individual variates"
            else membership_all[
                [
                    c
                    for c in (
                        BIOLOGICAL_TERMS if t == "Biological" else TECHNICAL_TERMS
                    )
                    if c in terms
                ]
            ].any(axis=1)
        )
        match_mask &= col == want
    matched = membership_all[match_mask].index.to_frame(index=False)

    sets_by_group = {}
    for _, row in group_keys.iterrows():
        sub_mask = pd.Series(True, index=matched.index)
        for c in group_cols:
            sub_mask &= matched[c] == row[c]
        sets_by_group[row["_label"]] = set(matched.loc[sub_mask, "feature"])
    labels = group_keys["_label"].tolist()
    overlap = pd.DataFrame(
        [[len(sets_by_group[a] & sets_by_group[b]) for b in labels] for a in labels],
        index=labels,
        columns=labels,
    )
    overlap_fig = px.imshow(
        overlap,
        text_auto=n_groups <= 60,
        color_continuous_scale="Viridis",
        aspect="auto",
        labels=dict(color="shared features"),
        title=f"Shared features for '{chosen_combo}' across {' x '.join(group_cols)}",
    )
    overlap_fig.update_layout(
        template="plotly_white", height=max(400, 40 * n_groups + 150)
    )
    event = st.plotly_chart(
        overlap_fig,
        width="stretch",
        key="upset_overlap_chart",
        on_select="rerun",
        selection_mode="points",
    )
    png_download(overlap_fig, "upset_overlap", f"overlap_{scope}_{chosen_combo}")

    # a click on a cell updates the two pickers below; they also work by hand.
    # Widgets ignore `index=` once their key already holds a value, so a new
    # click has to be pushed into session_state *before* the widgets are
    # created below -- but only once per click, or it would fight the user's
    # own dropdown choice on every later, unrelated rerun.
    for key in ("upset_overlap_row", "upset_overlap_col"):
        if key in st.session_state and st.session_state[key] not in labels:
            del st.session_state[key]  # stale value from a previous scope

    points = event.selection.points if event else []
    clicked = (points[0]["y"], points[0]["x"]) if points else None
    if clicked and (clicked[0] not in labels or clicked[1] not in labels):
        clicked = None  # stale selection from a previous scope/combination
    if clicked and st.session_state.get("upset_overlap_last_click") != clicked:
        st.session_state["upset_overlap_row"] = clicked[0]
        st.session_state["upset_overlap_col"] = clicked[1]
        st.session_state["upset_overlap_last_click"] = clicked

    c7, c8 = st.columns(2)
    row_label = c7.selectbox("Row group", labels, key="upset_overlap_row")
    col_label = c8.selectbox("Column group", labels, key="upset_overlap_col")
    clicked_features = (
        sets_by_group[row_label]
        if row_label == col_label
        else sets_by_group[row_label] & sets_by_group[col_label]
    )
    header = (
        f"**{len(clicked_features):,} features** in **{row_label}**"
        if row_label == col_label
        else f"**{len(clicked_features):,} features** shared between "
        f"**{row_label}** and **{col_label}**"
    )
    st.markdown(header)
    if not clicked_features:
        return
    feat_df = feature_meta[feature_meta["feature"].isin(clicked_features)]
    with st.expander(f"Feature names ({len(feat_df):,})"):
        st.dataframe(feat_df, width="stretch", hide_index=True)

    melted = feat_df.melt(
        id_vars="feature",
        value_vars=["Compartment", "Channel", "Feature_type"],
        var_name="attribute",
        value_name="value",
    )
    melted["value"] = melted["value"].fillna("Other")
    counts = (
        melted.groupby(["attribute", "value"], observed=True)
        .size()
        .reset_index(name="count")
    )
    detail_fig = px.bar(
        counts,
        x="value",
        y="count",
        color="attribute",
        facet_col="attribute",
        title="Selected features by Compartment / Channel / Feature type",
    )
    detail_fig.update_xaxes(matches=None, showticklabels=True)
    detail_fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    detail_fig.update_layout(template="plotly_white", showlegend=False, height=400)
    st.plotly_chart(detail_fig, width="stretch", key="upset_overlap_detail_chart")
    png_download(
        detail_fig, "upset_overlap_detail", f"overlap_detail_{row_label}_{col_label}"
    )


def _annotation_legends(tracks: list[tuple[str, list[str], str]]) -> None:
    """Small colored-swatch legends for the heatmap's annotation strips (Plotly
    heatmaps have no native per-category legend for a discrete color track)."""
    cols = st.columns(len(tracks))
    for col, (title, categories, palette_key) in zip(cols, tracks):
        colors = palette_for(palette_key, categories) or {}
        swatches = "".join(
            f'<span style="display:inline-block;width:10px;height:10px;'
            f"background:{colors.get(c, '#cccccc')};margin-right:4px;"
            f'border:1px solid #999;"></span>{c}<br>'
            for c in categories
        )
        with col:
            st.markdown(f"**{title}**<br>{swatches}", unsafe_allow_html=True)


def lm_model_variates_section(filters: Filters) -> None:
    """Where is treatment significant on its own vs alongside covariates?

    Feature x (patient, treatment) yes/no heatmap -- ported from
    ``6.plot_variate_importance.r``'s ``cooccurrence_heatmap()`` -- with
    Patient, Tumor type and Treatment column annotations and a Feature type
    row annotation. A toggle switches the criterion between ``treatment``
    being a hit on its own vs together with at least one biological covariate
    (cell count / organoid count terms);
    ``hit = pvalue_fdr < threshold & coefficient > min``, matching
    ``5.calculate_variate_importance.py``. A single-model drill-down table
    (its 'unique variate signature') is below the heatmap.
    """
    label, df = _linear_model_data("lmvar")
    if df is None:
        return
    df = apply_global_filters(df, filters)
    if df.empty:
        st.warning("No rows left after subsetting.")
        return

    c1, c2 = st.columns(2)
    fdr_max = c1.slider("FDR threshold", 0.001, 0.25, 0.05, key="lmvar_fdr")
    coef_min = c2.slider(
        "Min coefficient (increases only, matches 5.calculate_variate_importance.py)",
        0.0,
        1.0,
        0.1,
        key="lmvar_coef",
    )
    # the covariate group is fixed to biological: cell counts (treatment excluded)
    group_terms = [
        t for t in BIOLOGICAL_TERMS if t != "treatment" and t in df["term"].unique()
    ]
    if not group_terms:
        st.warning("No biological covariates in this dataset.")
        return
    criterion = st.radio(
        "Criterion",
        ["Treatment only", "Treatment + biological"],
        horizontal=True,
        key="lmvar_criterion",
    )

    hits = df.assign(hit=(df["pvalue_fdr"] < fdr_max) & (df["coefficient"] > coef_min))
    model_keys = ["patient", "treatment", "feature"]
    treatment_hit = (
        hits.loc[hits["term"] == "treatment"]
        .set_index(model_keys)["hit"]
        .rename("treatment_hit")
    )
    group_hit = (
        hits.loc[hits["term"].isin(group_terms)]
        .groupby(model_keys)["hit"]
        .any()
        .rename("group_hit")
    )
    models = pd.concat([treatment_hit, group_hit], axis=1).fillna(False).reset_index()
    # "Treatment only" excludes models where the covariate group is also a hit;
    # "Treatment + group" is every treatment hit, whether or not the group is
    # also significant -- a superset of "Treatment only", not an intersection
    models["treatment_only"] = models["treatment_hit"] & ~models["group_hit"]
    selected_col = (
        "treatment_only" if criterion == "Treatment only" else "treatment_hit"
    )

    model_meta = df.drop_duplicates(model_keys)[
        model_keys + ["drug", "Feature_type", "Channel", "Compartment", "tumor_type"]
    ]
    models = models.merge(model_meta, on=model_keys, how="left")
    sig = models.loc[models[selected_col]]
    n_hits = len(sig)
    st.caption(
        f"{n_hits:,} of {len(models):,} (patient, treatment, feature) models are "
        f"significant under **{criterion}**"
    )
    if sig.empty:
        st.warning("No models meet this criterion at the current thresholds.")
        return
    if n_hits > 20_000:
        st.info(
            f"{n_hits:,} significant models is a lot to render as a heatmap -- "
            "raise the FDR/coefficient thresholds to narrow it down."
        )
        return

    with st.expander("Heatmap: color bars & clustering", expanded=False):
        c6, c7 = st.columns(2)
        cluster_rows = c6.checkbox(
            "Cluster rows (features)", value=True, key="lmvar_cluster_rows"
        )
        cluster_cols = c6.checkbox(
            "Cluster columns (patient x treatment)",
            value=True,
            key="lmvar_cluster_cols",
        )
        row_tracks = c7.multiselect(
            "Row color bars",
            ROW_TRACK_NAMES,
            default=ROW_TRACK_NAMES,
            key="lmvar_row_tracks",
        )
        col_tracks = c7.multiselect(
            "Column color bars",
            COL_TRACK_NAMES,
            default=COL_TRACK_NAMES,
            key="lmvar_col_tracks",
        )

    fig = annotated_significance_heatmap(
        sig,
        criterion,
        legend_title=criterion,
        cluster_rows=cluster_rows,
        cluster_cols=cluster_cols,
        row_tracks=row_tracks,
        col_tracks=col_tracks,
    )
    if fig is None:
        st.warning("Nothing to plot.")
    else:
        st.plotly_chart(fig, width="content", key="lmvar_heatmap")
        png_download(fig, "lmvar_heatmap", f"{label}_{criterion.replace(' ', '_')}")
        all_legends = {
            "Patient": ("Patient", sorted(sig["patient"].unique()), "patient"),
            "Tumor type": (
                "Tumor type",
                sorted(sig["tumor_type"].dropna().unique()),
                "tumor_type",
            ),
            "Treatment": ("Treatment", sorted(sig["drug"].unique()), "treatment"),
            "Feature type": (
                "Feature type",
                sorted(sig["Feature_type"].fillna("Other").unique()),
                "feature_type",
            ),
            "Channel": (
                "Channel",
                sorted(sig["Channel"].fillna("Other").unique()),
                "channel",
            ),
            "Compartment": (
                "Compartment",
                sorted(sig["Compartment"].fillna("Other").unique()),
                "compartment",
            ),
        }
        shown = [all_legends[t] for t in col_tracks + row_tracks if t in all_legends]
        if shown:
            _annotation_legends(shown)

    st.divider()
    st.caption("Drill into one model's full term signature")
    subpanel_background("lm_variates_drilldown")
    c3, c4, c5 = st.columns(3)
    patient = c3.selectbox(
        "Patient", sorted(df["patient"].unique()), key="lmvar_patient"
    )
    df_p = df[df["patient"] == patient]
    treatment = c4.selectbox(
        "Treatment", sorted(df_p["treatment"].unique()), key="lmvar_treatment"
    )
    df_pt = df_p[df_p["treatment"] == treatment]
    feature = c5.selectbox(
        "Feature", sorted(df_pt["feature"].unique()), key="lmvar_feature"
    )
    model = df_pt[df_pt["feature"] == feature].sort_values("term")
    if model.empty:
        st.warning("No rows for this patient / treatment / feature combination.")
        return
    model = model.assign(
        hit=(model["pvalue_fdr"] < fdr_max) & (model["coefficient"] > coef_min)
    )
    unique_variates = model.loc[model["hit"], "term"].tolist()
    model["hit"] = model["hit"].astype(str)
    st.markdown(
        "**Unique variates (significant terms) for this model:** "
        + (
            ", ".join(f"`{t}`" for t in unique_variates)
            if unique_variates
            else "_none_"
        )
    )
    st.dataframe(
        model[["term", "coefficient", "pvalue_fdr", "hit", "rsquared", "rsquared_adj"]],
        width="stretch",
        hide_index=True,
    )


LINEAR_MODELING_SECTIONS = {
    "Volcano (effects vs significance)": lm_volcano_section,
    "Effect sizes": lm_effects_section,
    "Model fit": lm_fit_section,
    "UpSet (variate combinations)": lm_upset_section,
    "Unique variates per model": lm_model_variates_section,
}


# ---------------------------------------------------------------------------
# 5.differential_analysis (each treatment vs its patient's DMSO control)
# Plot types follow the R figures of PRs 55-58 (scripts 1, 3, 5, 7, 9).
# ---------------------------------------------------------------------------
COMPARTMENTS = {"organoid": "Organoid", "single_cell": "Single cell"}
MORPHOLOGY_COLUMN = {
    "organoid": "organoid_mean_abs_log2fc",
    "single_cell": "single_cell_mean_abs_log2fc",
}
TUMOR_TYPES = ["cNF", "pNF", "MPNST", "Other"]
MEK_MOA = "MEK1/2 inhibitor"
MEK_LABEL = "MEK inhibitor"
Q_THRESHOLD = 0.05
RED = "#d73027"
SEQUENTIAL = [
    [0, "#ffffff"],
    [1, "#54278f"],
]  # morphology / counts (R: white -> purple)
DIVERGING = [[0, "#b2182b"], [0.5, "#ffffff"], [1, "#2166ac"]]  # viability (R)
LEVEL_LABELS = {
    "within_patient": "Within patient",
    "within_tumor_type": "Within tumor type",
    "across_patients": "Across patients",
    "across_tumor_types_shared": "Across tumor types (shared)",
    "across_tumor_types_differ": "Across tumor types (differ)",
}
POOLED_GROUP_LABELS = {
    "across_patients": "All patients",
    "across_tumor_types_shared": "Tumor type means",
    "across_tumor_types_differ": "Between tumor types",
}
CONTRASTS = ["MEKi vs DMSO", "MEKi vs Staurosporine", "MEKi vs Digoxin"]


def _differential(key: str) -> pd.DataFrame | None:
    """Load one 5.differential_analysis table, or say which script writes it."""
    ds = registry()["differential"].get(key)
    if ds is None:
        missing_notice(
            humanize_label(key),
            DIFFERENTIAL_PRODUCED_BY[key],
            DIFFERENTIAL_RESULTS / f"{key}.parquet",
        )
        return None
    return load_dataset(str(ds.path))


def _with_treatment_label(df: pd.DataFrame) -> pd.DataFrame:
    """``treatment dose unit`` (e.g. ``Binimetinib 10 uM``), so each dose gets its
    own column in the heatmaps instead of being pooled under the drug name."""
    return df.assign(
        treatment_label=(
            df["treatment"].astype(str)
            + " "
            + df["dose"].astype(str)
            + " "
            + df["dose_unit"].astype(str)
        )
    )


def _treatment_order(df: pd.DataFrame) -> list[str]:
    """Treatment labels grouped by mechanism, then drug, then dose."""
    return (
        df.drop_duplicates("treatment_label")
        .sort_values(["moa", "treatment", "dose"])["treatment_label"]
        .tolist()
    )


def _tile_panels(
    df: pd.DataFrame,
    *,
    panel: str,
    row: str,
    col: str,
    value: str,
    panels: list[str],
    cols: list[str],
    title: str,
    colorscale,
    colorbar_title: str,
    zmin=None,
    zmax=None,
    zmid=None,
    label_values: bool = False,
    transform=None,
    tickvals=None,
    ticktext=None,
) -> go.Figure:
    """Tile heatmap with one panel per group, like R's
    ``facet_grid(group ~ ., scales = "free_y", space = "free_y")``. Missing tiles
    show the grey plot background (R's ``na.value = "grey85"``)."""
    row_order = {p: list(dict.fromkeys(df.loc[df[panel] == p, row])) for p in panels}
    fig = make_subplots(
        rows=len(panels),
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.06,
        row_heights=[len(row_order[p]) for p in panels],
        subplot_titles=panels,
    )
    for i, p in enumerate(panels, start=1):
        grid = (
            df.loc[df[panel] == p]
            .pivot_table(
                index=row, columns=col, values=value, aggfunc="first", dropna=False
            )
            .reindex(index=row_order[p], columns=cols)
            .to_numpy(dtype=float)
        )
        text = None
        if label_values:
            text = [[("" if np.isnan(v) else f"{v:.0f}") for v in r] for r in grid]
        fig.add_trace(
            go.Heatmap(
                z=transform(grid) if transform else grid,
                x=cols,
                y=row_order[p],
                customdata=grid,
                colorscale=colorscale,
                zmin=zmin,
                zmax=zmax,
                zmid=zmid,
                text=text,
                texttemplate="%{text}" if label_values else None,
                xgap=1,
                ygap=1,
                showscale=i == 1,
                colorbar=dict(
                    title=colorbar_title, tickvals=tickvals, ticktext=ticktext
                ),
                hovertemplate="%{y}<br>%{x}<br>%{customdata:.3g}<extra></extra>",
            ),
            row=i,
            col=1,
        )
    fig.update_yaxes(autorange="reversed")
    fig.update_xaxes(tickangle=-45)
    fig.update_layout(
        title=title,
        template="plotly_white",
        plot_bgcolor="#d9d9d9",
        showlegend=False,
        height=180 + 24 * sum(len(row_order[p]) for p in panels) + 40 * len(panels),
    )
    return fig


def _show(fig: go.Figure, key: str, filename: str) -> None:
    st.plotly_chart(fig, width="stretch", key=f"{key}_chart")
    png_download(fig, key, filename)


# ---------------------------------------------------------------------------
# Viability heatmap (shared by both compartments; R script 3)
# ---------------------------------------------------------------------------
def viability_heatmap_section(filters: Filters) -> None:
    per_patient = _differential("viability_log2fc_by_patient_tumor_type")
    per_type = _differential("viability_log2fc_by_tumor_type")
    if per_patient is None or per_type is None:
        return
    choice = st.radio(
        "Patients",
        ["All patients", "Patients with viability data"],
        horizontal=True,
        key="viab_heat_choice",
    )
    patient = _with_treatment_label(apply_global_filters(per_patient, filters))
    tumor = _with_treatment_label(apply_global_filters(per_type, filters))
    if choice == "Patients with viability data":
        patient = patient.loc[patient["viability_percent"].notna()]
        tumor = tumor.loc[tumor["n_patients"] > 0]
    rows = pd.concat(
        [
            patient.assign(
                row=patient["patient_tumor"], value=patient["viability_percent"]
            ),
            tumor.assign(
                row=tumor["tumor_type"] + " mean",
                value=tumor["mean_viability_percent"],
            ),
        ],
        ignore_index=True,
    )
    if rows.empty:
        st.warning("No rows left after subsetting.")
        return
    rows = rows.assign(is_mean=rows["row"].str.endswith(" mean")).sort_values(
        ["is_mean", "row"]
    )
    fig = _tile_panels(
        rows,
        panel="tumor_type",
        row="row",
        col="treatment_label",
        value="value",
        panels=[t for t in TUMOR_TYPES if t in set(rows["tumor_type"])],
        cols=_treatment_order(rows),
        title=f"Viability (% of DMSO): {choice.lower()}",
        colorscale=DIVERGING,
        colorbar_title="Viability<br>(% of DMSO)",
        zmid=100,
        label_values=True,
    )
    _show(fig, "viab_heat", "viability_percent_heatmap")


# ---------------------------------------------------------------------------
# Per-compartment sections
# ---------------------------------------------------------------------------
def _morphology_heatmap(summary: pd.DataFrame, comp: str, filters: Filters) -> None:
    """Tile heatmap of morphology per patient and treatment (R script 1)."""
    df = _with_treatment_label(apply_global_filters(summary, filters))
    if df.empty:
        st.warning("No rows left after subsetting.")
        return
    df = df.assign(row=df["patient_tumor"], value=df[MORPHOLOGY_COLUMN[comp]])
    fig = _tile_panels(
        df.sort_values("patient_tumor"),
        panel="tumor_type",
        row="row",
        col="treatment_label",
        value="value",
        panels=[t for t in TUMOR_TYPES if t in set(df["tumor_type"])],
        cols=_treatment_order(df),
        title=f"{COMPARTMENTS[comp]}: morphology difference from DMSO",
        colorscale=SEQUENTIAL,
        colorbar_title="Mean |difference|<br>from DMSO",
        zmin=0,
    )
    _show(fig, f"{comp}_morph_heat", f"{comp}_morphology_heatmap")


def _morphology_vs_viability(
    summary: pd.DataFrame, comp: str, filters: Filters
) -> None:
    """Scatter per patient x treatment, then per-treatment means with labelled
    outliers (R script 1, two figures)."""
    df = _with_treatment_label(apply_global_filters(summary, filters)).assign(
        morph=lambda d: d[MORPHOLOGY_COLUMN[comp]],
        treatment_type=lambda d: np.where(d["moa"] == MEK_MOA, MEK_LABEL, "Other"),
    )
    pts = df.dropna(subset=["viability_log2fc"])
    if pts.empty:
        st.info("No viability measurements for this compartment.")
        return
    symbols = {"Other": "circle", MEK_LABEL: "triangle-up"}
    common = dict(
        x="viability_log2fc",
        y="morph",
        color="tumor_type",
        symbol="treatment_type",
        color_discrete_map=TUMOR_TYPE_PALETTE,
        symbol_map=symbols,
        category_orders={"tumor_type": TUMOR_TYPES},
        labels={
            "viability_log2fc": "Viability log2FC vs DMSO",
            "morph": "Morphology mean |difference| vs DMSO",
            "tumor_type": "Tumor type",
            "treatment_type": "Treatment",
        },
    )
    fig = px.scatter(
        pts,
        hover_data=["patient_tumor", "treatment_label"],
        title=f"{COMPARTMENTS[comp]}: one point per patient x treatment",
        **common,
    )
    fig.add_vline(x=0, line_dash="dash", line_color="grey")
    fig.update_layout(template="plotly_white", height=520)
    _show(fig, f"{comp}_morph_scatter", f"{comp}_morphology_vs_viability")

    means = (
        pts.groupby(["treatment_label", "tumor_type", "treatment_type"], observed=True)
        .agg(
            viability_log2fc=("viability_log2fc", "mean"),
            morph=("morph", "mean"),
        )
        .reset_index()
    )
    per_treatment = pts.groupby("treatment_label").agg(
        mean_viab=("viability_log2fc", "mean"), mean_morph=("morph", "mean")
    )
    # outliers: strongly reduced viability, or top three morphology effects
    outliers = per_treatment.index[
        (per_treatment["mean_viab"] < -0.5)
        | (per_treatment["mean_morph"].rank(ascending=False, method="first") <= 3)
    ]
    means["label"] = ""
    top = (
        means.loc[means["treatment_label"].isin(outliers)]
        .groupby("treatment_label")["morph"]
        .idxmax()
    )
    means.loc[top, "label"] = means.loc[top, "treatment_label"]
    fig = px.scatter(
        means,
        text="label",
        hover_data=["treatment_label"],
        title=f"{COMPARTMENTS[comp]}: one point per treatment x tumor type, mean over patients",
        **common,
    )
    fig.update_traces(textposition="top center", textfont_size=10)
    fig.add_vline(x=0, line_dash="dash", line_color="grey")
    fig.update_layout(template="plotly_white", height=520)
    _show(fig, f"{comp}_treat_scatter", f"{comp}_treatment_means_vs_viability")


def _plate_position(comp: str, filters: Filters) -> None:
    """Null histogram and per-feature scatter (R script 5), DMSO pair boxplots
    (R script 5), and the DMSO column table."""
    glob = _differential("plate_position_global_test")
    null = _differential("plate_position_null")
    feat = _differential("plate_position_feature_tests")
    pairs = _differential("dmso_column_pair_correlations")
    tests = _differential("dmso_column_tests")
    if any(x is None for x in (glob, null, feat, pairs, tests)):
        return
    g = glob.loc[glob["compartment"] == comp]
    if g.empty:
        st.info("No plate-position test for this compartment.")
        return
    observed = float(g["observed_rms_effect"].iloc[0])
    fig = px.histogram(
        null.loc[null["compartment"] == comp],
        x="null_rms_effect",
        nbins=60,
        labels={"null_rms_effect": "RMS outer - inner effect across features"},
        title=f"{COMPARTMENTS[comp]}: grey = inner/outer labels shuffled, red = observed "
        f"(p = {float(g['p'].iloc[0]):.3f})",
    )
    fig.update_traces(marker=dict(color="#b3b3b3", line=dict(color="white", width=0.5)))
    fig.add_vline(x=observed, line_color=RED, line_width=2)
    fig.update_layout(template="plotly_white", height=380, yaxis_title="Permutations")
    _show(fig, f"{comp}_pp_null", f"{comp}_plate_position_null")

    f = feat.loc[feat["compartment"] == comp].assign(
        neg_log10_p=lambda d: -np.log10(d["p"].clip(lower=1e-300)),
        significance=lambda d: np.where(d["q"] < Q_THRESHOLD, "q < 0.05", "q >= 0.05"),
    )
    fig = px.scatter(
        f,
        x="outer_minus_inner",
        y="neg_log10_p",
        color="significance",
        color_discrete_map={"q < 0.05": RED, "q >= 0.05": "#999999"},
        category_orders={"significance": ["q >= 0.05", "q < 0.05"]},
        hover_data=["feature"],
        labels={
            "outer_minus_inner": "Mean outer - inner difference (normalized feature units)",
            "neg_log10_p": "-log10 permutation p",
        },
        title=f"{COMPARTMENTS[comp]}: inner vs outer wells, per feature",
    )
    fig.add_vline(x=0, line_dash="dash", line_color="grey")
    fig.update_layout(template="plotly_white", height=460)
    _show(fig, f"{comp}_pp_feat", f"{comp}_plate_position_features")

    st.markdown("**DMSO column check**")
    st.dataframe(
        tests.loc[tests["compartment"] == comp].rename(
            columns={"within_minus_between": "within minus between"}
        ),
        hide_index=True,
        width="stretch",
    )
    dmso = pairs.loc[pairs["compartment"] == comp]
    fig = px.box(
        dmso,
        x="patient_tumor",
        y="correlation",
        color="pair_type",
        points="all",
        labels={"correlation": "Pearson correlation between DMSO well profiles"},
        title=f"{COMPARTMENTS[comp]}: DMSO well pairs by patient",
    )
    fig.update_layout(template="plotly_white", height=460)
    _show(fig, f"{comp}_dmso_patient", f"{comp}_dmso_by_patient")
    fig = px.box(
        dmso,
        x="pair_type",
        y="correlation",
        color="patient_tumor",
        points="all",
        labels={"correlation": "Pearson correlation between DMSO well profiles"},
        title=f"{COMPARTMENTS[comp]}: DMSO well pairs by well pair",
    )
    fig.update_layout(template="plotly_white", height=460)
    _show(fig, f"{comp}_dmso_pair", f"{comp}_dmso_by_well_pair")


def _plate_maps(wells: pd.DataFrame) -> go.Figure:
    """One 8 x 12 plate per patient, outlined DMSO wells (R script 7, facet_wrap)."""
    patients = (
        wells[["patient_tumor", "tumor_type"]]
        .drop_duplicates()
        .assign(
            order=lambda d: d["tumor_type"].map(
                {t: i for i, t in enumerate(TUMOR_TYPES)}
            )
        )
        .sort_values(["order", "patient_tumor"])
    )
    n_cols = 4
    n_rows = -(-len(patients) // n_cols)
    labels = [
        f"{p} ({t})" for p, t in zip(patients["patient_tumor"], patients["tumor_type"])
    ]
    fig = make_subplots(
        rows=n_rows,
        cols=n_cols,
        subplot_titles=labels,
        horizontal_spacing=0.04,
        vertical_spacing=0.08,
    )
    row_letters = list("ABCDEFGH")
    columns = list(range(1, 13))
    for k, patient in enumerate(patients["patient_tumor"]):
        sub = wells.loc[wells["patient_tumor"] == patient]
        grid = (
            sub.pivot_table(
                index="plate_row",
                columns="plate_column",
                values="correlation_to_dmso",
                aggfunc="mean",
            )
            .reindex(index=row_letters, columns=columns)
            .to_numpy(dtype=float)
        )
        r, c = divmod(k, n_cols)
        fig.add_trace(
            go.Heatmap(
                z=grid,
                x=columns,
                y=list(range(1, 9)),
                colorscale="RdBu_r",
                zmid=0,
                zmin=-1,
                zmax=1,
                showscale=k == 0,
                colorbar=dict(title="Pearson r<br>to DMSO"),
                hovertemplate="%{y}%{x}<br>r = %{z:.2f}<extra></extra>",
            ),
            row=r + 1,
            col=c + 1,
        )
        xref = "x" if k == 0 else f"x{k + 1}"
        yref = "y" if k == 0 else f"y{k + 1}"
        for _, well in sub.loc[sub["treatment"] == "DMSO"].iterrows():
            row_index = row_letters.index(well["plate_row"]) + 1
            col_index = int(well["plate_column"])
            fig.add_shape(
                type="rect",
                x0=col_index - 0.5,
                x1=col_index + 0.5,
                y0=row_index - 0.5,
                y1=row_index + 0.5,
                line=dict(color="black", width=2),
                xref=xref,
                yref=yref,
            )
    fig.update_yaxes(autorange="reversed", showticklabels=False)
    fig.update_xaxes(showticklabels=False)
    fig.update_layout(
        title=f"{COMPARTMENTS[wells['compartment'].iloc[0]]}: well correlation to plate DMSO "
        "(black outline = DMSO well)",
        template="plotly_white",
        height=260 * n_rows,
    )
    return fig


def _well_correlation(comp: str, filters: Filters) -> None:
    wells = _differential("well_dmso_correlation")
    if wells is None:
        return
    wells = apply_global_filters(wells.loc[wells["compartment"] == comp], filters)
    if wells.empty:
        st.warning("No wells left after subsetting.")
        return
    wells = _with_treatment_label(wells)
    _show(_plate_maps(wells), f"{comp}_plate_maps", f"{comp}_plate_maps")

    order = _treatment_order(wells)
    fig = px.box(
        wells,
        x="treatment_label",
        y="correlation_to_dmso",
        points=False,
        category_orders={"treatment_label": order},
        title=f"{COMPARTMENTS[comp]}: well correlation to plate DMSO by treatment",
    )
    fig.update_traces(fillcolor="#e5e5e5", line_color="#7f7f7f")
    strip = px.strip(
        wells,
        x="treatment_label",
        y="correlation_to_dmso",
        color="tumor_type",
        category_orders={"treatment_label": order, "tumor_type": TUMOR_TYPES},
        color_discrete_map=TUMOR_TYPE_PALETTE,
    )
    fig.add_traces(strip.data)
    fig.add_hline(y=0, line_dash="dash", line_color="grey")
    fig.update_layout(
        template="plotly_white",
        height=520,
        xaxis_title=None,
        yaxis_title="Pearson r to plate DMSO<br>(DMSO wells: leave-one-out)",
        legend_title="Tumor type",
    )
    _show(fig, f"{comp}_wells_box", f"{comp}_well_correlation_by_treatment")


def _mek(comp: str, filters: Filters) -> None:
    """Count tiles, UpSet and category bars, and recurrence bars (R script 9)."""
    tests = _differential("mek_contrast_level_tests")
    if tests is None:
        return
    lv = tests.loc[tests["compartment"] == comp].assign(
        significant=lambda d: d["q"] < Q_THRESHOLD,
        level_label=lambda d: d["level"].map(LEVEL_LABELS),
        group_label=lambda d: np.select(
            [d["level"] == "within_patient", d["level"] == "within_tumor_type"],
            [d["group"] + " (" + d["tumor_type"].astype(str) + ")", d["group"]],
            default=d["level"].map(POOLED_GROUP_LABELS),
        ),
    )
    lv["tumor_type"] = lv["tumor_type"].astype(str)
    type_rank = {t: i for i, t in enumerate(TUMOR_TYPES)}
    level_rank = {k: i for i, k in enumerate(LEVEL_LABELS)}

    counts = (
        lv.groupby(
            ["level", "level_label", "tumor_type", "group_label", "contrast"],
            observed=True,
        )
        .agg(n_significant=("significant", "sum"))
        .reset_index()
        .assign(
            level_rank=lambda d: d["level"].map(level_rank),
            type_rank=lambda d: d["tumor_type"].map(lambda t: type_rank.get(t, 99)),
        )
        .sort_values(["level_rank", "type_rank", "group_label"])
    )
    max_count = max(float(counts["n_significant"].max()), 1.0)
    nice = [v for v in (0, 1, 5, 10, 25, 50, 100, 250, 500, 1000) if v <= max_count]
    fig = _tile_panels(
        counts,
        panel="level_label",
        row="group_label",
        col="contrast",
        value="n_significant",
        panels=[
            LEVEL_LABELS[k]
            for k in LEVEL_LABELS
            if LEVEL_LABELS[k] in set(counts["level_label"])
        ],
        cols=CONTRASTS,
        title=f"{COMPARTMENTS[comp]}: features with q < {Q_THRESHOLD}",
        colorscale=SEQUENTIAL,
        colorbar_title="Significant<br>features",
        zmin=0,
        label_values=True,
        transform=np.sqrt,
        tickvals=[float(np.sqrt(v)) for v in nice],
        ticktext=[str(v) for v in nice],
    )
    _show(fig, f"{comp}_mek_counts", f"{comp}_mek_significant_counts")

    level_label = st.selectbox(
        "Level for UpSet and category plots",
        list(LEVEL_LABELS.values()),
        key=f"{comp}_mek_level",
    )
    raw_level = {v: k for k, v in LEVEL_LABELS.items()}[level_label]
    sig = lv.loc[(lv["level"] == raw_level) & lv["significant"]]
    if sig.empty:
        st.info("No features reach q < 0.05 at this level.")
    else:
        element = np.where(
            sig["group"] == "all", sig["feature"], sig["group"] + " | " + sig["feature"]
        )
        membership = (
            sig.assign(element=element, hit=True)
            .pivot_table(
                index="element", columns="contrast", values="hit", aggfunc="max"
            )
            .reindex(columns=CONTRASTS)
            .fillna(False)
            .astype(bool)
        )
        combos = _upset_combinations(membership, CONTRASTS)
        if combos.empty:
            st.info("No contrast combinations at this level.")
        else:
            fig = upset_plot(
                _set_sizes(membership, CONTRASTS),
                combos,
                CONTRASTS,
                top_n=15,
                title=f"{COMPARTMENTS[comp]} | {level_label}",
            )
            if fig is not None:
                _show(fig, f"{comp}_mek_upset", f"{comp}_mek_upset_{raw_level}")

        cat = sig.assign(
            feature_category=sig["feature_category"].fillna("Other"),
            row_label=np.where(
                sig["group"] == "all", sig["feature_object"], sig["group_label"]
            ),
            direction=np.where(
                raw_level == "across_tumor_types_differ",
                "Differs by tumor type",
                np.where(sig["mean"] > 0, "Higher in MEKi", "Lower in MEKi"),
            ),
        )
        bars = (
            cat.groupby(
                ["direction", "contrast", "row_label", "feature_category"],
                observed=True,
            )
            .size()
            .rename("n_features")
            .reset_index()
        )
        fig = px.bar(
            bars,
            x="n_features",
            y="row_label",
            color="feature_category",
            facet_col="contrast",
            facet_row="direction",
            orientation="h",
            barmode="stack",
            labels={"n_features": f"Features with q < {Q_THRESHOLD}", "row_label": ""},
            title=f"{COMPARTMENTS[comp]} | {level_label}: feature categories",
        )
        fig.update_layout(
            template="plotly_white",
            height=max(420, 26 * bars["row_label"].nunique() + 160),
        )
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        _show(fig, f"{comp}_mek_cat", f"{comp}_mek_categories_{raw_level}")

    # recurrence: how many groups each feature is significant in, by pooled test
    within = lv.loc[lv["level"].isin(["within_patient", "within_tumor_type"])]
    n_groups = (
        within.groupby(["level", "contrast", "feature"], observed=True)["significant"]
        .sum()
        .rename("n_groups")
        .reset_index()
    )
    pooled_of = {
        "within_patient": "across_patients",
        "within_tumor_type": "across_tumor_types_shared",
    }
    pooled = lv.loc[
        lv["group"] == "all", ["level", "contrast", "feature", "significant"]
    ].rename(columns={"level": "pooled_level", "significant": "pooled_sig"})
    rec = (
        n_groups.assign(pooled_level=n_groups["level"].map(pooled_of))
        .merge(pooled, on=["pooled_level", "contrast", "feature"], how="left")
        .loc[lambda d: d["n_groups"] > 0]
    )
    if rec.empty:
        return
    rec = rec.assign(
        unit=rec["level"].map(
            {"within_patient": "Patients", "within_tumor_type": "Tumor types"}
        ),
        pooled=np.where(
            rec["pooled_sig"].fillna(False).astype(bool), "q < 0.05", "q >= 0.05"
        ),
        n_label=rec["n_groups"].astype(int).astype(str),
    )
    rec_counts = (
        rec.groupby(["unit", "contrast", "n_label", "pooled"])
        .size()
        .rename("n_features")
        .reset_index()
    )
    fig = px.bar(
        rec_counts,
        x="n_label",
        y="n_features",
        color="pooled",
        facet_row="unit",
        facet_col="contrast",
        color_discrete_map={"q >= 0.05": "#999999", "q < 0.05": RED},
        category_orders={"pooled": ["q >= 0.05", "q < 0.05"]},
        labels={
            "n_label": "Number of patients / tumor types in which the feature is significant",
            "n_features": "Features",
            "pooled": "Pooled test (across patients / tumor type means)",
        },
        title=f"{COMPARTMENTS[comp]}: recurrence of significant features",
    )
    fig.update_layout(template="plotly_white", height=560)
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    _show(fig, f"{comp}_mek_recur", f"{comp}_mek_recurrence")


def _compartment_section(comp: str, filters: Filters) -> None:
    summary = _differential("log2fc_summary")
    if summary is None:
        return
    st.markdown("#### Morphology heatmap")
    _morphology_heatmap(summary, comp, filters)
    st.markdown("#### Morphology vs viability")
    _morphology_vs_viability(summary, comp, filters)
    st.markdown("#### Plate position")
    _plate_position(comp, filters)
    st.markdown("#### Well correlation to DMSO")
    _well_correlation(comp, filters)
    st.markdown("#### MEK signatures")
    _mek(comp, filters)


DIFFERENTIAL_SECTIONS = {
    "Viability heatmap": viability_heatmap_section,
    "Organoid": lambda filters: _compartment_section("organoid", filters),
    "Single cell": lambda filters: _compartment_section("single_cell", filters),
}
