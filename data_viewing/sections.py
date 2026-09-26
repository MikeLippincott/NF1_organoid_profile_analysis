"""One render function per analysis type, grouped into one function per module."""

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
from data_io import EDA_RESULTS, Dataset, load_dataset, parquet_columns, registry
from plots import (
    NONE,
    apply_global_filters,
    explorer,
    local_subset,
    missing_notice,
    png_download,
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
    label = st.selectbox("Dataset", list(datasets), key=f"{key}_dataset")
    ds = datasets[label]
    df = load_dataset(str(ds.path), columns)
    if prepare:
        df = prepare(df, label)
    defaults = defaults_for(label, df) if defaults_for else {}
    explorer(df, f"{key}_{label}", filters, kinds=kinds, defaults=defaults, title=label)


def _first(df: pd.DataFrame, *names: str) -> str | None:
    return next((n for n in names if n in df.columns), None)


# ---------------------------------------------------------------------------
# 1.EDA
# ---------------------------------------------------------------------------
def umap_section(filters: Filters) -> None:
    def defaults(label, df):
        return {
            "kind": "scatter",
            "x": "UMAP1",
            "y": "UMAP2",
            "color": _first(df, "treatment", "patient_tumor"),
            "facet": None,
        }

    dataset_section(
        registry()["umap"],
        "umap",
        filters,
        "1.EDA/scripts/0.generate_umap.py",
        EDA_RESULTS / "umap",
        defaults,
        kinds=["scatter", "histogram", "heatmap"],
    )


def pca_section(filters: Filters) -> None:
    datasets = registry()["pca"]
    if not datasets:
        missing_notice("PCA", "1.EDA/scripts/2.generate_pca.py", EDA_RESULTS / "pca")
        return
    label = st.selectbox("Dataset", list(datasets), key="pca_dataset")
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
    fig = explorer(
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
    if mode.startswith("Replicate"):
        _pairs_heatmap(pair_files, filters)
    else:
        _per_patient_heatmap(matrix_files)


def _pairs_heatmap(pair_files, filters: Filters) -> None:
    if not pair_files:
        st.info("No `*_correlation_pairs.parquet` files found.")
        return
    pairs_path = st.selectbox(
        "Pairs file", pair_files, format_func=lambda p: p.stem, key="corr_pairs_file"
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
        selection[name] = col.selectbox(
            name, sorted(pairs[name].unique()), key=f"corr_{name}"
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
    fig = px.imshow(
        sub,
        x=[f"{i}: {v}" for i, v in enumerate(labels)],
        y=[f"{i}: {v}" for i, v in enumerate(labels)],
        color_continuous_scale="RdBu_r",
        zmin=-1,
        zmax=1,
        aspect="auto",
        title=" | ".join(f"{k}={v}" for k, v in selection.items()),
    )
    fig.update_layout(height=750, template="plotly_white")
    fig.update_xaxes(showticklabels=len(idx) <= 60)
    fig.update_yaxes(showticklabels=len(idx) <= 60)
    st.plotly_chart(fig, width="stretch", key="corr_pairs_chart")
    png_download(fig, "corr_pairs", "correlation_heatmap")


def _per_patient_heatmap(matrix_files) -> None:
    if not matrix_files:
        st.info("No per-patient correlation matrix file found.")
        return
    df = pd.read_parquet(matrix_files[0])
    c1, c2 = st.columns(2)
    variant = c1.selectbox(
        "Normalization variant", sorted(df["variant"].unique()), key="corr_variant"
    )
    patient = c2.selectbox(
        "Patient",
        sorted(df.loc[df["variant"] == variant, "patient"].unique()),
        key="corr_patient",
    )
    row = df[(df["variant"] == variant) & (df["patient"] == patient)].iloc[0]
    n = int(row["n_samples"])
    matrix = np.asarray(row["correlation"]).reshape(n, n)
    treatments = np.asarray(row["treatment"])
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
    fig = px.imshow(
        matrix[np.ix_(idx, idx)],
        color_continuous_scale="RdBu_r",
        zmin=-1,
        zmax=1,
        aspect="auto",
        title=f"{patient} - {variant} ({len(idx)} cells, sorted by treatment)",
    )
    fig.update_layout(height=750, template="plotly_white")
    fig.update_xaxes(showticklabels=False)
    fig.update_yaxes(showticklabels=False)
    st.plotly_chart(fig, width="stretch", key="corr_pp_chart")
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
            "color": _first(df, "dose"),
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


EDA_SECTIONS = {
    "UMAP": umap_section,
    "PCA": pca_section,
    "Correlation heatmaps": correlation_section,
    "Cell counts": cell_counts_section,
    "Area & volume": area_volume_section,
    "Neighbors": neighbors_section,
    "Intensity": intensity_section,
    "Count vs viability": count_viability_section,
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
def _linear_model_data(key: str, extra_defaults=None):
    datasets = registry()["linear_modeling"]
    if not datasets:
        missing_notice(
            "linear-model results",
            "4.linear_modeling/scripts/0.linear_modeling.py",
            "4.linear_modeling/results/linear_modeling",
        )
        return None, None
    label = st.selectbox("Model results", list(datasets), key=f"{key}_dataset")
    df = load_dataset(str(datasets[label].path))
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


LINEAR_MODELING_SECTIONS = {
    "Volcano (effects vs significance)": lm_volcano_section,
    "Effect sizes": lm_effects_section,
    "Model fit": lm_fit_section,
}
