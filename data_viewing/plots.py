"""Generic interactive plot explorer: subset, color and facet by any column."""

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from data_io import GLOBAL_FILTER_COLUMNS, _natural_key
from palettes import order_for, palette_for

NONE = "(none)"
MAX_FACETS = 30
CONTINUOUS_MIN_UNIQUE = 12
PNG_DPI = 600
CSS_DPI = 96  # plotly lays figures out at 96 px per inch
PLOT_KINDS = ["scatter", "box", "violin", "histogram", "bar", "heatmap"]


# ---------------------------------------------------------------------------
# Column helpers
# ---------------------------------------------------------------------------
def numeric_columns(df: pd.DataFrame) -> list[str]:
    return [
        c
        for c in df.columns
        if pd.api.types.is_numeric_dtype(df[c])
        and not pd.api.types.is_bool_dtype(df[c])
    ]


def categorical_columns(df: pd.DataFrame, max_unique: int = 200) -> list[str]:
    """Columns usable for color/facet/subset: non-numeric or low-cardinality."""
    cols = []
    for c in df.columns:
        if pd.api.types.is_list_like(df[c].iloc[0]) if len(df) else False:
            continue
        n = df[c].nunique(dropna=True)
        is_numeric = pd.api.types.is_numeric_dtype(
            df[c]
        ) and not pd.api.types.is_bool_dtype(df[c])
        if (not is_numeric and n <= max_unique) or (
            is_numeric and n <= CONTINUOUS_MIN_UNIQUE
        ):
            cols.append(c)
    return cols


def _order_categories(df: pd.DataFrame, col: str) -> None:
    """Sort category levels naturally (so dose 0.1 < 1 < 10) for stable legends/facets."""
    if col not in df.columns or pd.api.types.is_numeric_dtype(df[col]):
        return
    levels = sorted(df[col].dropna().astype(str).unique(), key=_natural_key)
    df[col] = pd.Categorical(df[col].astype(str), categories=levels, ordered=True)


# ---------------------------------------------------------------------------
# Subsetting
# ---------------------------------------------------------------------------
def apply_global_filters(
    df: pd.DataFrame, filters: dict[str, list[str]]
) -> pd.DataFrame:
    """Apply sidebar filters wherever the column exists (silently skipped otherwise)."""
    for col, selected in filters.items():
        if selected and col in df.columns:
            df = df[df[col].astype(str).isin(selected)]
    return df


def local_subset(df: pd.DataFrame, key: str) -> pd.DataFrame:
    """Section-level subset: keep/exclude values of any categorical column."""
    with st.expander("Subset rows by any metadata column", expanded=False):
        cols = [c for c in categorical_columns(df) if c not in GLOBAL_FILTER_COLUMNS]
        cols = _global_columns_in(df) + cols
        chosen = st.multiselect("Columns to subset", cols, key=f"{key}_sub_cols")
        for col in chosen:
            values = sorted(df[col].dropna().astype(str).unique(), key=_natural_key)
            mode = st.radio(
                f"{col}",
                ["keep", "exclude"],
                horizontal=True,
                key=f"{key}_sub_mode_{col}",
            )
            picked = st.multiselect(
                f"{col} values",
                values,
                key=f"{key}_sub_vals_{col}",
                label_visibility="collapsed",
            )
            if picked:
                mask = df[col].astype(str).isin(picked)
                df = df[mask] if mode == "keep" else df[~mask]
    return df


def _global_columns_in(df: pd.DataFrame) -> list[str]:
    return [c for c in GLOBAL_FILTER_COLUMNS if c in df.columns]


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------
def png_download(fig: go.Figure, key: str, filename: str) -> None:
    """Render a 600 dpi PNG on demand (only PNG is offered, per project convention)."""
    if not st.button("Prepare PNG (600 dpi)", key=f"{key}_png_prep"):
        return
    width = int(fig.layout.width or 1000)
    height = int(fig.layout.height or 700)
    try:
        with st.spinner("Rendering PNG..."):
            png = fig.to_image(
                format="png", width=width, height=height, scale=PNG_DPI / CSS_DPI
            )
    except Exception as err:  # kaleido needs a Chrome install
        st.error(
            f"PNG export failed: {err}. Run `plotly_get_chrome` once to install Chrome."
        )
        return
    st.download_button(
        "Download PNG",
        png,
        file_name=f"{filename}.png",
        mime="image/png",
        key=f"{key}_png_dl",
    )


def missing_notice(dataset_label: str, produced_by: str, path) -> None:
    st.info(
        f"No results for **{dataset_label}** yet. Run `{produced_by}` to generate "
        f"`{path}`."
    )


# ---------------------------------------------------------------------------
# Generic explorer
# ---------------------------------------------------------------------------
def explorer(
    df: pd.DataFrame,
    key: str,
    filters: dict[str, list[str]],
    kinds: list[str] | None = None,
    defaults: dict | None = None,
    title: str = "",
    height: int = 650,
) -> go.Figure | None:
    """Interactive plot with plot-type, axis, color, facet and subset controls.

    Parameters
    ----------
    df : data with canonical metadata columns.
    key : unique widget-key prefix for this section.
    filters : shared sidebar filters (column -> selected values).
    kinds : allowed plot types (defaults to all).
    defaults : initial values for ``kind``, ``x``, ``y``, ``color``, ``facet``.
    """
    kinds = kinds or PLOT_KINDS
    defaults = defaults or {}
    df = apply_global_filters(df, filters)
    df = local_subset(df, key)
    if df.empty:
        st.warning("No rows left after subsetting.")
        return None
    st.caption(f"{len(df):,} rows after subsetting")

    numeric = numeric_columns(df)
    categorical = categorical_columns(df)
    all_cols = [c for c in df.columns if not pd.api.types.is_list_like(df[c].iloc[0])]

    def pick(label, options, default, col, allow_none=True):
        opts = ([NONE] if allow_none else []) + list(options)
        idx = opts.index(default) if default in opts else 0
        return col.selectbox(label, opts, index=idx, key=f"{key}_{label}")

    c1, c2, c3 = st.columns(3)
    kind = pick(
        "Plot type", kinds, defaults.get("kind", kinds[0]), c1, allow_none=False
    )
    x = pick(
        "X",
        all_cols,
        defaults.get("x"),
        c2,
        allow_none=kind in ("histogram", "box", "violin", "bar"),
    )
    y_options = numeric if kind in ("box", "violin", "bar") else all_cols
    y = pick("Y", y_options, defaults.get("y"), c3, allow_none=kind == "histogram")

    c4, c5, c6 = st.columns(3)
    color = pick("Color", all_cols, defaults.get("color"), c4)
    facet = pick("Facet", categorical, defaults.get("facet"), c5)
    facet_wrap = c6.slider(
        "Facets per row", 1, 8, 3, key=f"{key}_wrap", disabled=facet == NONE
    )

    with st.expander("Appearance", expanded=False):
        a1, a2, a3, a4 = st.columns(4)
        opacity = a1.slider("Opacity", 0.05, 1.0, 0.6, key=f"{key}_opacity")
        size = a2.slider("Point size", 1, 15, 4, key=f"{key}_size")
        max_points = a3.number_input(
            "Max scatter points", 1000, 500_000, 50_000, step=5000, key=f"{key}_maxpts"
        )
        agg = a4.selectbox(
            "Aggregation (bar)", ["mean", "median", "sum", "count"], key=f"{key}_agg"
        )
        b1, b2, b3 = st.columns(3)
        as_cat = b1.checkbox(
            "Treat color as categorical", value=True, key=f"{key}_ascat"
        )
        log_x = b2.checkbox("Log X", key=f"{key}_logx")
        log_y = b3.checkbox("Log Y", key=f"{key}_logy")
        nbins = st.slider("Histogram bins", 5, 200, 50, key=f"{key}_bins")
        background = st.checkbox(
            "Scatter: show all points in every facet (grey), color only that facet's points",
            key=f"{key}_background",
            disabled=facet == NONE,
        )

    plot_df = df.copy()
    x = None if x == NONE else x
    y = None if y == NONE else y
    color = None if color == NONE else color
    facet = None if facet == NONE else facet

    if facet:
        n_facets = plot_df[facet].nunique()
        if n_facets > MAX_FACETS:
            top = plot_df[facet].value_counts().head(MAX_FACETS).index
            st.warning(
                f"{facet} has {n_facets} levels; showing the {MAX_FACETS} largest."
            )
            plot_df = plot_df[plot_df[facet].isin(top)]
    for col in (x, color, facet):
        if col and (col != color or as_cat or col in ("dose",)):
            _order_categories(plot_df, col)
    if color and as_cat and pd.api.types.is_numeric_dtype(plot_df[color]):
        plot_df[color] = plot_df[color].astype(str)

    common = dict(
        color=color, facet_col=facet, facet_col_wrap=facet_wrap if facet else 0
    )
    common.update(_palette_kwargs(plot_df, x, color, facet))
    if not facet:
        common.pop("facet_col_wrap")
    log = dict(log_x=log_x, log_y=log_y)

    try:
        fig = _build(
            kind,
            plot_df,
            x,
            y,
            common,
            log,
            opacity,
            size,
            max_points,
            agg,
            nbins,
            background and facet is not None,
        )
    except Exception as err:
        st.error(f"Could not draw this combination: {err}")
        return None

    fig.update_layout(title=title, height=height, template="plotly_white")
    if facet:
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    st.plotly_chart(fig, width="stretch", key=f"{key}_chart")
    png_download(fig, key, key)
    return fig


def _palette_kwargs(df: pd.DataFrame, x, color, facet) -> dict:
    """R-theme colors for the color column and R-defined orders for x/color/facet."""
    kwargs: dict = {}
    if color and not (
        pd.api.types.is_numeric_dtype(df[color])
        and df[color].nunique() > CONTINUOUS_MIN_UNIQUE
    ):
        levels = [str(v) for v in df[color].dropna().unique()]
        colors = palette_for(color, levels)
        if colors:
            kwargs["color_discrete_map"] = colors
    orders = {}
    for col in (x, color, facet):
        if col and not pd.api.types.is_numeric_dtype(df[col]):
            levels = [str(v) for v in df[col].dropna().unique()]
            order = order_for(col, levels)
            if order:
                orders[col] = order
    if orders:
        kwargs["category_orders"] = orders
    return kwargs


def _build(
    kind, df, x, y, common, log, opacity, size, max_points, agg, nbins, background=False
):
    if kind == "scatter":
        if x is None or y is None:
            raise ValueError("choose both X and Y")
        if len(df) > max_points:
            st.caption(f"Showing a random {max_points:,} of {len(df):,} points.")
            df = df.sample(int(max_points), random_state=0)
        fig = px.scatter(
            df, x=x, y=y, opacity=opacity, render_mode="webgl", **common, **log
        )
        fig.update_traces(marker=dict(size=size))
        if background:
            _add_grey_background(fig, df, x, y, size)
        return fig
    if kind in ("box", "violin"):
        if y is None:
            raise ValueError("choose a numeric Y")
        fn = px.box if kind == "box" else px.violin
        kwargs = dict(box=True) if kind == "violin" else {}
        return fn(df, x=x, y=y, **kwargs, **common, log_y=log["log_y"])
    if kind == "histogram":
        col = x or y
        if col is None:
            raise ValueError("choose X")
        return px.histogram(
            df,
            x=col,
            nbins=nbins,
            barmode="overlay",
            opacity=0.65,
            **common,
            log_x=log["log_x"],
        )
    if kind == "bar":
        if y is None:
            raise ValueError("choose a numeric Y")
        group = [c for c in (x, common["color"], common["facet_col"]) if c]
        if not group:
            raise ValueError("choose X, color or facet to group by")
        grouped = (
            df.groupby(group, observed=True)[y].agg(agg).reset_index()
            if agg != "count"
            else df.groupby(group, observed=True)[y].count().reset_index()
        )
        return px.bar(grouped, x=x, y=y, barmode="group", **common, log_y=log["log_y"])
    if kind == "heatmap":
        if x is None or y is None:
            raise ValueError("choose X and Y")
        return px.density_heatmap(
            df,
            x=x,
            y=y,
            z=None,
            facet_col=common["facet_col"],
            facet_col_wrap=common.get("facet_col_wrap", 0) or None,
            nbinsx=nbins,
            nbinsy=nbins,
            color_continuous_scale="Viridis",
        )
    raise ValueError(f"unknown plot type {kind}")


def _add_grey_background(
    fig: go.Figure, df: pd.DataFrame, x: str, y: str, size: int
) -> None:
    """Draw every point in grey behind the colored points of each facet panel."""
    panels = {(t.xaxis, t.yaxis) for t in fig.data}
    grey = [
        go.Scattergl(
            x=df[x],
            y=df[y],
            mode="markers",
            xaxis=xaxis,
            yaxis=yaxis,
            marker=dict(size=size, color="lightgrey", opacity=0.25),
            showlegend=False,
            hoverinfo="skip",
        )
        for xaxis, yaxis in sorted(panels, key=str)
    ]
    n_colored = len(fig.data)
    fig.add_traces(grey)
    # background traces first so the colored points draw on top
    n_grey = len(grey)
    fig.data = fig.data[n_colored:] + fig.data[:n_colored]
    assert len(fig.data) == n_colored + n_grey
