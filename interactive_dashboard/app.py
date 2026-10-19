"""NF1 organoid profile viewer.

Run from anywhere inside the repo with:
    streamlit run interactive_dashboard/app.py
"""

import plotly.io as pio
import streamlit as st
from background import module_background, section_background
from data_io import (
    BUCKET_SYNC_ERROR,
    GLOBAL_FILTER_COLUMNS,
    all_filterable_paths,
    global_filter_options,
)
from memory_trace import sidebar_memory, trace_memory
from sections import (  # noqa: F401 -- VIABILITY_SECTIONS kept for the commented-out tab below
    DIFFERENTIAL_QC_SECTIONS,
    DIFFERENTIAL_SECTIONS,
    EDA_SECTIONS,
    LINEAR_MODELING_SECTIONS,
    OVERVIEW_SECTIONS,
    VIABILITY_SECTIONS,
)

pio.templates.default = "plotly_white"
st.set_page_config(page_title="NF1 organoid profile data viewer", layout="wide")

if BUCKET_SYNC_ERROR:
    st.error(
        f"**Couldn't load data from the HF Bucket** -- every tab below will show "
        f'"no results" until this is fixed.\n\n{BUCKET_SYNC_ERROR}',
        icon="🚨",
    )

# ---------------------------------------------------------------------------
# Modules, grouped into two pages: Biological (treatment-effect findings) and
# Technical (statistical/ML methodology). Each module keeps its own tab,
# section radio, and background text -- only which page it renders under changes.
# ---------------------------------------------------------------------------
BIOLOGICAL_MODULES = {
    "Overview": OVERVIEW_SECTIONS,
    "EDA": EDA_SECTIONS,
    "Linear modeling": LINEAR_MODELING_SECTIONS,
    "Treatment vs DMSO": DIFFERENTIAL_SECTIONS,
}
TECHNICAL_MODULES = {
    # "Viability prediction" is hidden for now (its results need another
    # pass); the module import and its sections.py code are kept intact -- restore
    # this tab by uncommenting the line below.
    # "Viability prediction": VIABILITY_SECTIONS,
    "QC checks": DIFFERENTIAL_QC_SECTIONS,
}


def _render_modules(modules: dict, filters: dict[str, list[str]]) -> None:
    """One tab per module, one selectable section per analysis type."""
    tabs = st.tabs(list(modules))
    for tab, (module, sections) in zip(tabs, modules.items()):
        with tab:
            module_background(module)
            # only the selected section renders, so heavy tables load on demand
            section = st.radio(
                "Analysis", list(sections), horizontal=True, key=f"section_{module}"
            )
            st.subheader(section)
            section_background(section)
            with trace_memory(f"{module} / {section}"):
                sections[section](filters)


def biological_page() -> None:
    st.title("Biological")
    _render_modules(BIOLOGICAL_MODULES, filters)


def technical_page() -> None:
    st.title("Technical")
    _render_modules(TECHNICAL_MODULES, filters)


pg = st.navigation(
    [
        st.Page(biological_page, title="Biological", default=True),
        st.Page(technical_page, title="Technical"),
    ]
)

# ---------------------------------------------------------------------------
# Shared sidebar filters (applied wherever a table has the column, on either page)
# ---------------------------------------------------------------------------
st.sidebar.header("Subset (all tabs)")
st.sidebar.caption(
    "Filters apply to every plot whose table has that column. "
    "Leave empty to keep everything."
)
options = global_filter_options(all_filterable_paths())
filters: dict[str, list[str]] = {}
for column in GLOBAL_FILTER_COLUMNS:
    if column in options:
        filters[column] = st.sidebar.multiselect(
            column, options[column], key=f"global_{column}"
        )
if st.sidebar.button("Clear all filters"):
    for column in GLOBAL_FILTER_COLUMNS:
        st.session_state.pop(f"global_{column}", None)
    st.rerun()

pg.run()

# shown last so the readout reflects memory after the section has rendered
sidebar_memory()
