"""NF1 organoid profile viewer.

Run from anywhere inside the repo with:
    streamlit run data_viewing/app.py
"""

import plotly.io as pio
import streamlit as st
from data_io import GLOBAL_FILTER_COLUMNS, all_filterable_paths, global_filter_options
from sections import EDA_SECTIONS, LINEAR_MODELING_SECTIONS, VIABILITY_SECTIONS

pio.templates.default = "plotly_white"
st.set_page_config(page_title="NF1 organoid data viewer", layout="wide")
st.title("NF1 organoid profile data viewer")

# ---------------------------------------------------------------------------
# Shared sidebar filters (applied wherever a table has the column)
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

# ---------------------------------------------------------------------------
# One tab per module, one selectable section per analysis type
# ---------------------------------------------------------------------------
MODULES = {
    "1.EDA": EDA_SECTIONS,
    "3.viability_prediction_models": VIABILITY_SECTIONS,
    "4.linear_modeling": LINEAR_MODELING_SECTIONS,
}

tabs = st.tabs(list(MODULES))
for tab, (module, sections) in zip(tabs, MODULES.items()):
    with tab:
        # only the selected section renders, so heavy tables load on demand
        section = st.radio(
            "Analysis", list(sections), horizontal=True, key=f"section_{module}"
        )
        st.subheader(section)
        sections[section](filters)
