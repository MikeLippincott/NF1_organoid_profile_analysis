# data_viewing

Interactive Streamlit app for re-drawing the repo's plots on the fly.

```bash
source uv_setup.sh          # once; installs streamlit, plotly, kaleido
bash data_viewing/run_app.sh [port]   # default port 8501
# or: streamlit run data_viewing/app.py
```

## Layout
- One tab per module: `1.EDA`, `3.viability_prediction_models`, `4.linear_modeling`.
- One section per analysis type (pick it with the radio buttons; only the selected
  section loads its data).
- Every plot lets you choose plot type, X, Y, **color** and **facet** by any
  column, and **subset** rows by any metadata column (keep or exclude values).
- The sidebar holds shared subset filters (patient, treatment, dose, class, target,
  therapeutic category, well, image mode, modality) that apply to every plot whose
  table has that column.
- Metadata column names are harmonized (`Metadata_Experiment_Treatment` and
  `Metadata_treatment` both become `treatment`; see `data_io.CANONICAL_COLUMNS`).
- "Prepare PNG (600 dpi)" saves the current figure as a PNG. The first export may need
  Chrome: run `plotly_get_chrome` once.

## Data
The app only reads precomputed parquet results (nothing is recomputed):

| Module | Read from |
|---|---|
| 1.EDA | `1.EDA/results/**` |
| 3.viability_prediction_models | `3.viability_prediction_models/model_results/combined_*.parquet` |
| 4.linear_modeling | `4.linear_modeling/results/linear_modeling/*.parquet` |

A section whose results do not exist yet shows which script to run.

## Files
- `run_app.sh`: launcher (finds the repo root and `.venv`).
- `app.py`: page, sidebar filters, tabs.
- `sections.py`: one function per analysis type.
- `plots.py`: generic plot explorer (subset / color / facet controls, PNG export).
- `data_io.py`: paths, dataset registry, metadata harmonization.
