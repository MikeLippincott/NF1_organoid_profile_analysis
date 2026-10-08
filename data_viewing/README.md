---
title: NF1 organoid profile data viewer
emoji: 📉
colorFrom: yellow
colorTo: red
sdk: streamlit
app_file: streamlit_app.py
pinned: false
---

# NF1 organoid profile data viewer

Interactive Streamlit app for exploring the NF1 organoid profiling analysis: drug/plate layout, EDA (UMAP, PCA, correlation, cell counts, area & volume, neighbors, intensity), linear-modeling results (volcano, effect sizes, model fit, UpSet, variate significance), and differential analysis against DMSO (morphology and viability, viability by tumor type, plate-position checks, MEK signatures).

## Running the app

Run everything from this `data_viewing/` directory with [`just`](https://github.com/casey/just):

```bash
just --list                # list the recipes
just run_app_locally       # start the app (http://localhost:8501)
```

Bare `just` only lists the recipes. It never uploads.

Without `just`, the same thing is:

```bash
pip install -r requirements.txt
bash run_app.sh [port]     # default port 8501
# or: streamlit run app.py
```

`streamlit_app.py` is an alias entry point for hosts that look for that exact filename (e.g. Streamlit Community Cloud); it re-execs `app.py`, which holds the real logic.

## Justfile recipes

| Recipe | What it does |
|---|---|
| `just` / `just --list` | List the recipes. Runs nothing else. |
| `just sync_data` | Syncs with the HF Bucket, choosing by system: on a local machine it **uploads** (stage results from the main analysis repo into `data/`, then push to the bucket; needs `HF_TOKEN` with write access); on Streamlit Community Cloud (detected by `/mount/src`) it only **downloads**. |
| `just upload_data_to_hf_bucket` | Upload, whatever the system. |
| `just download_data_from_hf_bucket` | Download the bucket into `data/`, whatever the system. |
| `just run_app_locally` | Run `run_app.sh` (Streamlit on port 8501). Blocks until you stop it. |
| `just all_data` | Same as `sync_data`. |
| `just all` | `all_data` then `run_app_locally`. |

Notes:
- Upload reads whatever is in the main repo's `*/results/` folders. A module with no results yet is skipped with a `skip (missing)` line, so the app shows "No results" for it.
- The system check is a heuristic. If Streamlit moves its mount point, update `on_streamlit` in the `justfile`.
- `all` does not return until the app is stopped (Ctrl+C).
- Recipes call `uv run python3 sync_data.py ...`, so run `./uv_setup.sh` in the repo root once first so `huggingface_hub` is available in `.venv`.

## Layout
- One tab per module: `Overview`, `EDA`, `Linear modeling`, `Treatment vs DMSO`. (`Viability prediction` is implemented in `sections.py` but hidden from the tab bar for now -- see the comment in `app.py`.)
- One section per analysis type (pick it with the radio buttons; only the selected section loads its data). `Treatment vs DMSO` has three: `Viability heatmap` (shared), and `Organoid` and `Single cell`, each with its own morphology heatmap, morphology vs viability scatters, plate-position, well-correlation and MEK plots, using the same plot types as the R figures in `5.differential_analysis/scripts/`.
- Every plot lets you choose plot type, X, Y and **color** by any column, and **facet** by any categorical column. **Shape** (scatter only) is limited to `dose` and `tumor_type`. You can also **subset** rows by any metadata column (keep or exclude values).
- Scatter plots that are faceted have "show all points in every facet (grey)" checked by default, under **Appearance**. Untick it to show only each facet's own points.
- The sidebar holds shared subset filters (patient, treatment, dose, class, target, therapeutic category, MoA, tumor type, well, image mode, modality) that apply to every plot whose table has that column.
- Metadata column names are harmonized (`Metadata_Experiment_Treatment` and `Metadata_treatment` both become `treatment`; see `data_io.CANONICAL_COLUMNS`).
- "Prepare PNG (600 dpi)" saves the current figure as a PNG.

## Data
This repo ships **code only** -- the precomputed results live in an HF Bucket, not in git (`data/` is `.gitignore`'d on purpose). `sync_data.py` moves data between the main analysis repo, `data/`, and that bucket:

```bash
just upload_data_to_hf_bucket      # main repo -> data/ -> bucket  (= uv run python3 sync_data.py upload)
just download_data_from_hf_bucket  # bucket -> data/               (= uv run python3 sync_data.py download)
```

`download` is also what the app does automatically the first time it starts with no local `data/` -- `data_io.py` calls the same `sync_bucket()` into `DATA_DIR` (default: `data/`), cached for the life of the process so it only happens once. Set `HF_BUCKET` (`namespace/bucket-name`) to point at a different bucket; a private one also needs `HF_TOKEN`. On an HF Space that mounts the bucket as a storage volume instead (Space Settings -> Storage Buckets), point `DATA_DIR` at that mount path and the download is skipped entirely, since the files are already there.

| Module | Read from |
|---|---|
| Overview | `data/platemaps/*` |
| EDA | `data/eda/{umap,pca,correlation,cell_counts,neighbors,intensity}` |
| Viability prediction | `data/viability_models/combined_*.parquet` |
| Linear modeling | `data/linear_modeling/{models,variate_importance}` |
| Treatment vs DMSO | `data/differential_analysis/*.parquet` (all tables from `5.differential_analysis/results/`) |

A section whose results aren't present shows which script, in the main repo, produces them -- it never crashes the app.

## Files
- `justfile`: the recipes above (run `just --list`).
- `run_app.sh`: local launcher (called by `run_app_locally`).
- `requirements.txt`: runtime deps.
- `sync_data.py`: moves the trimmed data above between the main repo, `data/`, and the HF Bucket (both directions); `MANIFEST` lists what is copied.
- `streamlit_app.py`: alias entry point; re-execs `app.py`.
- `app.py`: page, sidebar filters, tabs.
- `sections.py`: one function per analysis type.
- `plots.py`: generic plot explorer (subset / color / facet / shape controls, heatmaps, PNG export).
- `data_io.py`: paths, dataset registry, metadata harmonization.
- `palettes.py`: colors/orders ported from the R plotting theme, and `humanize_label()` for dropdown display names.
- `background.py`: the "About" text shown in each tab/section.
