# NF1_organoid_profile_analysis
This repo contains analysis code of profiles generated in multiple image-based profiling repos.

- 3D image-based profiles of NF1 organoids were generated from: [NF1 3D Organoid profiling pipeline](https://github.com/WayScience/NF1_3D_organoid_profiling_pipeline)
- 2D image-based profiles of NF1 organoids were generated from: [NF1 2D organoid profiling pipeline](https://github.com/WayScience/NF1_2D_organoid_profiling_pipeline)
- Cell types will be predicted using: [Cell type prediction in NF1 organoids](https://github.com/WayScience/NF1_organoid_cell_segmentation)

## Repo modules
- `0.download_data`: Download image-based profiles from the internet (not yet available)
- `1.EDA`: Exploratory data analysis of image-based profiles
- `2.2d_vs_3d_analysis`: Comparison of 2D and 3D profiles (patient correlation, mAP, drug hits, entropy, kBET, sparse CCA)
- `3.viability_prediction_models`: Viability prediction models trained on the profiles
- `4.linear_modeling`: Per-patient linear models of every feature on treatment plus count and technical covariates, and the variance, variate-importance and hit analyses built on them (see [`4.linear_modeling/README.md`](4.linear_modeling/README.md))

## Computational environment
Notebooks in this repo are split between Python and R, each managed by a separate environment.

### Python (uv)
Python notebooks use a `uv`-managed virtual environment defined in `pyproject.toml`/`uv.lock`.

```bash
source uv_setup.sh
```

This creates `.venv` and registers a `python3` Jupyter kernel.
Select this kernel when running the Python notebooks.

### R (uvr)
R notebooks use a `uvr`-managed R library defined in `uvr.toml`/`uvr.lock` (R `>=4.3.0`, set by `r_version`).
`uv_setup.sh` also sets this up: it runs `uvr sync` to install the packages into `.uvr/library`, checks the uvr-managed R against `r_version`, and registers an IRkernel Jupyter kernel named after the `uvr.toml` project (`NF1_organoid_profile_analysis`) that runs that R with the project library.
Select this kernel when running the R notebooks.

R scripts run in the same environment with `uvr run`:

```bash
uvr run 4.linear_modeling/scripts/6.plot_variate_importance.r
```

Add a package with `uvr add <package>` so it is recorded in `uvr.toml` and `uvr.lock`.

### R (mamba, older modules)
The R notebooks in `1.EDA` and `2.2d_vs_3d_analysis` were run with the mamba environment in `environments/r_env.yml` (`mamba env create -f environments/r_env.yml`) and a generic `ir` kernel.
