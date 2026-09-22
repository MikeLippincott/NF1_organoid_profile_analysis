# 4. Linear modeling

This module fits one ordinary least squares (OLS) model per **(patient, treatment/dose, feature)** on the 3D NF1 organoid profiles and uses the saved fit statistics to ask which morphology changes come from treatment rather than from counts or technical (plate / position) covariates.
Each model compares one treatment with DMSO within a single patient:

```text
feature ~ treatment + cell_count + organoid_count + cell_per_organoid_count                      (2.linear_modeling)
        + manhattan_distance_from_center + cell_x_position + cell_y_position
        + cell_z_position + cell_z_depth                                                          (3.linear_modeling_technical_vars)
```

For every model and term the results store the coefficient, the p-value, the Benjamini-Hochberg FDR (`pvalue_fdr`, corrected separately per term across all features, treatments and patients), R2 / adjusted R2, and the variance shares (`term_pct_of_total_var` from type II sums of squares, `residual_pct`).
The steps after 3 only read these saved statistics; no model is refit.

## Pipeline

Every step exists as a notebook (`notebooks/`) and a mirrored script (`scripts/`) with the same base name.
Python steps run in the `uv` environment and R (ggplot2) steps run in the `uvr` environment (see the [root README](../README.md#computational-environment)).

| Step | Language | What it does | Main outputs |
|------|----------|--------------|--------------|
| `0.preprocessing_non_fs_agg` | Python | Aggregates the normalized organoid and single-cell profiles to one row per well | `data/*_norm_aggregated_profile.parquet` |
| `1.calculate_well_manhattan_distance` | Python | Manhattan distance of every well from the plate center | `results/well_manhattan_distance/` |
| `2.linear_modeling` | Python | Treatment + count covariate models | `results/linear_modeling/<profile>.parquet` |
| `3.linear_modeling_technical_vars` | Python | Adds the plate-position and cell-position covariates | `results/linear_modeling/<profile>_technical_model.parquet` |
| `4.variance_decomposition` | Python | Where the variance goes, and why so much of it is residual | `results/decomposed_variance/`, `results/residual_diagnostics/`, `results/variate_vs_residual/`, `results/top_models/` |
| `5.calculate_variate_importance` | Python | Venn / UpSet tables of which features are hits for which term | `results/variate_importance/` |
| `6.plot_variate_importance` | R | Venn / UpSet plots, treatment-only co-occurrence heatmaps, title-free headline panels | `figures/variate_importance/` |
| `7.calculate_variate_class_upsets_and_clustermap` | Python | Variate-class membership, patient UpSets, patient x treatment clustermaps, treatment-only features | `results/variate_class_plots/<tag>/` |
| `8.plot_variate_class_upsets_and_clustermap` | R | Plots for step 7 | `figures/variate_class_plots/<tag>/` |
| `9.explore_linear_model_haystacks` | Python | Exhaustive exploration of the technical models: effect sizes, hit landscape, reproducibility, dose response, the "needles" | `results/explore_linear_models/` |
| `10.plot_explore_linear_model_haystacks` | R | Plots for step 9 | `figures/explore_linear_models/explore_linear_models.pdf` |
| `11.headline_results_figure` | R | Multi-panel summary figure from the outputs of 6 and 9 | `figures/headline_results_figure/headline_results_figure.png` |

Step 7 is configurable from the command line (`python scripts/7.calculate_variate_class_upsets_and_clustermap.py --help`), and its run tag (for example `sc_T-O-C-N-M-X-Y-Z-D`) is the first argument of step 8.
Step 8 reads the top treatment-only features from step 7's tables, so it always matches step 7's `--top-features`.

## Hit definitions

Each step uses the hit definition that suits its question, so the counts are not interchangeable:

* `5.calculate_variate_importance`: `pvalue_fdr < 0.05` and `coefficient > 0.1` (increases only).
* `7.calculate_variate_class_upsets_and_clustermap`: `pvalue_fdr < --fdr-max` and standardized effect `> --effect-min` (set on the command line), with count coefficients standardized by the covariate SD.
* `9.explore_linear_model_haystacks`: `pvalue_fdr < 0.05`, `R2 > 0.5`, `adj. R2 > 0` and `|coefficient| > 0.01`.

## Running

```bash
source uv_setup.sh                                   # once: uv + uvr environments and Jupyter kernels
bash 4.linear_modeling/run_all_linear_modeling_local.sh
```

`run_all_linear_modeling_local.sh` regenerates `scripts/` from `notebooks/` with `jupyter nbconvert`, then runs every step in order.
