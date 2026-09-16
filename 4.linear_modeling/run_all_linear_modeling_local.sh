#!/bin/bash

# set -euo pipefail

git_root=$(git rev-parse --show-toplevel)
if [ -z "$git_root" ]; then
    echo "Error: Could not find the git root directory."
    exit 1
fi

module_dir="$git_root/4.linear_modeling"
venv_python="$git_root/.venv/bin/python"

if [[ ! -x "$venv_python" ]]; then
    echo "Error: uv-managed virtual environment not found at $git_root/.venv (run uv_setup.sh first)."
    exit 1
fi

# regenerate scripts/ from notebooks/ so the mirrored .py/.r files never drift
uv run jupyter nbconvert --to=script \
    --FilesWriter.build_directory="$module_dir/scripts/" \
    "$module_dir"/notebooks/*.ipynb

# deactivate any existing conda environment
conda deactivate >/dev/null 2>&1 || true
# deactivate any existing venv environment
deactivate >/dev/null 2>&1 || true

# 0-2: fit the linear models (pooled, and per-patient)
# uv run python "$module_dir/script/0.linear_modeling.py"
# uv run python "$module_dir/scripts/2.variate_importance.py"
# uv run python "$module_dir/scripts/3.per_patient_variate_importance.py"

# 3-5: derive analysis results from the fits above (no new model fitting)
uv run python "$module_dir/scripts/4.find_significant_features.py"
uv run python "$module_dir/scripts/5.explore_leads.py"
uv run python "$module_dir/scripts/6.significant_variates_per_patient_feature.py"

# 19: nested variance decomposition (patient/treatment/well/fov/residual) —
# operates on the raw profiles directly, cross-references 1.variate_importance.py's
# term_summary.parquet, so it just needs to run after step 1 above
uv run python "$module_dir/scripts/20.variance_decomposition_by_replicate_level.py"

# 6-18: plotting only, one idea per script — each reads already-saved results
conda run -n gff_figure_env Rscript "$module_dir/scripts/7.plot_linear_modeling_results.r"
uv run python "$module_dir/scripts/8.plot_leads_breadth_vs_magnitude.py"
uv run python "$module_dir/scripts/9.plot_leads_moa_hit_rate.py"
uv run python "$module_dir/scripts/10.plot_leads_channel_feature_heatmap.py"
uv run python "$module_dir/scripts/11.plot_variate_importance_term_summary.py"
uv run python "$module_dir/scripts/12.plot_variate_importance_feature_heatmap.py"
uv run python "$module_dir/scripts/13.plot_variate_importance_organoid_vs_sc.py"
uv run python "$module_dir/scripts/14.plot_per_patient_term_summary.py"
uv run python "$module_dir/scripts/15.plot_per_patient_treatment_importance_by_patient.py"
uv run python "$module_dir/scripts/16.plot_per_patient_top_variate_by_patient.py"
uv run python "$module_dir/scripts/17.plot_per_patient_coefficient_overview.py"
uv run python "$module_dir/scripts/18.plot_per_patient_volcano.py"
uv run python "$module_dir/scripts/19.plot_per_patient_treatment_heatmap.py"

echo "4.linear_modeling pipeline complete."
