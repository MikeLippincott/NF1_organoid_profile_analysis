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

if ! command -v uvr >/dev/null 2>&1; then
    echo "Error: uvr not found (needed to run the R plotting notebooks)."
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

uv run python "$module_dir/scripts/calculate_well_manhattan_distance.py"
uv run python "$module_dir/scripts/00.preprocessing_non_fs_agg.py"
uv run python "$module_dir/scripts/0.linear_modeling.py"
uv run python "$module_dir/scripts/1.linear_modeling_technical_vars.py"
uv run python "$module_dir/scripts/2.variance_decomposition.py"

# each calculation step is followed by its ggplot2 plotting step (R, run in the uvr project at the git root)
# uv run python "$module_dir/scripts/3.variate_importance.py"
uvr run "$module_dir/scripts/6.plot_variate_importance.r"

# uv run python "$module_dir/scripts/4.explore_linear_model_haystacks.py"
uvr run "$module_dir/scripts/7.plot_explore_linear_model_haystacks.r"

# uv run python "$module_dir/scripts/5.variate_class_upsets_and_clustermap.py"
 uvr run "$module_dir/scripts/8.plot_variate_class_upsets_and_clustermap.r"

echo "4.linear_modeling pipeline complete."
