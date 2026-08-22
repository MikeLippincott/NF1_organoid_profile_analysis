#!/bin/bash

git_root=$(git rev-parse --show-toplevel)

ENV_PATH="$git_root/.venv"


if [[ -e "$ENV_PATH" ]]; then
    echo "Using existing virtual environment at $ENV_PATH"
else
    mkdir -p "$ENV_PATH"
    echo "Created virtual environment directory at $ENV_PATH"
fi

# run twice to ensure we are in a clean environment
# and not accidentally using an existing one
# try but will not fail if not in a conda environment
loops=2
for _ in $(seq 1 $loops); do
    if command -v conda >/dev/null 2>&1; then
        conda deactivate >/dev/null 2>&1 || true
    fi
    deactivate >/dev/null 2>&1 || true
done
unset VIRTUAL_ENV

rm -f uv.lock
rm -rf .venv

uv venv
uv sync

# shellcheck disable=SC1091
source .venv/bin/activate
# Use RELATIVE path - simple and reliable
uv pip install -e ./utils

# Link R/Rscript from the gff_figure_env conda env (see environments/r_env.yml)
# into the uv venv's bin dir, so `uv run Rscript ...` resolves to an R
# installation that actually has the required packages (ggplot2, dplyr,
# arrow, tidyr, readr, RColorBrewer, ...).
if command -v conda >/dev/null 2>&1; then
    R_ENV_PATH=$(conda env list | awk '$1 == "gff_figure_env" {print $NF}')
    if [[ -n "$R_ENV_PATH" && -x "$R_ENV_PATH/bin/Rscript" ]]; then
        ln -sf "$R_ENV_PATH/bin/R" "$ENV_PATH/bin/R"
        ln -sf "$R_ENV_PATH/bin/Rscript" "$ENV_PATH/bin/Rscript"
        echo "Linked R/Rscript from gff_figure_env ($R_ENV_PATH) into $ENV_PATH/bin"
    else
        echo "Warning: gff_figure_env not found or missing Rscript; run" \
             "'mamba env create -f environments/r_env.yml' first" \
             "(see environments/local_create_envs.sh)."
    fi
else
    echo "Warning: conda not found; cannot link R/Rscript into the uv venv."
fi

echo "Virtual environment setup complete. To activate it, run:"
echo "source $ENV_PATH/bin/activate"
