#!/bin/bash

set -euo pipefail

git_root=$(git rev-parse --show-toplevel)
if [ -z "$git_root" ]; then
    echo "Error: Could not find the git root directory."
    exit 1
fi

jupyter nbconvert --to=script --FilesWriter.build_directory="$git_root"/4.linear_modeling/scripts/ "$git_root"/4.linear_modeling/notebooks/*.ipynb

cd "$git_root"/4.linear_modeling/scripts

uv run python caluclate_well_manhattan_distance.py
uv run python 0.linear_modeling.py
uv run python 1.linear_modeling_technical_vars.py

echo "Linear modeling complete."
