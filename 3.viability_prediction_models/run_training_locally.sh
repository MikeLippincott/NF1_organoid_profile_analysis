#!/bin/bash

set -euo pipefail

git_root=$(git rev-parse --show-toplevel)
if [ -z "$git_root" ]; then
    echo "Error: Could not find the git root directory."
    exit 1
fi

jupyter nbconvert --to=script --FilesWriter.build_directory="$git_root"/3.viability_prediction_models/scripts/ "$git_root"/3.viability_prediction_models/notebooks/*.ipynb

cd "$git_root"/3.viability_prediction_models/scripts

uv run python 0.pre-processing_profiles_for_viability_models.py
uv run python 1.viability_prediction.py

echo "Training complete."
