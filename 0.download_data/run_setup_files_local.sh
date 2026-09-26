#!/bin/bash

set -euo pipefail

git_root=$(git rev-parse --show-toplevel)
if [ -z "$git_root" ]; then
    echo "Error: Could not find the git root directory."
    exit 1
fi

jupyter nbconvert --to=script --FilesWriter.build_directory="$git_root"/0.download_data/scripts/ "$git_root"/0.download_data/notebooks/*.ipynb

cd "$git_root"/0.download_data/scripts

uv run python setup_files.py

echo "Data setup complete."
