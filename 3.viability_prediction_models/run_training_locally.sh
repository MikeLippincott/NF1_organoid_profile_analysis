#!/bin/bash

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/ notebooks/*.ipynb

cd scripts || exit 1

uv run python viability_prediction.py

cd ../ || exit 1

echo "Training complete."
