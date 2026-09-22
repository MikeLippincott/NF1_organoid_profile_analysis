#!/bin/bash

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/ notebooks/*.ipynb

cd scripts || exit 1

echo "Running pre-processing..."
uv run python 0.pre-processing_profiles_for_viability_models.py

echo "Running viability prediction model training..."
uv run python 1.viability_prediction.py

echo "Running model results visualization..."
uvr run 2.visualize_model_results.r

echo "Visualizing single-cell data..."
uv run python 3.single-cell_visualizations.py

echo "Running multi-panel summary figure..."
uvr run 4.multi_panel_summary_figure.r

cd ../ || exit 1

echo "Training complete."
