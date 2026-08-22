#!/bin/bash

# Run all analysis scripts using uv env
# set -e

git_root=$(git rev-parse --show-toplevel)
# convert to scripts
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/ notebooks/*.ipynb

# Activate uv environment
# shellcheck disable=SC1091
source "$git_root/.venv/bin/activate"

# Run the analysis scripts and plot scripts
uv run python "$git_root"/4.analysis/scripts/0.harmonize_metadata.py
uv run python "$git_root"/4.analysis/scripts/1.calculate_descriptive_stats.py
uv run python "$git_root"/4.analysis/scripts/2.test_feature_significance.py
uv run Rscript "$git_root"/4.analysis/scripts/3.plot_descriptive_stats_summary.r
uv run Rscript "$git_root"/4.analysis/scripts/4.plot_volume_area_distributions.r
uv run python "$git_root"/4.analysis/scripts/5.calculate_count_viability_join.py
uv run Rscript "$git_root"/4.analysis/scripts/6.plot_count_viability_scatter.r
uv run Rscript "$git_root"/4.analysis/scripts/7.plot_volume_area_vs_count.r
uv run python "$git_root"/4.analysis/scripts/8.calculate_neighbor_features.py
uv run Rscript "$git_root"/4.analysis/scripts/9.plot_neighbor_features.r
uv run python "$git_root"/4.analysis/scripts/10.calculate_intensity_summary.py
uv run Rscript "$git_root"/4.analysis/scripts/11.plot_intensity_facets.r
uv run python "$git_root"/4.analysis/scripts/12.calculate_area_volume_by_patient_treatment.py
uv run Rscript "$git_root"/4.analysis/scripts/13.plot_area_vs_volume.r
uv run python "$git_root"/4.analysis/scripts/14.calculate_cells_per_fov_viability.py
uv run Rscript "$git_root"/4.analysis/scripts/15.plot_cells_per_fov_viability.r

echo "All analysis scripts completed successfully!"


