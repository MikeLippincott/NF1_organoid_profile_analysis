#!/bin/bash

git_root=$(git rev-parse --show-toplevel)
if [ -z "$git_root" ]; then
    echo "Error: Could not find the git root directory."
    exit 1
fi

patient_array_file_path="$git_root/data/patient_IDs.txt"
# read the patient IDs from the file into an array
if [[ -f "$patient_array_file_path" ]]; then
    readarray -t patient_array < "$patient_array_file_path"
else
    echo "Error: File $patient_array_file_path does not exist."
    exit 1
fi


# deactivate any existing conda environment
conda deactivate
# deactivate any existing venv environment
deactivate 2>/dev/null

uv run --project "$git_root" jupyter nbconvert --to=script --FilesWriter.build_directory="$git_root"/1.EDA/scripts/ "$git_root"/1.EDA/notebooks/*.ipynb



# Run the Python and R scripts using uv
# uv run --project "$git_root" python "$git_root"/1.EDA/scripts/0.generate_umap.py
# conda run -n gff_figure_env Rscript "$git_root"/1.EDA/scripts/1.plot_umap.r
# uv run --project "$git_root" python "$git_root"/1.EDA/scripts/2.generate_pca.py
# conda run -n gff_figure_env Rscript "$git_root"/1.EDA/scripts/3.plot_pca.r
for patient in "${patient_array[@]}"; do
    echo "Processing patient: $patient"
    conda run -n gff_figure_env Rscript "$git_root"/1.EDA/scripts/4.consensus_profiles.r --patient "$patient"
done
# uv run --project "$git_root" python "$git_root"/1.EDA/scripts/5.calculate_correlation_matrix.py
# conda run -n gff_figure_env Rscript "$git_root"/1.EDA/scripts/6.plot_correlation_matrix.r
# uv run --project "$git_root" python "$git_root"/1.EDA/scripts/7.generate_cell_counts.py
# conda run -n gff_figure_env Rscript "$git_root"/1.EDA/scripts/8.plot_cell_counts.r
conda run -n gff_figure_env Rscript "$git_root"/1.EDA/scripts/9.plot_volume_area_distributions.r
uv run --project "$git_root" python "$git_root"/1.EDA/scripts/10.calculate_count_viability_join.py
conda run -n gff_figure_env Rscript "$git_root"/1.EDA/scripts/11.plot_count_viability_scatter.r
conda run -n gff_figure_env Rscript "$git_root"/1.EDA/scripts/12.plot_volume_area_vs_count.r
uv run --project "$git_root" python "$git_root"/1.EDA/scripts/13.calculate_neighbor_features.py
conda run -n gff_figure_env Rscript "$git_root"/1.EDA/scripts/14.plot_neighbor_features.r
uv run --project "$git_root" python "$git_root"/1.EDA/scripts/15.calculate_intensity_summary.py
conda run -n gff_figure_env Rscript "$git_root"/1.EDA/scripts/16.plot_intensity_facets.r
uv run --project "$git_root" python "$git_root"/1.EDA/scripts/17.calculate_area_volume_by_patient_treatment.py
conda run -n gff_figure_env Rscript "$git_root"/1.EDA/scripts/18.plot_area_vs_volume.r
