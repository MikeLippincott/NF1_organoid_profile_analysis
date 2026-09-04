list_of_packages <- c("ggplot2", "dplyr", "arrow", "tidyr", "RColorBrewer", "ggrastr", "patchwork")
for (package in list_of_packages) {
    suppressPackageStartupMessages(
        suppressWarnings(
            library(package, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
        )
    )
}

find_git_root <- function() {
    cwd <- getwd()
    if (dir.exists(file.path(cwd, ".git"))) {
        return(cwd)
    }
    current_path <- cwd
    while (dirname(current_path) != current_path) {
        parent_path <- dirname(current_path)
        if (dir.exists(file.path(parent_path, ".git"))) {
            return(parent_path)
        }
        current_path <- parent_path
    }
    stop("No Git root directory found.")
}

root_dir <- find_git_root()
source(file.path(root_dir, "utils", "r_plot_themes.r"))

results_dir <- file.path(root_dir, "1.EDA", "results", "neighbors")
figures_dir <- file.path(root_dir, "1.EDA", "figures", "neighbors")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

plot_theme <- theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, size = 13),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 9),
    legend.position = "none"
)

# --- local crowding: 3D nuc (single-cell) and 3D org (organoid), plus 3D org density ---
# Patient/treatment are shown as separate marginal panels (not a full
# patient x treatment grid, which would be ~250 tiny panels) so it's possible
# to see which factor accounts for more of the spread in each metric.
nuc_3d <- read_parquet(file.path(results_dir, "nuclei_neighbors_3D.parquet"))
# 3D organoid neighbor counts are derived in notebook 13 (bounding-sphere
# touch test within each organoid's own image), since no Organoid_Neighbors_*
# column exists in the 3D profiles the way it does in 2D.
org_3d <- read_parquet(file.path(results_dir, "organoid_neighbors_3D.parquet"))
# Organoid cell density (cells / organoid volume in um^3), also derived in
# notebook 13 since it isn't a profile column either.
org_density_3d <- read_parquet(file.path(results_dir, "organoid_density_3D.parquet"))

nuc_3d <- nuc_3d %>%
    mutate(Metadata_Experiment_Treatment = factor(Metadata_Experiment_Treatment, levels = intersect(custom_treatment_order, unique(Metadata_Experiment_Treatment))))
org_3d <- org_3d %>%
    mutate(Metadata_Experiment_Treatment = factor(Metadata_Experiment_Treatment, levels = intersect(custom_treatment_order, unique(Metadata_Experiment_Treatment))))
org_density_3d <- org_density_3d %>%
    mutate(Metadata_Experiment_Treatment = factor(Metadata_Experiment_Treatment, levels = intersect(custom_treatment_order, unique(Metadata_Experiment_Treatment))))

patient_treatment_panel <- function(df, y_col, y_label, title_prefix, y_trans = "identity") {
    p_patient <- (
        ggplot(df, aes(x = Metadata_Biology_PatientTumor, y = .data[[y_col]], fill = Metadata_Biology_PatientTumor))
        + rasterise(geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.2), dpi = 300)
        + scale_y_continuous(trans = y_trans)
        + labs(title = paste0(title_prefix, " - by patient (treatments pooled)"),
               x = "Patient", y = y_label)
        + plot_theme
    )
    p_treatment <- (
        ggplot(df, aes(x = Metadata_Experiment_Treatment, y = .data[[y_col]], fill = Metadata_Experiment_Treatment))
        + rasterise(geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.2), dpi = 300)
        + scale_fill_manual(values = custom_treatment_palette, na.value = "grey70")
        + scale_y_continuous(trans = y_trans)
        + labs(title = paste0(title_prefix, " - by treatment (patients pooled)"),
               x = "Treatment", y = y_label)
        + plot_theme
    )
    p_patient + p_treatment
}

p_crowding_nuc <- patient_treatment_panel(
    nuc_3d, "Metadata_Neighbors_NeighborsCountAdjacent", "Adjacent neighbor count",
    "Nuclei number of neighbors"
)
p_crowding_organoid <- patient_treatment_panel(
    org_3d, "Organoid_Neighbors_NumberOfNeighbors_Adjacent", "Adjacent neighbor count",
    "Organoids number of neighbors"
)
p_density_organoid <- patient_treatment_panel(
    org_density_3d, "Organoid_CellDensity_CellsPerUm3", "Cell density (cells / um^3, pseudo-log scale)",
    "Organoids cell density", y_trans = scales::pseudo_log_trans(sigma = 1e-4)
)

pdf_path <- file.path(figures_dir, "3D_neighbor_organoid_story.pdf")
pdf(pdf_path, width = 16, height = 9, onefile = TRUE)
print(p_crowding_nuc)
print(p_crowding_organoid)
print(p_density_organoid)
invisible(dev.off())
