list_of_packages <- c("ggplot2", "dplyr", "arrow", "RColorBrewer")
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

results_dir <- file.path(root_dir, "4.analysis", "results")
figures_dir <- file.path(root_dir, "4.analysis", "figures", "count_viability")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_parquet(file.path(results_dir, "count_viability_joined.parquet"))
df$Treatment <- factor(df$Treatment, levels = intersect(custom_treatment_order, unique(df$Treatment)))

plot_theme <- theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 9),
    strip.text = element_text(size = 10)
)

# --- per-patient: cell count vs viability, faceted by patient (3D only, one row per organoid-modality) ---
df_3d <- df %>% filter(modality == "3D")
p_3d_patient <- ggplot(df_3d, aes(x = min_max_viability, y = mean_cell_count, color = Treatment)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(values = custom_treatment_palette, na.value = "grey70") +
    facet_wrap(~Metadata_patient_tumor, scales = "free_y") +
    labs(title = "3D: mean cell count vs. viability, per patient",
         x = "Viability (min-max normalized)", y = "Mean cells per organoid") +
    plot_theme
ggsave(filename = file.path(figures_dir, "3D_count_vs_viability_per_patient.png"),
       plot = p_3d_patient, width = 12, height = 8, dpi = 600, units = "in")

p_3d_pooled <- ggplot(df_3d, aes(x = min_max_viability, y = mean_cell_count, color = Treatment)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
    scale_color_manual(values = custom_treatment_palette, na.value = "grey70") +
    labs(title = "3D pooled (all patients): mean cell count vs. viability",
         x = "Viability (min-max normalized)", y = "Mean cells per organoid") +
    plot_theme
ggsave(filename = file.path(figures_dir, "3D_count_vs_viability_pooled.png"),
       plot = p_3d_pooled, width = 9, height = 6, dpi = 600, units = "in")

# --- 2D: same, faceted by projection too ---
df_2d <- df %>% filter(modality == "2D")
p_2d_patient <- ggplot(df_2d, aes(x = min_max_viability, y = mean_cell_count, color = Treatment)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(values = custom_treatment_palette, na.value = "grey70") +
    facet_grid(projection ~ Metadata_patient_tumor, scales = "free_y") +
    labs(title = "2D: mean cell count (well-level avg.) vs. viability, per patient x projection method",
         x = "Viability (min-max normalized)", y = "Mean cells per organoid") +
    plot_theme + theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1), strip.text = element_text(size = 8))
ggsave(filename = file.path(figures_dir, "2D_count_vs_viability_per_patient.png"),
       plot = p_2d_patient, width = 16, height = 8, dpi = 600, units = "in")

p_2d_pooled <- ggplot(df_2d, aes(x = min_max_viability, y = mean_cell_count, color = Treatment)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
    scale_color_manual(values = custom_treatment_palette, na.value = "grey70") +
    facet_wrap(~projection) +
    labs(title = "2D pooled (all patients): mean cell count vs. viability, by projection method",
         x = "Viability (min-max normalized)", y = "Mean cells per organoid") +
    plot_theme
ggsave(filename = file.path(figures_dir, "2D_count_vs_viability_pooled.png"),
       plot = p_2d_pooled, width = 12, height = 5, dpi = 600, units = "in")

cat("Wrote 4 figures to", figures_dir, "\n")
