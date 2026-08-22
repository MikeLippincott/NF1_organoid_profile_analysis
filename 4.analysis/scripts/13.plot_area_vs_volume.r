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

results_dir <- file.path(root_dir, "4.analysis", "results", "area_vs_volume")
figures_dir <- file.path(root_dir, "4.analysis", "figures", "area_vs_volume")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_parquet(file.path(results_dir, "area_vs_volume_by_patient_treatment.parquet"))
df$Metadata_treatment <- factor(df$Metadata_treatment,
                                  levels = intersect(custom_treatment_order, unique(df$Metadata_treatment)))
patients <- sort(unique(df$Metadata_patient))
patient_colors <- setNames(tab20_palette_for_patients[seq_along(patients)], patients)

plot_theme <- theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, size = 13),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(size = 8),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 8)
)

# --- one figure per 2D projection method: faceted by patient, colored by treatment ---
for (proj in unique(df$projection)) {
    d_proj <- df %>% filter(projection == proj)
    p <- ggplot(d_proj, aes(x = mean_volume, y = mean_area, color = Metadata_treatment)) +
        geom_point(size = 2, alpha = 0.85) +
        scale_color_manual(values = custom_treatment_palette, na.value = "grey70") +
        facet_wrap(~Metadata_patient, scales = "free") +
        labs(
            title = paste0("Organoid area (2D, ", proj, ") vs. volume (3D), by patient and treatment"),
            x = "Mean organoid volume, 3D (z-scored)",
            y = paste0("Mean organoid area, 2D ", proj, " (z-scored)"),
            color = "Treatment"
        ) +
        plot_theme
    ggsave(filename = file.path(figures_dir, paste0("area_vs_volume_by_patient_", proj, ".png")),
           plot = p, width = 12, height = 9, dpi = 600, units = "in")
}

# --- pooled: one figure faceted by projection, colored by patient ---
p_pooled <- ggplot(df, aes(x = mean_volume, y = mean_area, color = Metadata_patient)) +
    geom_point(size = 2, alpha = 0.85) +
    scale_color_manual(values = patient_colors) +
    facet_wrap(~projection) +
    labs(
        title = "Organoid area (2D) vs. volume (3D), pooled across treatments, by projection method",
        x = "Mean organoid volume, 3D (z-scored)",
        y = "Mean organoid area, 2D (z-scored)",
        color = "Patient"
    ) +
    plot_theme
ggsave(filename = file.path(figures_dir, "area_vs_volume_pooled_by_projection.png"),
       plot = p_pooled, width = 14, height = 6, dpi = 600, units = "in")

cat("Wrote", length(unique(df$projection)) + 1, "figures to", figures_dir, "\n")
