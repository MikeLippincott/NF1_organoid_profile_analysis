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

results_dir <- file.path(root_dir, "1.EDA", "results", "area_vs_volume")
figures_dir <- file.path(root_dir, "1.EDA", "figures", "area_vs_volume")
figures_dir_rel <- file.path("1.EDA", "figures", "area_vs_volume")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

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

metric_palette <- c("Area (2D)" = "#3D7DCC", "Volume (3D)" = "#E06C75")

load_kind <- function(kind) {
    area <- read_parquet(file.path(results_dir, paste0("area_2D_", kind, ".parquet")))
    volume <- read_parquet(file.path(results_dir, paste0("volume_3D_", kind, ".parquet")))
    area$Metadata_treatment <- factor(area$Metadata_treatment,
        levels = intersect(custom_treatment_order, unique(area$Metadata_treatment)))
    volume$Metadata_treatment <- factor(volume$Metadata_treatment,
        levels = intersect(custom_treatment_order, unique(volume$Metadata_treatment)))
    list(area = area, volume = volume)
}

# stack area (2D, one projection) and volume (3D) into one long table - no
# join needed, since we're comparing the two distributions' shapes, not
# pairing individual records
stack_area_volume <- function(area_df, volume_df, projection_name) {
    bind_rows(
        area_df %>%
            filter(projection == projection_name) %>%
            transmute(Metadata_patient_tumor, Metadata_treatment, Metadata_dose,
                      value = area, metric = "Area (2D)"),
        volume_df %>%
            transmute(Metadata_patient_tumor, Metadata_treatment, Metadata_dose,
                      value = volume, metric = "Volume (3D)")
    )
}

# --- per kind (organoid / single-cell) x per 2D projection method: area and
# volume distributions side-by-side, faceted by patient_tumor and pooled by
# treatment ---
n_figures <- 0
for (kind in c("organoid", "cell")) {
    label <- if (kind == "organoid") "Organoid" else "Single-cell"
    d <- load_kind(kind)

    for (proj in unique(d$area$projection)) {
        long_df <- stack_area_volume(d$area, d$volume, proj)
        long_df$metric <- factor(long_df$metric, levels = c("Area (2D)", "Volume (3D)"))

        p_by_patient <- (
            ggplot(long_df, aes(x = Metadata_patient_tumor, y = value, fill = metric))
            + geom_violin(alpha = 0.6, trim = TRUE, position = position_dodge(width = 0.8), width = 0.8)
            + geom_boxplot(width = 0.15, alpha = 0.85, outlier.size = 0.3, outlier.alpha = 0.3,
                            position = position_dodge(width = 0.8))
            + scale_fill_manual(values = metric_palette)
            + labs(
                title = paste0(label, ": area (2D, ", proj, ") vs. volume (3D) distributions, by patient"),
                x = "Patient", y = "Value (z-scored)", fill = "Metric"
            )
            + plot_theme
        )
        ggsave(filename = file.path(figures_dir, paste0(kind, "_area_vs_volume_", proj, "_by_patient.png")),
               plot = p_by_patient, width = 12, height = 7, dpi = 600, units = "in")
        n_figures <- n_figures + 1

        p_by_treatment <- (
            ggplot(long_df, aes(x = Metadata_treatment, y = value, fill = metric))
            + geom_violin(alpha = 0.6, trim = TRUE, position = position_dodge(width = 0.8), width = 0.8)
            + geom_boxplot(width = 0.15, alpha = 0.85, outlier.size = 0.3, outlier.alpha = 0.3,
                            position = position_dodge(width = 0.8))
            + scale_fill_manual(values = metric_palette)
            + labs(
                title = paste0(label, " pooled (all patients): area (2D, ", proj, ") vs. volume (3D) distributions, by treatment"),
                x = "Treatment", y = "Value (z-scored)", fill = "Metric"
            )
            + plot_theme
        )
        ggsave(filename = file.path(figures_dir, paste0(kind, "_area_vs_volume_", proj, "_by_treatment.png")),
               plot = p_by_treatment, width = 14, height = 7, dpi = 600, units = "in")
        n_figures <- n_figures + 1
    }
}
cat("Wrote", n_figures, "figures to", figures_dir_rel, "\n")
