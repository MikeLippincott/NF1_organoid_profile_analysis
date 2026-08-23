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

load_kind <- function(kind) {
    area <- read_parquet(file.path(results_dir, paste0("area_2D_", kind, ".parquet")))
    volume <- read_parquet(file.path(results_dir, paste0("volume_3D_", kind, ".parquet")))
    area$Metadata_treatment <- factor(area$Metadata_treatment,
        levels = intersect(custom_treatment_order, unique(area$Metadata_treatment)))
    volume$Metadata_treatment <- factor(volume$Metadata_treatment,
        levels = intersect(custom_treatment_order, unique(volume$Metadata_treatment)))
    list(area = area, volume = volume)
}

patient_palette_for <- function(patient_tumors) {
    setNames(tab20_palette_for_patients[seq_along(unique(patient_tumors))], unique(patient_tumors))
}

# full cross-join within patient_tumor x treatment x dose (no shared organoid/cell ID)
cross_join_area_volume <- function(area_df, volume_df, projection_name) {
    area_df %>%
        filter(projection == projection_name) %>%
        inner_join(volume_df, by = c("Metadata_patient_tumor", "Metadata_treatment", "Metadata_dose"),
                   relationship = "many-to-many")
}

# --- per kind (organoid / single-cell) x per 2D projection method:
# heatmap faceted by patient_tumor, scatter faceted by patient_tumor (color =
# treatment), scatter faceted by treatment (color = patient_tumor) - all with
# dose as shape where points are drawn, and a y = x reference line ---
n_figures <- 0
for (kind in c("organoid", "cell")) {
    label <- if (kind == "organoid") "Organoid" else "Single-cell"
    d <- load_kind(kind)

    for (proj in unique(d$area$projection)) {
        joined <- cross_join_area_volume(d$area, d$volume, proj)
        patient_palette <- patient_palette_for(joined$Metadata_patient_tumor)
        axis_range <- range(c(joined$area, joined$volume), na.rm = TRUE)

        p_heatmap <- (
            ggplot(joined, aes(x = volume, y = area))
            + geom_bin2d(bins = 60)
            + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 0.5)
            + scale_fill_viridis_c(trans = "log10")
            + coord_equal(xlim = axis_range, ylim = axis_range)
            + facet_wrap(~Metadata_patient_tumor)
            + labs(
                title = paste0(label, " area (2D, ", proj, ") vs. volume (3D), every point, by patient_tumor"),
                x = paste0(label, " volume, 3D (z-scored)"),
                y = paste0(label, " area, 2D ", proj, " (z-scored)"),
                fill = "Count (log10)"
            )
            + plot_theme
        )
        ggsave(filename = file.path(figures_dir, paste0(kind, "_area_vs_volume_", proj, "_heatmap_by_patient.png")),
               plot = p_heatmap, width = 14, height = 10, dpi = 600, units = "in")
        n_figures <- n_figures + 1

        p_by_patient_color_treatment <- (
            ggplot(joined, aes(x = volume, y = area, color = Metadata_treatment, shape = factor(Metadata_dose)))
            + geom_point(size = 0.6, alpha = 0.15)
            + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5)
            + scale_color_manual(values = custom_treatment_palette, na.value = "grey70")
            + coord_equal(xlim = axis_range, ylim = axis_range)
            + facet_wrap(~Metadata_patient_tumor)
            + labs(
                title = paste0(label, " area (2D, ", proj, ") vs. volume (3D), every point, by patient_tumor, colored by treatment"),
                x = paste0(label, " volume, 3D (z-scored)"),
                y = paste0(label, " area, 2D ", proj, " (z-scored)"),
                color = "Treatment", shape = "Dose"
            )
            + plot_theme
            + guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))
        )
        ggsave(filename = file.path(figures_dir, paste0(kind, "_area_vs_volume_", proj, "_by_patient_color_treatment.png")),
               plot = p_by_patient_color_treatment, width = 14, height = 10, dpi = 600, units = "in")
        n_figures <- n_figures + 1

        p_by_treatment_color_patient <- (
            ggplot(joined, aes(x = volume, y = area, color = Metadata_patient_tumor, shape = factor(Metadata_dose)))
            + geom_point(size = 0.6, alpha = 0.15)
            + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5)
            + scale_color_manual(values = patient_palette)
            + coord_equal(xlim = axis_range, ylim = axis_range)
            + facet_wrap(~Metadata_treatment)
            + labs(
                title = paste0(label, " area (2D, ", proj, ") vs. volume (3D), every point, by treatment, colored by patient_tumor"),
                x = paste0(label, " volume, 3D (z-scored)"),
                y = paste0(label, " area, 2D ", proj, " (z-scored)"),
                color = "Patient", shape = "Dose"
            )
            + plot_theme
            + guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))
        )
        ggsave(filename = file.path(figures_dir, paste0(kind, "_area_vs_volume_", proj, "_by_treatment_color_patient.png")),
               plot = p_by_treatment_color_patient, width = 18, height = 12, dpi = 600, units = "in")
        n_figures <- n_figures + 1

        rm(joined)
        gc()
    }
}
cat("Wrote", n_figures, "per-projection figures to", figures_dir, "\n")

# --- pooled across patients: one heatmap + one scatter-by-treatment + one
# scatter-by-patient, each faceted by projection method, per kind ---
n_pooled_figures <- 0
for (kind in c("organoid", "cell")) {
    label <- if (kind == "organoid") "Organoid" else "Single-cell"
    d <- load_kind(kind)

    joined_all <- bind_rows(lapply(unique(d$area$projection), function(proj) {
        cross_join_area_volume(d$area, d$volume, proj) %>% mutate(projection = proj)
    }))
    patient_palette <- patient_palette_for(joined_all$Metadata_patient_tumor)
    axis_range <- range(c(joined_all$area, joined_all$volume), na.rm = TRUE)

    p_pooled_heatmap <- (
        ggplot(joined_all, aes(x = volume, y = area))
        + geom_bin2d(bins = 60)
        + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 0.5)
        + scale_fill_viridis_c(trans = "log10")
        + coord_equal(xlim = axis_range, ylim = axis_range)
        + facet_wrap(~projection)
        + labs(
            title = paste0(label, " area (2D) vs. volume (3D), every point, pooled across patients, by projection method"),
            x = paste0(label, " volume, 3D (z-scored)"),
            y = paste0(label, " area, 2D (z-scored)"),
            fill = "Count (log10)"
        )
        + plot_theme
    )
    ggsave(filename = file.path(figures_dir, paste0(kind, "_area_vs_volume_pooled_heatmap_by_projection.png")),
           plot = p_pooled_heatmap, width = 16, height = 6, dpi = 600, units = "in")
    n_pooled_figures <- n_pooled_figures + 1

    p_pooled_color_treatment <- (
        ggplot(joined_all, aes(x = volume, y = area, color = Metadata_treatment, shape = factor(Metadata_dose)))
        + geom_point(size = 0.5, alpha = 0.1)
        + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5)
        + scale_color_manual(values = custom_treatment_palette, na.value = "grey70")
        + coord_equal(xlim = axis_range, ylim = axis_range)
        + facet_wrap(~projection)
        + labs(
            title = paste0(label, " area (2D) vs. volume (3D), every point, pooled, colored by treatment"),
            x = paste0(label, " volume, 3D (z-scored)"),
            y = paste0(label, " area, 2D (z-scored)"),
            color = "Treatment", shape = "Dose"
        )
        + plot_theme
        + guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))
    )
    ggsave(filename = file.path(figures_dir, paste0(kind, "_area_vs_volume_pooled_by_projection_color_treatment.png")),
           plot = p_pooled_color_treatment, width = 16, height = 6, dpi = 600, units = "in")
    n_pooled_figures <- n_pooled_figures + 1

    p_pooled_color_patient <- (
        ggplot(joined_all, aes(x = volume, y = area, color = Metadata_patient_tumor, shape = factor(Metadata_dose)))
        + geom_point(size = 0.5, alpha = 0.1)
        + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5)
        + scale_color_manual(values = patient_palette)
        + coord_equal(xlim = axis_range, ylim = axis_range)
        + facet_wrap(~projection)
        + labs(
            title = paste0(label, " area (2D) vs. volume (3D), every point, pooled, colored by patient_tumor"),
            x = paste0(label, " volume, 3D (z-scored)"),
            y = paste0(label, " area, 2D (z-scored)"),
            color = "Patient", shape = "Dose"
        )
        + plot_theme
        + guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))
    )
    ggsave(filename = file.path(figures_dir, paste0(kind, "_area_vs_volume_pooled_by_projection_color_patient.png")),
           plot = p_pooled_color_patient, width = 16, height = 6, dpi = 600, units = "in")
    n_pooled_figures <- n_pooled_figures + 1

    rm(joined_all)
    gc()
}

cat("Wrote", n_figures + n_pooled_figures, "total figures to", figures_dir, "\n")
