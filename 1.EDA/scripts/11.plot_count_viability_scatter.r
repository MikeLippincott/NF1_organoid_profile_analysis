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

results_dir <- file.path(root_dir, "1.EDA", "results")
figures_dir <- file.path(root_dir, "1.EDA", "figures", "count_viability")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_parquet(file.path(results_dir, "count_viability_joined.parquet"))
df$Treatment <- factor(df$Treatment, levels = intersect(custom_treatment_order, unique(df$Treatment)))

plot_theme <- (
    theme_bw()
    + theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 9),
        strip.text = element_text(size = 10)
    )
)

# --- per-patient: mean cells per organoid vs. viability, faceted by patient ---
p_patient <- (
    ggplot(df, aes(x = min_max_viability, y = mean_cell_count, color = Treatment, shape = factor(Metadata_dose)))
    + geom_point(size = 2, alpha = 0.8)
    + scale_color_manual(values = custom_treatment_palette, na.value = "grey70")
    + facet_wrap(~Metadata_patient_tumor, scales = "free_y")
    + labs(
        title = "3D: mean cells per organoid vs. viability, per patient",
        x = "Viability (min-max normalized)", y = "Mean cells per organoid", shape = "Dose"
    )
    + plot_theme
)
ggsave(
    filename = file.path(figures_dir, "3D_count_vs_viability_per_patient.png"),
    plot = p_patient, width = 12, height = 8, dpi = 600, units = "in"
)

p_pooled <- (
    ggplot(df, aes(x = min_max_viability, y = mean_cell_count, color = Treatment, shape = factor(Metadata_dose)))
    + geom_point(size = 2, alpha = 0.8)
    + geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5)
    + scale_color_manual(values = custom_treatment_palette, na.value = "grey70")
    + labs(
        title = "3D pooled (all patients): mean cells per organoid vs. viability",
        x = "Viability (min-max normalized)", y = "Mean cells per organoid", shape = "Dose"
    )
    + plot_theme
)
ggsave(
    filename = file.path(figures_dir, "3D_count_vs_viability_pooled.png"),
    plot = p_pooled, width = 9, height = 6, dpi = 600, units = "in"
)

cat("Wrote 2 figures to", figures_dir, "\n")


df_norm <- read_parquet(file.path(results_dir, "count_norm_viability_joined.parquet"))
df_norm$Treatment <- factor(df_norm$Treatment, levels = intersect(custom_treatment_order, unique(df_norm$Treatment)))

patient_palette <- setNames(
    tab20_palette_for_patients[seq_along(unique(df_norm$Metadata_patient_tumor))],
    unique(df_norm$Metadata_patient_tumor)
)

df_norm_3d <- df_norm %>% filter(modality == "3D")
df_norm_2d <- df_norm %>% filter(modality == "2D")

count_metrics <- list(
    list(col = "total_cell_count_norm", lab = "Total cells per treatment (FOV-normalized)", prefix = "cell_count_norm"),
    list(col = "total_organoid_count_norm", lab = "Total organoids per treatment (FOV-normalized)", prefix = "organoid_count_norm")
)

# --- 3D: FOV-normalized cell/organoid counts vs. viability, 3 facet variants ---
for (m in count_metrics) {
    p_by_patient <- (
        ggplot(df_norm_3d, aes(x = min_max_viability, y = .data[[m$col]], color = Treatment, shape = factor(Metadata_dose)))
        + geom_point(size = 2, alpha = 0.8)
        + scale_color_manual(values = custom_treatment_palette, na.value = "grey70")
        + facet_wrap(~Metadata_patient_tumor, scales = "free_y")
        + labs(
            title = paste0("3D: ", m$lab, " vs. viability, by patient"),
            x = "Viability (min-max normalized)", y = m$lab, shape = "Dose"
        )
        + plot_theme
    )
    ggsave(
        filename = file.path(figures_dir, paste0("3D_", m$prefix, "_vs_viability_by_patient.png")),
        plot = p_by_patient, width = 12, height = 8, dpi = 600, units = "in"
    )

    p_by_treatment <- (
        ggplot(df_norm_3d, aes(x = min_max_viability, y = .data[[m$col]], color = Metadata_patient_tumor, shape = factor(Metadata_dose)))
        + geom_point(size = 2, alpha = 0.8)
        + scale_color_manual(values = patient_palette)
        + facet_wrap(~Treatment, scales = "free_y")
        + labs(
            title = paste0("3D: ", m$lab, " vs. viability, by treatment"),
            x = "Viability (min-max normalized)", y = m$lab, color = "Patient", shape = "Dose"
        )
        + plot_theme
        + theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1))
    )
    ggsave(
        filename = file.path(figures_dir, paste0("3D_", m$prefix, "_vs_viability_by_treatment.png")),
        plot = p_by_treatment, width = 14, height = 10, dpi = 600, units = "in"
    )

    p_by_both <- (
        ggplot(df_norm_3d, aes(x = min_max_viability, y = .data[[m$col]], shape = factor(Metadata_dose)))
        + geom_point(size = 1.5, alpha = 0.7, color = "#3D7DCC")
        + facet_grid(Treatment ~ Metadata_patient_tumor, scales = "free_y")
        + labs(
            title = paste0("3D: ", m$lab, " vs. viability, by patient and treatment"),
            x = "Viability (min-max normalized)", y = m$lab, shape = "Dose"
        )
        + plot_theme
        + theme(
            axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
            strip.text.y = element_text(size = 6, angle = 0)
        )
    )
    ggsave(
        filename = file.path(figures_dir, paste0("3D_", m$prefix, "_vs_viability_by_patient_and_treatment.png")),
        plot = p_by_both, width = 16, height = 20, dpi = 600, units = "in"
    )
}
cat("Wrote", 3 * length(count_metrics), "3D figures to", figures_dir, "\n")

# --- 2D: FOV-normalized cell/organoid counts vs. viability, all 2D methods
# (max_projection, middle_slice, middle_n_slice), 3 facet variants ---
n_2d_figures <- 0
for (m in count_metrics) {
    p_by_patient <- (
        ggplot(df_norm_2d, aes(x = min_max_viability, y = .data[[m$col]], color = Treatment, shape = factor(Metadata_dose)))
        + geom_point(size = 1.5, alpha = 0.7)
        + scale_color_manual(values = custom_treatment_palette, na.value = "grey70")
        + facet_grid(projection ~ Metadata_patient_tumor, scales = "free_y")
        + labs(
            title = paste0("2D: ", m$lab, " vs. viability, by patient x projection method"),
            x = "Viability (min-max normalized)", y = m$lab, shape = "Dose"
        )
        + plot_theme
        + theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1), strip.text = element_text(size = 7))
    )
    ggsave(
        filename = file.path(figures_dir, paste0("2D_", m$prefix, "_vs_viability_by_patient.png")),
        plot = p_by_patient, width = 18, height = 8, dpi = 600, units = "in"
    )
    n_2d_figures <- n_2d_figures + 1

    p_by_treatment <- (
        ggplot(df_norm_2d, aes(x = min_max_viability, y = .data[[m$col]], color = Metadata_patient_tumor, shape = factor(Metadata_dose)))
        + geom_point(size = 1.5, alpha = 0.7)
        + scale_color_manual(values = patient_palette)
        + facet_grid(projection ~ Treatment, scales = "free_y")
        + labs(
            title = paste0("2D: ", m$lab, " vs. viability, by treatment x projection method"),
            x = "Viability (min-max normalized)", y = m$lab, color = "Patient", shape = "Dose"
        )
        + plot_theme
        + theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1), strip.text = element_text(size = 7))
    )
    ggsave(
        filename = file.path(figures_dir, paste0("2D_", m$prefix, "_vs_viability_by_treatment.png")),
        plot = p_by_treatment, width = 20, height = 8, dpi = 600, units = "in"
    )
    n_2d_figures <- n_2d_figures + 1

    for (proj in unique(df_norm_2d$projection)) {
        df_proj <- df_norm_2d %>% filter(projection == proj)
        p_both <- (
            ggplot(df_proj, aes(x = min_max_viability, y = .data[[m$col]], shape = factor(Metadata_dose)))
            + geom_point(size = 1.2, alpha = 0.7, color = "#3D7DCC")
            + facet_grid(Treatment ~ Metadata_patient_tumor, scales = "free_y")
            + labs(
                title = paste0("2D (", proj, "): ", m$lab, " vs. viability, by patient and treatment"),
                x = "Viability (min-max normalized)", y = m$lab, shape = "Dose"
            )
            + plot_theme
            + theme(
                axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
                strip.text.y = element_text(size = 6, angle = 0)
            )
        )
        ggsave(
            filename = file.path(figures_dir, paste0("2D_", proj, "_", m$prefix, "_vs_viability_by_patient_and_treatment.png")),
            plot = p_both, width = 16, height = 20, dpi = 600, units = "in"
        )
        n_2d_figures <- n_2d_figures + 1
    }
}
cat("Wrote", n_2d_figures, "2D figures to", figures_dir, "\n")
