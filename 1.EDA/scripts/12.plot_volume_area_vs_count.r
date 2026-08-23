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
figures_dir <- file.path(root_dir, "1.EDA", "figures", "volume_area_vs_count")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

plot_theme <- (
    theme_bw()
    + theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 9)
    )
)

# --- FOV-normalized total cell count per patient x treatment x dose (x
# projection, for 2D), from the canonical (non-sammed/non-nucleocentric)
# single-cell profile type per modality, matching sc_norm.parquet / *_sc.parquet.
# There is no true per-organoid cell count in the current pipeline output, so
# this treatment-level count is broadcast onto every organoid below. ---
raw_counts <- read_parquet(file.path(results_dir, "cell_counts", "cell_counts.parquet"))

canonical_sc_profile_types <- data.frame(
    Metadata_profile_type = c(
        "sc_norm_norm_profile_3D",
        "sc_profiles_2D_max_projection",
        "sc_profiles_2D_middle_slice",
        "sc_profiles_2D_middle_n_slice"
    ),
    modality = c("3D", "2D", "2D", "2D"),
    projection = c("none", "max_projection", "middle_slice", "middle_n_slice")
)

total_cell_counts <- raw_counts %>%
    inner_join(canonical_sc_profile_types, by = "Metadata_profile_type") %>%
    group_by(Metadata_Biology_PatientTumor, Metadata_Experiment_Treatment, Metadata_Experiment_Dose, modality, projection) %>%
    summarise(
        total_cells = sum(Metadata_n_cells),
        total_fovs = sum(Metadata_n_fovs),
        .groups = "drop"
    ) %>%
    mutate(total_cell_count_norm = total_cells / total_fovs) %>%
    rename(
        Metadata_patient_tumor = Metadata_Biology_PatientTumor,
        Metadata_treatment = Metadata_Experiment_Treatment,
        Metadata_dose = Metadata_Experiment_Dose
    ) %>%
    select(Metadata_patient_tumor, Metadata_treatment, Metadata_dose, modality, projection, total_cell_count_norm)

# --- 3D: per-organoid volume vs. treatment-level FOV-normalized total cell count
# (broadcast to every organoid in that patient x treatment x dose) ---
patients_3d <- setdiff(list.dirs(file.path(root_dir, "data", "profiles_3D"), recursive = FALSE, full.names = FALSE), c("all_patients", "NF0037_T1_CQ1"))
vol_rows <- list()
for (patient in patients_3d) {
    f <- file.path(root_dir, "data", "profiles_3D", patient, "5.normalized_profiles", "organoid_norm.parquet")
    if (!file.exists(f)) next
    df <- read_parquet(f, col_select = c("Metadata_Experiment_Treatment", "Metadata_Experiment_Dose", "Organoid_NoChannel_AreaSizeShape_Volume"))
    df$Metadata_patient_tumor <- patient
    vol_rows[[patient]] <- df
}
vol_df <- bind_rows(vol_rows)
colnames(vol_df)[colnames(vol_df) == "Metadata_Experiment_Treatment"] <- "Metadata_treatment"
colnames(vol_df)[colnames(vol_df) == "Metadata_Experiment_Dose"] <- "Metadata_dose"

counts_3d <- total_cell_counts %>% filter(modality == "3D") %>%
    select(Metadata_patient_tumor, Metadata_treatment, Metadata_dose, total_cell_count_norm)

joined_3d <- vol_df %>%
    inner_join(counts_3d, by = c("Metadata_patient_tumor", "Metadata_treatment", "Metadata_dose"))
joined_3d$Metadata_treatment <- factor(joined_3d$Metadata_treatment,
                                          levels = intersect(custom_treatment_order, unique(joined_3d$Metadata_treatment)))

p_3d_patient <- (
    ggplot(joined_3d, aes(x = total_cell_count_norm, y = Organoid_NoChannel_AreaSizeShape_Volume, color = Metadata_patient_tumor))
    + geom_point(size = 1, alpha = 0.5)
    + geom_smooth(method = "lm", se = FALSE, linewidth = 0.6)
    + scale_color_manual(values = setNames(
        tab20_palette_for_patients[seq_along(unique(joined_3d$Metadata_patient_tumor))],
        unique(joined_3d$Metadata_patient_tumor)
    ))
    + labs(
        title = "3D: organoid volume vs. total cell count (FOV-normalized), colored by patient",
        x = "Total cells per treatment (FOV-normalized)", y = "Organoid volume (z-scored)", color = "Patient"
    )
    + plot_theme
)
ggsave(
    filename = file.path(figures_dir, "3D_volume_vs_count_by_patient.png"),
    plot = p_3d_patient, width = 10, height = 7, dpi = 600, units = "in"
)

p_3d_treatment <- (
    ggplot(joined_3d, aes(x = total_cell_count_norm, y = Organoid_NoChannel_AreaSizeShape_Volume, color = Metadata_treatment))
    + geom_point(size = 1, alpha = 0.4)
    + geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6)
    + scale_color_manual(values = custom_treatment_palette, na.value = "grey70")
    + labs(
        title = "3D pooled (all patients): organoid volume vs. total cell count (FOV-normalized), colored by treatment",
        x = "Total cells per treatment (FOV-normalized)", y = "Organoid volume (z-scored)", color = "Treatment"
    )
    + plot_theme
)
ggsave(
    filename = file.path(figures_dir, "3D_volume_vs_count_pooled_by_treatment.png"),
    plot = p_3d_treatment, width = 10, height = 7, dpi = 600, units = "in"
)

p_3d_pooled_density <- (
    ggplot(joined_3d, aes(x = total_cell_count_norm, y = Organoid_NoChannel_AreaSizeShape_Volume))
    + geom_bin2d(bins = 60)
    + geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.6)
    + scale_fill_viridis_c(trans = "log10")
    + labs(
        title = "3D pooled (all patients): organoid volume vs. total cell count (FOV-normalized)",
        x = "Total cells per treatment (FOV-normalized)", y = "Organoid volume (z-scored)", fill = "Count (log10)"
    )
    + plot_theme
    + theme(legend.position = "right")
)
ggsave(
    filename = file.path(figures_dir, "3D_volume_vs_count_pooled_density.png"),
    plot = p_3d_pooled_density, width = 10, height = 7, dpi = 600, units = "in"
)

# --- 3D: same relationship, faceted by patient, by treatment, and by both ---
p_3d_facet_patient <- (
    ggplot(joined_3d, aes(x = total_cell_count_norm, y = Organoid_NoChannel_AreaSizeShape_Volume, color = Metadata_treatment))
    + geom_point(size = 1, alpha = 0.5)
    + scale_color_manual(values = custom_treatment_palette, na.value = "grey70")
    + facet_wrap(~Metadata_patient_tumor, scales = "free")
    + labs(
        title = "3D: organoid volume vs. total cell count (FOV-normalized), by patient",
        x = "Total cells per treatment (FOV-normalized)", y = "Organoid volume (z-scored)", color = "Treatment"
    )
    + plot_theme
)
ggsave(
    filename = file.path(figures_dir, "3D_volume_vs_count_facet_by_patient.png"),
    plot = p_3d_facet_patient, width = 14, height = 10, dpi = 600, units = "in"
)

p_3d_facet_treatment <- (
    ggplot(joined_3d, aes(x = total_cell_count_norm, y = Organoid_NoChannel_AreaSizeShape_Volume, color = Metadata_patient_tumor))
    + geom_point(size = 1, alpha = 0.5)
    + scale_color_manual(values = setNames(
        tab20_palette_for_patients[seq_along(unique(joined_3d$Metadata_patient_tumor))],
        unique(joined_3d$Metadata_patient_tumor)
    ))
    + facet_wrap(~Metadata_treatment, scales = "free")
    + labs(
        title = "3D: organoid volume vs. total cell count (FOV-normalized), by treatment",
        x = "Total cells per treatment (FOV-normalized)", y = "Organoid volume (z-scored)", color = "Patient"
    )
    + plot_theme
    + theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1))
)
ggsave(
    filename = file.path(figures_dir, "3D_volume_vs_count_facet_by_treatment.png"),
    plot = p_3d_facet_treatment, width = 18, height = 14, dpi = 600, units = "in"
)

p_3d_facet_both <- (
    ggplot(joined_3d, aes(x = total_cell_count_norm, y = Organoid_NoChannel_AreaSizeShape_Volume))
    + geom_point(size = 0.8, alpha = 0.4, color = "#3D7DCC")
    + facet_grid(Metadata_treatment ~ Metadata_patient_tumor, scales = "free")
    + labs(
        title = "3D: organoid volume vs. total cell count (FOV-normalized), by patient and treatment",
        x = "Total cells per treatment (FOV-normalized)", y = "Organoid volume (z-scored)"
    )
    + plot_theme
    + theme(
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
        strip.text.y = element_text(size = 6, angle = 0)
    )
)
ggsave(
    filename = file.path(figures_dir, "3D_volume_vs_count_facet_by_patient_and_treatment.png"),
    plot = p_3d_facet_both, width = 16, height = 20, dpi = 600, units = "in"
)

# --- 2D: per-organoid area vs. treatment-level FOV-normalized total cell count
# (broadcast to every organoid in that patient x treatment x dose x projection),
# pooled density plot faceted by projection method ---
patients_2d <- setdiff(list.dirs(file.path(root_dir, "data", "profiles_2D"), recursive = FALSE, full.names = FALSE), c("all_patients", "NF0037_T1_CQ1"))
projection_prefix <- c(max_projection = "max_projected", middle_slice = "middle_slice", middle_n_slice = "middle_n_slice")
area_rows <- list()
for (proj in names(projection_prefix)) {
    prefix <- projection_prefix[[proj]]
    for (patient in patients_2d) {
        f <- file.path(root_dir, "data", "profiles_2D", patient, "5.normalized", paste0(prefix, "_organoid.parquet"))
        if (!file.exists(f)) next
        df <- read_parquet(f, col_select = c("Metadata_treatment", "Metadata_dose", "Organoid_AreaShape_Area"))
        df$Metadata_patient_tumor <- patient
        df$projection <- proj
        area_rows[[paste(proj, patient)]] <- df
    }
}
area_df <- bind_rows(area_rows)

counts_2d <- total_cell_counts %>% filter(modality == "2D") %>%
    select(Metadata_patient_tumor, Metadata_treatment, Metadata_dose, projection, total_cell_count_norm)

joined_2d <- area_df %>%
    inner_join(counts_2d, by = c("Metadata_patient_tumor", "Metadata_treatment", "Metadata_dose", "projection"))
joined_2d$Metadata_treatment <- factor(joined_2d$Metadata_treatment,
                                          levels = intersect(custom_treatment_order, unique(joined_2d$Metadata_treatment)))

p_2d_pooled <- (
    ggplot(joined_2d, aes(x = total_cell_count_norm, y = Organoid_AreaShape_Area))
    + geom_bin2d(bins = 60)
    + geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.6)
    + scale_fill_viridis_c(trans = "log10")
    + facet_wrap(~projection)
    + labs(
        title = "2D pooled (all patients): organoid area vs. total cell count (FOV-normalized), by projection method",
        x = "Total cells per treatment (FOV-normalized)", y = "Organoid area (z-scored)", fill = "Count (log10)"
    )
    + plot_theme
    + theme(legend.position = "right")
)
ggsave(
    filename = file.path(figures_dir, "2D_area_vs_count_pooled_by_projection.png"),
    plot = p_2d_pooled, width = 14, height = 6, dpi = 600, units = "in"
)

p_2d_patient <- (
    ggplot(joined_2d, aes(x = total_cell_count_norm, y = Organoid_AreaShape_Area, color = Metadata_patient_tumor))
    + geom_point(size = 1, alpha = 0.4)
    + geom_smooth(method = "lm", se = FALSE, linewidth = 0.6)
    + scale_color_manual(values = setNames(
        tab20_palette_for_patients[seq_along(unique(joined_2d$Metadata_patient_tumor))],
        unique(joined_2d$Metadata_patient_tumor)
    ))
    + facet_wrap(~projection)
    + labs(
        title = "2D: organoid area vs. total cell count (FOV-normalized), colored by patient, by projection method",
        x = "Total cells per treatment (FOV-normalized)", y = "Organoid area (z-scored)", color = "Patient"
    )
    + plot_theme
)
ggsave(
    filename = file.path(figures_dir, "2D_area_vs_count_by_patient.png"),
    plot = p_2d_patient, width = 16, height = 6, dpi = 600, units = "in"
)

p_2d_treatment <- (
    ggplot(joined_2d, aes(x = total_cell_count_norm, y = Organoid_AreaShape_Area, color = Metadata_treatment))
    + geom_point(size = 1, alpha = 0.3)
    + geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6)
    + scale_color_manual(values = custom_treatment_palette, na.value = "grey70")
    + facet_wrap(~projection)
    + labs(
        title = "2D pooled (all patients): organoid area vs. total cell count (FOV-normalized), colored by treatment, by projection method",
        x = "Total cells per treatment (FOV-normalized)", y = "Organoid area (z-scored)", color = "Treatment"
    )
    + plot_theme
)
ggsave(
    filename = file.path(figures_dir, "2D_area_vs_count_pooled_by_treatment.png"),
    plot = p_2d_treatment, width = 16, height = 6, dpi = 600, units = "in"
)

cat("Wrote 5 figures to", figures_dir, "\n")

# --- 2D: same relationship, faceted by patient, by treatment, and by both
# (all 2D methods included via the projection dimension) ---
p_2d_facet_patient <- (
    ggplot(joined_2d, aes(x = total_cell_count_norm, y = Organoid_AreaShape_Area, color = Metadata_treatment))
    + geom_point(size = 0.8, alpha = 0.4)
    + scale_color_manual(values = custom_treatment_palette, na.value = "grey70")
    + facet_grid(projection ~ Metadata_patient_tumor, scales = "free")
    + labs(
        title = "2D: organoid area vs. total cell count (FOV-normalized), by patient x projection method",
        x = "Total cells per treatment (FOV-normalized)", y = "Organoid area (z-scored)", color = "Treatment"
    )
    + plot_theme
    + theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1), strip.text = element_text(size = 7))
)
ggsave(
    filename = file.path(figures_dir, "2D_area_vs_count_facet_by_patient.png"),
    plot = p_2d_facet_patient, width = 20, height = 8, dpi = 600, units = "in"
)

p_2d_facet_treatment <- (
    ggplot(joined_2d, aes(x = total_cell_count_norm, y = Organoid_AreaShape_Area, color = Metadata_patient_tumor))
    + geom_point(size = 0.8, alpha = 0.4)
    + scale_color_manual(values = setNames(
        tab20_palette_for_patients[seq_along(unique(joined_2d$Metadata_patient_tumor))],
        unique(joined_2d$Metadata_patient_tumor)
    ))
    + facet_grid(projection ~ Metadata_treatment, scales = "free")
    + labs(
        title = "2D: organoid area vs. total cell count (FOV-normalized), by treatment x projection method",
        x = "Total cells per treatment (FOV-normalized)", y = "Organoid area (z-scored)", color = "Patient"
    )
    + plot_theme
    + theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1), strip.text = element_text(size = 7))
)
ggsave(
    filename = file.path(figures_dir, "2D_area_vs_count_facet_by_treatment.png"),
    plot = p_2d_facet_treatment, width = 22, height = 8, dpi = 600, units = "in"
)

n_facet_both_figures <- 0
for (proj in unique(joined_2d$projection)) {
    df_proj <- joined_2d %>% filter(projection == proj)
    p_both <- (
        ggplot(df_proj, aes(x = total_cell_count_norm, y = Organoid_AreaShape_Area))
        + geom_point(size = 0.6, alpha = 0.3, color = "#3D7DCC")
        + facet_grid(Metadata_treatment ~ Metadata_patient_tumor, scales = "free")
        + labs(
            title = paste0("2D (", proj, "): organoid area vs. total cell count (FOV-normalized), by patient and treatment"),
            x = "Total cells per treatment (FOV-normalized)", y = "Organoid area (z-scored)"
        )
        + plot_theme
        + theme(
            axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
            strip.text.y = element_text(size = 6, angle = 0)
        )
    )
    ggsave(
        filename = file.path(figures_dir, paste0("2D_", proj, "_area_vs_count_facet_by_patient_and_treatment.png")),
        plot = p_both, width = 16, height = 20, dpi = 600, units = "in"
    )
    n_facet_both_figures <- n_facet_both_figures + 1
}

cat("Wrote", 2 + n_facet_both_figures, "2D facet figures to", figures_dir, "\n")
