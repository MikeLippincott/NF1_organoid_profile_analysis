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
figures_dir <- file.path(root_dir, "4.analysis", "figures", "volume_area_vs_count")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

cell_counts <- read_parquet(file.path(results_dir, "cell_counts.parquet"))

plot_theme <- theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 9)
)

# --- 3D: join organoid volume + Metadata_Object_ObjectID to cell_counts (organoid_id) ---
patients_3d <- setdiff(list.dirs(file.path(root_dir, "data", "profiles_3D"), recursive = FALSE, full.names = FALSE), "all_patients")
vol_rows <- list()
for (patient in patients_3d) {
    f <- file.path(root_dir, "data", "profiles_3D", patient, "5.normalized_profiles", "organoid_norm.parquet")
    if (!file.exists(f)) next
    df <- read_parquet(f, col_select = c("Metadata_Experiment_Well", "Metadata_Object_ObjectID", "Organoid_NoChannel_AreaSizeShape_Volume"))
    df$Metadata_patient_tumor <- patient
    vol_rows[[patient]] <- df
}
vol_df <- bind_rows(vol_rows)

counts_3d <- cell_counts %>% filter(modality == "3D") %>%
    select(Metadata_patient_tumor, organoid_id, Metadata_treatment, cell_count)

joined_3d <- vol_df %>%
    inner_join(counts_3d, by = c("Metadata_patient_tumor", "Metadata_Object_ObjectID" = "organoid_id"))
joined_3d$Metadata_treatment <- factor(joined_3d$Metadata_treatment,
                                          levels = intersect(custom_treatment_order, unique(joined_3d$Metadata_treatment)))

p_3d_patient <- ggplot(joined_3d, aes(x = cell_count, y = Organoid_NoChannel_AreaSizeShape_Volume, color = Metadata_patient_tumor)) +
    geom_point(size = 1, alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
    scale_color_manual(values = setNames(tab20_palette_for_patients[seq_along(unique(joined_3d$Metadata_patient_tumor))],
                                           unique(joined_3d$Metadata_patient_tumor))) +
    labs(title = "3D: organoid volume vs. cell count, colored by patient",
         x = "Cells per organoid", y = "Organoid volume (z-scored)", color = "Patient") +
    plot_theme
ggsave(filename = file.path(figures_dir, "3D_volume_vs_count_by_patient.png"),
       plot = p_3d_patient, width = 10, height = 7, dpi = 600, units = "in")

p_3d_treatment <- ggplot(joined_3d, aes(x = cell_count, y = Organoid_NoChannel_AreaSizeShape_Volume, color = Metadata_treatment)) +
    geom_point(size = 1, alpha = 0.4) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
    scale_color_manual(values = custom_treatment_palette, na.value = "grey70") +
    labs(title = "3D pooled (all patients): organoid volume vs. cell count, colored by treatment",
         x = "Cells per organoid", y = "Organoid volume (z-scored)", color = "Treatment") +
    plot_theme
ggsave(filename = file.path(figures_dir, "3D_volume_vs_count_pooled_by_treatment.png"),
       plot = p_3d_treatment, width = 10, height = 7, dpi = 600, units = "in")

# --- 2D: area vs. well-level average cell count (no true per-organoid id link - see script 4 deviation note) ---
patients_2d <- setdiff(list.dirs(file.path(root_dir, "data", "profiles_2D"), recursive = FALSE, full.names = FALSE), "all_patients")
projection_prefix <- c(max_projection = "max_projected", middle_slice = "middle_slice", middle_n_slice = "middle_n_slice")
area_rows <- list()
for (proj in names(projection_prefix)) {
    prefix <- projection_prefix[[proj]]
    for (patient in patients_2d) {
        f <- file.path(root_dir, "data", "profiles_2D", patient, "5.normalized", paste0(prefix, "_organoid.parquet"))
        if (!file.exists(f)) next
        df <- read_parquet(f, col_select = c("Metadata_Well", "Metadata_treatment", "Metadata_Organoid_Number_Object_Number", "Organoid_AreaShape_Area"))
        df$Metadata_patient_tumor <- patient
        df$projection <- proj
        area_rows[[paste(proj, patient)]] <- df
    }
}
area_df <- bind_rows(area_rows)

counts_2d <- cell_counts %>% filter(modality == "2D") %>%
    select(Metadata_patient_tumor, organoid_id, projection, cell_count) %>%
    rename(Metadata_Organoid_Number_Object_Number = organoid_id)

joined_2d <- area_df %>%
    inner_join(counts_2d, by = c("Metadata_patient_tumor", "Metadata_Organoid_Number_Object_Number", "projection"))
joined_2d$Metadata_treatment <- factor(joined_2d$Metadata_treatment,
                                          levels = intersect(custom_treatment_order, unique(joined_2d$Metadata_treatment)))

p_2d_pooled <- ggplot(joined_2d, aes(x = cell_count, y = Organoid_AreaShape_Area)) +
    geom_bin2d(bins = 60) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.6) +
    scale_fill_viridis_c(trans = "log10") +
    facet_wrap(~projection) +
    labs(title = "2D pooled (all patients): organoid area vs. well-level avg. cell count, by projection method",
         x = "Mean cells per organoid (well-level average)", y = "Organoid area (z-scored)", fill = "Count (log10)") +
    plot_theme + theme(legend.position = "right")
ggsave(filename = file.path(figures_dir, "2D_area_vs_count_pooled_by_projection.png"),
       plot = p_2d_pooled, width = 14, height = 6, dpi = 600, units = "in")

cat("Wrote 3 figures to", figures_dir, "\n")
