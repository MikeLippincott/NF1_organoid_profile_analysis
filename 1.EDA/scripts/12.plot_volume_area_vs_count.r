list_of_packages <- c("ggplot2", "dplyr", "arrow", "RColorBrewer", "ggrastr")
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

# Biological question: as organoids are exposed to more cells (a proxy for
# treatment effect on growth), does organoid/cell volume change? Restricted
# to 3D: the 2D projection-method variants of this analysis compared
# representation choices, not biology, and were dropped.
#
# Organoid_NoChannel_AreaSizeShape_Volume and Cell_NoChannel_AreaSizeShape_Volume
# are z-scored PER PATIENT upstream (each patient's population is normalized
# against its own mean/SD), so patients are NOT on a common scale: e.g.
# mean/SD of organoid volume range from (0.49, 2.15) to (3.67, 12.0) across
# our 12 patients. We therefore never pool/overlay this column across
# patients on one shared axis -- only ever facet by patient, where each
# panel is valid on its own terms.

# --- FOV-normalized total cell count per patient x treatment x dose, from the
# canonical (non-sammed/non-nucleocentric) 3D single-cell profile type,
# matching sc_norm.parquet. There is no true per-organoid cell count in the
# current pipeline output, so this treatment-level count is broadcast onto
# every organoid below. ---
raw_counts <- read_parquet(file.path(results_dir, "cell_counts", "cell_counts.parquet"))

counts_3d <- raw_counts %>%
    filter(Metadata_profile_type == "sc_norm_norm_profile_3D") %>%
    group_by(Metadata_Biology_PatientTumor, Metadata_Experiment_Treatment, Metadata_Experiment_Dose) %>%
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
    select(Metadata_patient_tumor, Metadata_treatment, Metadata_dose, total_cell_count_norm)

# --- per-organoid volume and its own per-organoid cell count (both live on
# the same row in organoid_norm.parquet -- no need to broadcast a
# patient x treatment x dose mean onto every organoid; the real, unaggregated
# per-organoid count gives each point its own x value instead of collapsing
# a whole treatment group onto one repeated number) ---
patients_3d <- setdiff(list.dirs(file.path(root_dir, "data", "profiles_3D"), recursive = FALSE, full.names = FALSE), c("all_patients", "NF0037_T1_CQ1"))
vol_rows <- list()
for (patient in patients_3d) {
    f <- file.path(root_dir, "data", "profiles_3D", patient, "5.normalized_profiles", "organoid_norm.parquet")
    if (!file.exists(f)) next
    df <- read_parquet(f, col_select = c(
        "Metadata_Experiment_Treatment", "Metadata_Experiment_Dose",
        "Metadata_Object_OrganoidSingleCellCount", "Organoid_NoChannel_AreaSizeShape_Volume"
    ))
    df$Metadata_patient_tumor <- patient
    vol_rows[[patient]] <- df
}
vol_df <- bind_rows(vol_rows)
colnames(vol_df)[colnames(vol_df) == "Metadata_Experiment_Treatment"] <- "Metadata_treatment"
colnames(vol_df)[colnames(vol_df) == "Metadata_Experiment_Dose"] <- "Metadata_dose"

joined_3d <- vol_df
joined_3d$Metadata_treatment <- factor(joined_3d$Metadata_treatment,
                                          levels = intersect(custom_treatment_order, unique(joined_3d$Metadata_treatment)))

# --- per-organoid scatter, one PDF page per patient (the only axis this
# z-score is valid on), faceted by Treatment within each page instead of
# color-coded (compound colors are hard to distinguish across ~20
# treatments). Y axis is shared across a patient's treatment facets (all
# z-scored against that one patient's mean/SD, so they're comparable);
# x axis is free per facet since cell counts vary a lot by treatment. ---
make_patient_pages <- function(df, x_col, y_col, title_prefix, x_lab, y_lab, point_size, point_alpha) {
    lapply(sort(unique(df$Metadata_patient_tumor)), function(patient) {
        sub_df <- df %>% filter(Metadata_patient_tumor == patient)
        (
            ggplot(sub_df, aes(x = .data[[x_col]], y = .data[[y_col]]))
            + rasterise(geom_point(size = point_size, alpha = point_alpha, color = "steelblue"), dpi = 300)
            + facet_wrap(~Metadata_treatment, scales = "free_x")
            + labs(
                title = paste0(title_prefix, " -- patient ", patient),
                x = x_lab, y = y_lab
            )
            + theme_manuscript()
        )
    })
}

p_facet_patient_pages <- make_patient_pages(
    joined_3d, "Metadata_Object_OrganoidSingleCellCount", "Organoid_NoChannel_AreaSizeShape_Volume",
    "3D: organoid volume vs. cells per organoid, by treatment",
    "Cells per organoid", "Organoid volume (z-scored within patient)",
    point_size = 1, point_alpha = 0.5
)

# --- single-cell volume vs. the same FOV-normalized count, one PDF page per
# patient, faceted by Treatment within each page ---
sc_vol_rows <- list()
for (patient in patients_3d) {
    f <- file.path(root_dir, "data", "profiles_3D", patient, "5.normalized_profiles", "sc_norm.parquet")
    if (!file.exists(f)) next
    df <- read_parquet(f, col_select = c("Metadata_Experiment_Treatment", "Metadata_Experiment_Dose", "Cell_NoChannel_AreaSizeShape_Volume"))
    df$Metadata_patient_tumor <- patient
    sc_vol_rows[[patient]] <- df
}
sc_vol_df <- bind_rows(sc_vol_rows)
colnames(sc_vol_df)[colnames(sc_vol_df) == "Metadata_Experiment_Treatment"] <- "Metadata_treatment"
colnames(sc_vol_df)[colnames(sc_vol_df) == "Metadata_Experiment_Dose"] <- "Metadata_dose"

joined_sc <- sc_vol_df %>%
    inner_join(counts_3d, by = c("Metadata_patient_tumor", "Metadata_treatment", "Metadata_dose"))
joined_sc$Metadata_treatment <- factor(joined_sc$Metadata_treatment,
                                          levels = intersect(custom_treatment_order, unique(joined_sc$Metadata_treatment)))

p_sc_facet_patient_pages <- make_patient_pages(
    joined_sc, "total_cell_count_norm", "Cell_NoChannel_AreaSizeShape_Volume",
    "3D: single-cell volume vs. total cell count (FOV-normalized), by treatment",
    "Total cells per treatment (FOV-normalized)", "Cell volume (z-scored within patient)",
    point_size = 0.5, point_alpha = 0.4
)

pdf(file.path(figures_dir, "3D_volume_vs_count.pdf"), width = 11, height = 8.5, onefile = TRUE)
for (p in p_facet_patient_pages) print(p)
for (p in p_sc_facet_patient_pages) print(p)
dev.off()
