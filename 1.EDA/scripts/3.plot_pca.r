list_of_packages <- c("ggplot2", "dplyr", "tidyr", "circlize", "RColorBrewer")
for (package in list_of_packages) {
    suppressPackageStartupMessages(
        suppressWarnings(
            library(
                package,
                character.only = TRUE,
                quietly = TRUE,
                warn.conflicts = FALSE
            )
        )
    )
}

# Get the current working directory and find Git root
find_git_root <- function() {
    # Get current working directory
    cwd <- getwd()

    # Check if current directory has .git
    if (dir.exists(file.path(cwd, ".git"))) {
        return(cwd)
    }

    # If not, search parent directories
    current_path <- cwd
    while (dirname(current_path) != current_path) {  # While not at root
        parent_path <- dirname(current_path)
        if (dir.exists(file.path(parent_path, ".git"))) {
            return(parent_path)
        }
        current_path <- parent_path
    }

    # If no Git root found, stop with error
    stop("No Git root directory found.")
}

# Find the Git root directory
root_dir <- find_git_root()
cat("Git root directory:", root_dir, "\n")
source(file.path(root_dir, "utils", "r_plot_themes.r"))

figures_path <- file.path(root_dir,"1.EDA/figures/pca")
if (!dir.exists(figures_path)) {
  dir.create(figures_path, recursive = TRUE)
}

pca_theme <- theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size = 12)
    )

width <- 8
height <- 8
options(repr.plot.width = width, repr.plot.height = height)
make_pca_plot <- function(df, explained_variance_df, title, treatment_col, save_path, facet_col = NULL, width = 10, height = 5) {
    # Extract explained variance for PC0 and PC1
    # Assumes explained_variance_df has columns: PC (e.g. "PC0", "PC1") and explained_variance (proportion, 0-1)
    var_pc0 <- explained_variance_df$PC0_explained_variance
    var_pc1 <- explained_variance_df$PC1_explained_variance

    x_label <- sprintf("PC0 (%.1f%%)", var_pc0 * 100)
    y_label <- sprintf("PC1 (%.1f%%)", var_pc1 * 100)

    p <- (
        ggplot(df, aes(x = PC0, y = PC1, color = .data[[treatment_col]]))
        + geom_point(alpha = 0.5, size = 1)
        + scale_color_manual(values = custom_treatment_palette)
        + labs(title = title, x = x_label, y = y_label)
        + theme_bw()
        + pca_theme
        + guides(
            color = guide_legend(
                title = "Treatment",
                override.aes = list(alpha = 1, size = 5),
                ncol = 1
            )
        )
    )
    if (!is.null(facet_col)) {
        p <- p + facet_wrap(as.formula(paste0("~", facet_col)), nrow = 3)
    }
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggsave(p, file = save_path, width = width, height = height, dpi = 300)
    return(p)
}

pca_results_dir <- file.path(root_dir, "1.EDA/results/pca")

slice_specs <- list(
    list(dir_name = "max_projection", file_prefix = "2D_max_projection", title_label = "2D MIP", save_prefix = "all_patients_2D_MIP"),
    list(dir_name = "middle_slice", file_prefix = "2D_middle_slice", title_label = "2D Middle Slice", save_prefix = "all_patients_2D_middle_slice"),
    list(dir_name = "middle_n_slice", file_prefix = "2D_middle_n_slice", title_label = "2D Middle N Slice", save_prefix = "all_patients_2D_middle_n_slice")
)

profile_specs_2d <- list(
    list(suffix = "fs_profiles", label = "FS", file_tag = "features", facet = TRUE),
    list(suffix = "agg_profiles", label = "Aggregate", file_tag = "agg_features", facet = FALSE),
    list(suffix = "consensus_profiles", label = "Consensus", file_tag = "consensus_features", facet = FALSE)
)

entity_specs <- list(
    list(entity = "sc", entity_label = "Single Cell"),
    list(entity = "organoid", entity_label = "Organoid")
)

for (slice in slice_specs) {
    for (entity in entity_specs) {
        for (profile in profile_specs_2d) {
            file_name <- paste0(slice$file_prefix, "_", entity$entity, "_", profile$suffix, "_embeddings.parquet")
            file_path <- file.path(pca_results_dir, file_name)
            explained_variance_file_name <- paste0(slice$file_prefix, "_", entity$entity, "_", profile$suffix, "_explained_variance.parquet")
            explained_variance_df <- arrow::read_parquet(file.path(pca_results_dir, explained_variance_file_name))
            if (!file.exists(file_path)) {
                cat("Missing file, skipping:", file_path, "\n")
                next
            }

            df <- arrow::read_parquet(file_path)
            title <- paste0("All patients: ", slice$title_label, " ", entity$entity_label, " ", profile$label, " Profiles")
            save_dir <- file.path(figures_path, "2D", slice$dir_name)

            save_path <- file.path(save_dir, paste0(slice$save_prefix, "_pca_", entity$entity, "_", profile$file_tag, ".png"))
            p <- make_pca_plot(df, explained_variance_df, title, "Metadata_treatment", save_path, width = width, height = height)
            print(p)

            if (profile$facet && "Metadata_patient" %in% colnames(df)) {
                save_path_faceted <- file.path(save_dir, paste0(slice$save_prefix, "_pca_", entity$entity, "_", profile$file_tag, "_facet_by_patient.png"))
                p_faceted <- make_pca_plot(df, explained_variance_df,title, "Metadata_treatment", save_path_faceted, facet_col = "Metadata_patient_tumor", width = width, height = height)
                print(p_faceted)
            }
        }
    }
}

# normalized profiles
normalized_profiles <- c(
    "organoid_norm",
    "sammed_organoid_norm",
    "sc_norm",
    "sammed_sc_norm",
    "nucleocentric_morphem_norm",
    "sammed_nucleocentric_norm"
)

profile_specs_3d <- list(
    list(file_prefix = "3D_1.feature_selected_profiles", suffix = "fs_profiles", label = "Feature Selected", file_tag = "fs_features", facet = TRUE),
    list(file_prefix = "3D_2.aggregated_profiles", suffix = "sc_agg_profiles", label = "Aggregated", file_tag = "agg_features", facet = FALSE),
    list(file_prefix = "3D_3.consensus_profiles", suffix = "sc_consensus_profiles", label = "Consensus", file_tag = "consensus_features", facet = FALSE)
)

for (norm_profile in normalized_profiles) {
    for (profile in profile_specs_3d) {
        file_name <- paste0(profile$file_prefix, "_", norm_profile, "_", profile$suffix, "_embeddings.parquet")
        file_path <- file.path(pca_results_dir, file_name)
        if (!file.exists(file_path)) {
            cat("Missing file, skipping:", file_path, "\n")
            next
        }



        df <- arrow::read_parquet(file_path)
        if (!"PC1" %in% colnames(df)) {
            warning("PC1 column not found in df — skipping plot for '", title, "'.")
            next
        }
        explained_variance_file_name <- paste0(profile$file_prefix, "_", norm_profile, "_", profile$suffix, "_explained_variance.parquet")
        explained_variance_df <- arrow::read_parquet(file.path(pca_results_dir, explained_variance_file_name))
        title <- paste0("All patients: 3D (", norm_profile, ") ", profile$label, " Profiles")
        save_dir <- file.path(figures_path, "3D", norm_profile)

        save_path <- file.path(save_dir, paste0("all_patients_3D_pca_", profile$file_tag, ".png"))
        p <- make_pca_plot(df, explained_variance_df,title, "Metadata_Experiment_Treatment", save_path, width = width, height = height)
        print(p)

        if (profile$facet && "Metadata_Biology_PatientTumor" %in% colnames(df)) {
            save_path_faceted <- file.path(save_dir, paste0("all_patients_3D_pca_", profile$file_tag, "_facet_by_patient.png"))
            p_faceted <- make_pca_plot(df, explained_variance_df,title, "Metadata_Experiment_Treatment", save_path_faceted, facet_col = "Metadata_Biology_PatientTumor", width = width, height = height)
            print(p_faceted)
        }
    }
}
