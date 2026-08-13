packages <- c("ggplot2", "dplyr", "arrow", "ComplexHeatmap", "circlize", "scales")
for (pkg in packages) {
    suppressPackageStartupMessages(
        suppressWarnings(
            library(pkg, character.only = TRUE)
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

# matplotlib's tab20 categorical palette
tab20_palette <- c(
    "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a",
    "#d62728", "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
    "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d",
    "#17becf", "#9edae5"
)

correlation_dir <- file.path(root_dir, "1.EDA", "results", "correlation")
figures_dir <- file.path(root_dir, "1.EDA", "figures", "correlation_heatmaps")
if (!dir.exists(figures_dir)) {
    dir.create(figures_dir, recursive = TRUE)
}

correlation_files <- list.files(correlation_dir, full.names = TRUE)
length(correlation_files)

# Each correlation file has its own metadata schema (2D files use
# e.g. Metadata_treatment, 3D files use Metadata_Experiment_Treatment), so
# find columns by pattern instead of hard-coding a name.
plot_correlation_heatmap <- function(file_path, figures_dir) {
    df <- arrow::read_parquet(file_path)

    metadata_cols <- grep("^Metadata_", colnames(df), value = TRUE)
    sample_cols <- grep("^Sample_", colnames(df), value = TRUE)

    mat <- as.matrix(df[, sample_cols])
    dimnames(mat) <- NULL

    treatment_col <- metadata_cols[grepl("treatment", metadata_cols, ignore.case = TRUE)]
    dose_col <- metadata_cols[grepl("dose", metadata_cols, ignore.case = TRUE) &
        !grepl("unit", metadata_cols, ignore.case = TRUE)]
    patient_col <- metadata_cols[grepl("patient|tumor", metadata_cols, ignore.case = TRUE)]

    annotations <- list()

    if (length(treatment_col) >= 1 && length(dose_col) >= 1) {
        # custom_treatment_palette (from utils/r_plot_themes.r) is keyed by
        # "Treatment_Dose", matching how these profiles encode treatment.
        treatment_dose_values <- paste0(
            as.character(df[[treatment_col[1]]]), "_", as.character(df[[dose_col[1]]])
        )
        annotations$Treatment <- treatment_dose_values
    }

    if (length(patient_col) >= 1) {
        patient_values <- as.character(df[[patient_col[1]]])
        unique_patients <- sort(unique(patient_values))
        patient_colors <- setNames(
            tab20_palette[seq_along(unique_patients)],
            unique_patients
        )
        annotations$Patient <- patient_values
    }

    top_annotation <- NULL
    if (length(annotations) > 0) {
        anno_cols <- list()
        if (!is.null(annotations$Treatment)) {
            anno_cols$Treatment <- custom_treatment_palette
        }
        if (!is.null(annotations$Patient)) {
            anno_cols$Patient <- patient_colors
        }
        top_annotation <- do.call(HeatmapAnnotation, c(
            annotations,
            list(
                col = anno_cols,
                show_legend = c(Treatment = TRUE, Patient = TRUE)[names(annotations)],
                annotation_legend_param = list(
                    Treatment = list(title = "Treatment", ncol = 1, labels_gp = gpar(fontsize = 8))
                ),
                annotation_name_side = "left",
                annotation_name_gp = gpar(fontsize = 10)
            )
        ))
    }

    col_fun <- colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))

    # A handful of samples have zero-variance feature vectors, which makes
    # their Pearson correlation undefined (NA). hclust can't cluster on NA
    # distances, so skip dendrograms for any matrix that contains them.
    can_cluster <- !anyNA(mat)

    heatmap_plot <- Heatmap(
        mat,
        name = "Corr",
        col = col_fun,
        na_col = "grey80",
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_row_dend = FALSE,
        show_column_dend = can_cluster,
        cluster_rows = can_cluster,
        cluster_columns = can_cluster,
        top_annotation = top_annotation,
        column_title = tools::file_path_sans_ext(basename(file_path)),
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        heatmap_legend_param = list(title = "Pearson\ncorrelation")
    )

    output_path <- file.path(
        figures_dir,
        paste0(tools::file_path_sans_ext(basename(file_path)), "_heatmap.png")
    )

    plot_width <- max(8, min(24, nrow(mat) / 40 + 3))
    plot_height <- max(6, min(20, nrow(mat) / 40))
    png(output_path, width = plot_width, height = plot_height, units = "in", res = 300)
    draw(heatmap_plot, merge_legend = TRUE)
    dev.off()

    output_path
}

saved_paths <- character(length(correlation_files))
for (i in seq_along(correlation_files)) {
    saved_paths[i] <- plot_correlation_heatmap(correlation_files[i], figures_dir)
    cat("Saved:", saved_paths[i], "\n")
}
