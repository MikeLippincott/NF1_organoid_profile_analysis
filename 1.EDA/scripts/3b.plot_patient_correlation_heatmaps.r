# Load libraries
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    BiocManager::install("ComplexHeatmap")
}
suppressPackageStartupMessages({
    library(arrow)
    library(dplyr)
    library(tidyr)
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
})



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
    stop("No Git root directory found.")
}

# Find the Git root directory
root_dir <- find_git_root()

results_dir <- file.path(root_dir, "1.EDA", "results", "correlation")
figures_dir <- file.path(root_dir, "1.EDA", "figures", "correlation")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Function to generate a correlation heatmap
plot_correlation_heatmap <- function(
    corr_long_df,
    patient_id,
    profile_type,
    output_path
) {
    # Filter to patient
    patient_data <- corr_long_df

    cat(sprintf(
        "Patient %s (%s): %d rows\n",
        patient_id, profile_type, nrow(patient_data)
    ))

    # Pivot to wide matrix: rows = 2D features, columns = 3D features
    corr_matrix <- patient_data %>%
        select(feature_2d, feature_3d, pearson_r) %>%
        pivot_wider(
            names_from = feature_3d,
            values_from = pearson_r
        ) %>%
        tibble::column_to_rownames("feature_2d") %>%
        as.matrix()

    # Replace NaN with 0 for clustering
    corr_matrix[is.na(corr_matrix)] <- 0

    # Color scale from blue (-1) to white (0) to red (1)
    col_fun <- colorRamp2(
        c(-1, 0, 1),
        c("#2166AC", "white", "#B2182B")
    )

    # Create the heatmap
    ht <- Heatmap(
        corr_matrix,
        name = "Pearson\nCorrelation",
        col = col_fun,

        # Clustering
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2",

        # Hide individual feature labels
        show_row_names = FALSE,
        show_column_names = FALSE,

        # Axis labels
        row_title = "2D Features",
        row_title_side = "left",
        column_title = "3D Features",
        column_title_side = "bottom",

        # Dendrogram styling
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        row_dend_width = unit(15, "mm"),
        column_dend_height = unit(15, "mm"),

        # Legend
        heatmap_legend_param = list(
            title = "Pearson\nCorrelation",
            at = c(-1, -0.5, 0, 0.5, 1),
            labels = c("-1", "-0.5", "0", "0.5", "1"),
            legend_height = unit(4, "cm")
        ),
        use_raster = TRUE,
        raster_quality = 3
    )

    # Save to file
    png(
        output_path,
        width = 12,
        height = 10,
        units = "in",
        res = 300
    )

    # Draw with main title at top
    draw(
        ht,
        column_title = sprintf(
            "Pearson Correlation: 2D vs 3D Features — %s — %s",
            patient_id, profile_type
        ),
        column_title_gp = gpar(fontsize = 14, fontface = "bold")
    )

    dev.off()
    cat(sprintf("  Saved: %s\n", output_path))
}

#Organoid heatmaps
organoid_files <- list.files(results_dir, pattern = "_organoid_agg_correlation\\.parquet$", full.names = TRUE)

for (f in organoid_files) {
    # Extract patient ID from filename
    patient_id <- gsub("_organoid_agg_correlation\\.parquet$", "", basename(f))

    corr_data <- arrow::read_parquet(f)

    output_path <- file.path(figures_dir, sprintf("%s_organoid_heatmap.png", patient_id))

    plot_correlation_heatmap(
        corr_long_df = corr_data,
        patient_id = patient_id,
        profile_type = "Organoid",
        output_path = output_path
    )
}

#SC heatmaps
sc_files <- list.files(results_dir, pattern = "_sc_agg_correlation\\.parquet$", full.names = TRUE)

for (f in sc_files) {
    patient_id <- gsub("_sc_agg_correlation\\.parquet$", "", basename(f))

    corr_data <- arrow::read_parquet(f)

    output_path <- file.path(figures_dir, sprintf("%s_sc_heatmap.png", patient_id))

    plot_correlation_heatmap(
        corr_long_df = corr_data,
        patient_id = patient_id,
        profile_type = "Single Cell",
        output_path = output_path
    )
}
