list_of_packages <- c("ggplot2", "dplyr", "tidyr", "patchwork", "arrow", "ggrepel", "ggrastr", "scales", "ComplexHeatmap", "circlize", "grid", "RColorBrewer")
for (package in list_of_packages) {
    suppressPackageStartupMessages(
        suppressWarnings(
            library(package, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
        )
    )
}

find_git_root <- function() {
    current_path <- getwd()
    while (!dir.exists(file.path(current_path, ".git"))) {
        parent_path <- dirname(current_path)
        if (parent_path == current_path) stop("No Git root directory found.")
        current_path <- parent_path
    }
    current_path
}
root_dir <- find_git_root()
source(file.path(root_dir, "utils/r_plot_themes.r"))
source(file.path(root_dir, "utils/r_plot_funcs.r"))

# run tag of the tables made by 7.calculate_variate_class_upsets_and_clustermap (profile + variate groups)
cli_args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(cli_args) >= 1) cli_args[1] else "sc_T-O-C-N-M-X-Y-Z-D"

results_path <- file.path(root_dir, "4.linear_modeling/results/variate_class_plots", tag)
figures_path <- file.path(root_dir, "4.linear_modeling/figures/variate_class_plots", tag)
dir.create(file.path(figures_path, "upset"), recursive = TRUE, showWarnings = FALSE)
read_result <- function(name) arrow::read_parquet(file.path(results_path, name))

membership <- read_result("membership.parquet")
patients <- sort(unique(membership$patient))
treatments <- sort(unique(membership$treatment))
tag_pdf <- file.path(figures_path, paste0(tag, ".pdf"))
print(tag)

# Task A: UpSet plots, sets = patients. One pdf per class, one page per treatment.
plot_patient_upset <- function(combos) {
    n <- nrow(combos)
    n_pat <- length(patients)
    combos$x <- seq_len(n)
    use_log <- max(combos$n_features) / max(min(combos$n_features), 1) > 30
    p_bar <- (
        ggplot(combos, aes(x = x))
        + geom_rect(aes(xmin = x - 0.35, xmax = x + 0.35, ymin = if (use_log) 0.8 else 0, ymax = n_features), fill = "#333333")
        + geom_text(aes(y = n_features, label = n_features), vjust = -0.3, size = 3)
        + labs(
            x = NULL,
            y = paste0("features with exactly this\npatient combination", if (use_log) " (log)" else "")
        )
        + theme_classic(base_size = 9)
        + theme(
            axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()
        )
    )
    p_bar <- if (use_log) {
        p_bar + scale_y_log10(expand = c(0, 0)) + coord_cartesian(xlim = c(0.4, n + 0.6), ylim = c(0.8, max(combos$n_features) * 3))
    } else {
        p_bar + scale_y_continuous(expand = c(0, 0)) + coord_cartesian(xlim = c(0.4, n + 0.6), ylim = c(0, max(combos$n_features) * 1.18))
    }
    dots <- do.call(rbind, lapply(seq_len(n), function(i) {
        data.frame(
            x = i, row = seq_len(n_pat), patient = patients,
            on = patients %in% strsplit(combos$patients[i], ";", fixed = TRUE)[[1]]
        )
    }))
    stripes <- data.frame(row = seq_len(n_pat), shade = seq_len(n_pat) %% 2 == 1)
    lines_df <- dots |> filter(on) |> group_by(x) |> filter(n() > 1) |> ungroup()
    p_mat <- (
        ggplot()
        + geom_rect(
            data = filter(stripes, shade),
            aes(xmin = 0.4, xmax = n + 0.6, ymin = row - 0.5, ymax = row + 0.5), fill = "#f0f0f0"
        )
        + geom_point(data = dots, aes(x, row), colour = "#d8d8d8", size = 3)
        + geom_line(data = lines_df, aes(x, row, group = x), linewidth = 0.8)
        + geom_point(data = filter(dots, on), aes(x, row), colour = "black", size = 3)
        + scale_y_reverse(breaks = seq_len(n_pat), labels = patients, limits = c(n_pat + 0.5, 0.5), expand = c(0, 0))
        + scale_x_continuous(limits = c(0.4, n + 0.6), expand = c(0, 0))
        + labs(x = NULL, y = NULL)
        + theme_void(base_size = 9)
        + theme(axis.text.y = element_text(size = 9, hjust = 1))
    )
    (p_bar / p_mat) + plot_layout(heights = c(2.4, n_pat * 0.42))
}

# a class's pdf is only rebuilt when it is missing
upset_summary <- read_result("upset_top_combinations.parquet")
for (cls in sort(unique(upset_summary$class))) {
    pdf_path <- file.path(figures_path, "upset", paste0(cls, ".pdf"))
    if (file.exists(pdf_path)) next
    cls_df <- upset_summary |> filter(class == cls)
    plots <- list()
    widths <- c()
    heights <- c()
    for (trt in sort(unique(cls_df$treatment))) {
        combos <- cls_df |> filter(treatment == trt) |> arrange(rank)
        plots[[length(plots) + 1]] <- plot_patient_upset(combos)
        widths <- c(widths, max(7.5, 0.75 * nrow(combos) + 3.5))
        heights <- c(heights, 3.2 + 0.36 * length(patients))
    }
    save_plots_pdf(plots, pdf_path, width = widths, height = heights)
}
cat(length(unique(upset_summary$class)), "class upset pdfs in", file.path(figures_path, "upset"), "\n")

# every remaining figure is a page of the run's single pdf
pages <- list()
page_w <- c()
page_h <- c()
add_page <- function(plot, width, height) {
    pages[[length(pages) + 1]] <<- plot
    page_w <<- c(page_w, width)
    page_h <<- c(page_h, height)
}

# Task B: clustermap of all patient x treatment class-count profiles. log1p, column min-max scaling; rows and
# columns are hierarchically clustered here (Ward.D2 on the scaled values). Row annotations use the manuscript
# theme palettes; the bar at the bottom shows which variates are significant in each class (column).
patient_palette <- setNames(tab20_palette_for_patients[seq_along(tumor_type_lookup)], names(tumor_type_lookup))
dose_label_palette <- c("1uM" = dose_palette[["1"]], "10uM" = dose_palette[["10"]], "10nM" = "#8da0cb")
# class letters -> model variates (as coded by 7.calculate_variate_class_upsets_and_clustermap)
variate_letters <- c(
    T = "treatment", O = "organoid_count", C = "cell_per_organoid_count", N = "cell_count",
    M = "manhattan_distance_from_center", X = "cell_x_position", Y = "cell_y_position",
    Z = "cell_z_position", D = "cell_z_depth"
)

plot_clustermap <- function(matrix_name, order_prefix, label) {
    mat <- read_result(matrix_name) |>
        pivot_longer(-c(patient, treatment), names_to = "class", values_to = "value") |>
        mutate(key = paste(patient, treatment, sep = "|"), value = log1p(value))
    col_order <- read_result(paste0(order_prefix, "_column_order.parquet"))$column_order
    row_order <- read_result(paste0(order_prefix, "_row_order.parquet"))$row_order
    mat <- mat |>
        filter(class %in% col_order, key %in% row_order) |>
        group_by(class) |>
        mutate(scaled = (value - min(value)) / (max(value) - min(value))) |>
        ungroup()
    scaled <- long_to_matrix(mat, "key", "class", "scaled", row_levels = row_order, col_levels = col_order)
    scaled[is.na(scaled)] <- 0  # constant columns have no range to scale

    patient_of_row <- sub("\\|.*$", "", row_order)
    treatment_of_row <- sub("^[^|]*\\|", "", row_order)
    left_annotation <- rowAnnotation(
        Patient = patient_of_row, `Tumor type` = tumor_type_lookup[patient_of_row],
        Treatment = sub("_.*$", "", treatment_of_row), Dose = sub("^.*_", "", treatment_of_row),
        col = list(Patient = patient_palette, `Tumor type` = tumor_type_palette,
                   Treatment = custom_treatment_palette, Dose = dose_label_palette),
        annotation_name_gp = gpar(fontsize = 9), simple_anno_size = unit(4, "mm"),
        annotation_legend_param = list(Treatment = list(ncol = 2))
    )
    # presence / absence of each variate in each class (column)
    in_class <- sapply(names(variate_letters), function(v) ifelse(sapply(strsplit(col_order, "+", fixed = TRUE), function(x) v %in% x), "yes", "no"))
    colnames(in_class) <- paste0(names(variate_letters), " ", variate_letters)
    bottom_annotation <- HeatmapAnnotation(
        df = as.data.frame(in_class),
        col = setNames(rep(list(c(no = "white", yes = "#333333")), ncol(in_class)), colnames(in_class)),
        show_legend = c(TRUE, rep(FALSE, ncol(in_class) - 1)),
        annotation_legend_param = list(list(title = "variate significant\nin class", at = c("yes", "no"))),
        annotation_name_side = "right", annotation_name_gp = gpar(fontsize = 8),
        simple_anno_size = unit(3, "mm")
    )
    simple_heatmap(
        scaled, label, heat_col_fun(c(0, 1), palette = "Viridis", rev = FALSE),
        cluster_rows = TRUE, cluster_columns = TRUE,
        clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
        left_annotation = left_annotation, bottom_annotation = bottom_annotation,
        show_row_names = FALSE, show_column_names = FALSE, use_raster = TRUE
    )
}
add_page(
    heatmap_grid_page(list(plot_clustermap("class_count_matrix.parquet", "clustermap_counts", "column-scaled\nlog1p(count)"))),
    10, 12
)
add_page(
    heatmap_grid_page(list(plot_clustermap("class_fraction_matrix.parquet", "clustermap_fractions", "column-scaled\nlog1p(fraction)"))),
    10, 12
)

# Task C: which features are treatment-only (class T)?
recurrence <- read_result("treatment_only_feature_recurrence.parquet") |>
    mutate(channel = ifelse(is.na(channel), "none", channel))
# the top features are the ones 7.calculate_variate_class_upsets_and_clustermap kept (its --top_features)
per_treatment_wide <- read_result("treatment_only_feature_by_treatment.parquet")
top <- recurrence |> filter(feature %in% per_treatment_wide$feature)
top_features <- top$feature
top$feature <- factor(top$feature, levels = rev(top$feature))
p_bar <- (
    ggplot(top, aes(x = n_patient_treatments, y = feature, fill = feature_type))
    + geom_col()
    + labs(
        x = "# (patient, treatment) pairs where feature is treatment-only", y = NULL,
        fill = "feature type"
    )
    + theme_classic(base_size = 10)
    + theme(axis.text.y = element_text(size = 7))
)
add_page(p_bar, 12, max(6, 0.28 * nrow(top)))

grid_df <- recurrence |> count(channel, feature_type, name = "n_features")
grid_mat <- long_to_matrix(grid_df, "channel", "feature_type", "n_features")
add_page(
    heatmap_grid_page(list(simple_heatmap(
        grid_mat, "# distinct\nfeatures", heat_col_fun(grid_mat), cell_fmt = "%d", cell_size = 10,
        column_names_rot = 45
    ))),
    9, 7
)

per_treatment <- per_treatment_wide |>
    pivot_longer(-feature, names_to = "treatment", values_to = "n_patients")
per_treatment_mat <- long_to_matrix(per_treatment, "feature", "treatment", "n_patients", row_levels = top_features)
add_page(
    heatmap_grid_page(list(simple_heatmap(
        per_treatment_mat, "# patients\n(feature is\ntreatment-only)", heat_col_fun(per_treatment_mat),
        base_size = 9
    ))),
    max(14, 0.4 * length(unique(per_treatment$treatment)) + 9), max(6, 0.28 * nrow(top))
)

# Task D: distribution of the top shared treatment-only features
shared <- read_result("top_shared_features_coefficients.parquet")
top_shared <- intersect(recurrence$feature, unique(as.character(shared$feature)))
shared <- shared |> mutate(feature = factor(feature, levels = rev(top_shared)))
sig_pal <- c("TRUE" = "#d62728", "FALSE" = "#7f7f7f")
p_box <- (
    ggplot(shared, aes(x = coefficient, y = feature))
    + geom_boxplot(fill = "lightgrey", outlier.shape = NA)
    + geom_jitter(aes(colour = treatment_only), height = 0.2, size = 0.6, alpha = 0.5)
    + geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4)
    + scale_colour_manual(values = sig_pal, name = "treatment-only significant")
    + labs(x = "treatment coefficient (all patient x treatment models)", y = NULL)
    + theme_classic(base_size = 10)
)
frac <- shared |> group_by(feature) |> summarise(fraction = mean(treatment_only), .groups = "drop")
p_frac <- (
    ggplot(frac, aes(x = fraction, y = feature))
    + geom_col(fill = "#d62728")
    + labs(x = "fraction of models that are treatment-only significant", y = NULL)
    + theme_classic(base_size = 10)
    + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
)
add_page((p_box + p_frac) + plot_layout(widths = c(2.2, 1)), 16, 5.5)

space <- read_result("top_shared_features_space_distribution.parquet") |>
    pivot_longer(-c(value, part), names_to = "set", values_to = "fraction") |>
    mutate(part = factor(part, levels = c("compartment", "channel", "feature_type")))
p_space <- (
    ggplot(space, aes(x = value, y = fraction, fill = set))
    + geom_col(position = position_dodge(width = 0.8), width = 0.75)
    + facet_wrap(~part, scales = "free_x", nrow = 1)
    + scale_fill_manual(values = setNames(c("#bdbdbd", "#d62728"), sort(unique(space$set))))
    + labs(x = NULL, y = "fraction of features", fill = NULL)
    + theme_classic(base_size = 10)
    + theme(axis.text.x = element_text(angle = 45, hjust = 1))
)
add_page(p_space, 18, 4.5)

for (by in c("patient", "treatment")) {
    med <- shared |> group_by(feature, .data[[by]]) |> summarise(median_coef = median(coefficient), .groups = "drop")
    n_significant <- shared |> filter(treatment_only) |> count(feature, .data[[by]], name = "n_significant")
    med <- med |> left_join(n_significant, by = c("feature", by)) |> mutate(n_significant = ifelse(is.na(n_significant), 0L, n_significant))
    lim <- quantile(abs(med$median_coef), 0.98, na.rm = TRUE)
    med_mat <- long_to_matrix(med, "feature", by, "median_coef", row_levels = top_shared)
    significant_labels <- if (by == "patient") {
        matrix(as.character(long_to_matrix(med, "feature", by, "n_significant", row_levels = top_shared)), nrow = length(top_shared),
               dimnames = dimnames(med_mat))
    } else NULL
    ht <- simple_heatmap(
        med_mat, "median treatment\ncoefficient",
        circlize::colorRamp2(c(-lim, 0, lim), c("#3b4cc0", "white", "#b40426")),
        cell_labels = significant_labels, cell_size = 9, column_names_rot = 90
    )
    # extra width to fit the long feature-name row labels (see simple_heatmap's row_names_max_width) next to the heatmap body and legend
    row_label_width <- max(nchar(top_shared)) * 0.09 + 0.5
    base_width <- if (by == "patient") 8 else max(9, 0.4 * length(treatments) + 5)
    add_page(heatmap_grid_page(list(ht)), base_width + row_label_width, 5)
}

save_plots_pdf(pages, tag_pdf, width = page_w, height = page_h)
cat(length(pages), "pages ->", tag_pdf, "\n")
