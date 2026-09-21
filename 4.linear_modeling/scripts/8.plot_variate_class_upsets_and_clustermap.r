list_of_packages <- c("ggplot2", "dplyr", "tidyr", "patchwork", "arrow", "ggrepel", "ggrastr", "scales")
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
source(file.path(root_dir, "utils/r_plot_funcs.r"))

# run tag of the tables made by 5.variate_class_upsets_and_clustermap (profile + variate groups)
cli_args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(cli_args) >= 1) cli_args[1] else "sc_T-O-C-N-M-X-Y-Z-D"
TOP_FEATURES <- 40  # features shown in the treatment-only plots
N_SHARED <- 10      # shared treatment-only features in task D

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
plot_patient_upset <- function(combos, title) {
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
            y = paste0("features with exactly this\npatient combination", if (use_log) " (log)" else ""),
            title = title
        )
        + theme_classic(base_size = 9)
        + theme(
            axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
            plot.title = element_text(size = 11)
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
        plots[[length(plots) + 1]] <- plot_patient_upset(combos, paste0(trt, " | class ", cls, " | sets = patients"))
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

# Task B: clustermap of all patient x treatment class-count profiles. log1p, column min-max scaling; the Ward
# row / column order comes from the tables (the dendrograms themselves are not drawn).
patient_pal <- setNames(hue_pal()(length(patients)), patients)
# stepped by 7 (coprime with 25) so neighbouring treatments get clearly different colours
treatment_pal <- setNames(hue_pal()(length(treatments))[((seq_along(treatments) - 1) * 7) %% length(treatments) + 1], treatments)

plot_clustermap <- function(matrix_name, order_prefix, label, title) {
    mat <- read_result(matrix_name) |>
        pivot_longer(-c(patient, treatment), names_to = "class", values_to = "value") |>
        mutate(key = paste(patient, treatment, sep = "|"), value = log1p(value))
    col_order <- read_result(paste0(order_prefix, "_column_order.parquet"))$column_order
    row_order <- read_result(paste0(order_prefix, "_row_order.parquet"))$row_order
    mat <- mat |>
        filter(class %in% col_order, key %in% row_order) |>
        group_by(class) |>
        mutate(scaled = (value - min(value)) / (max(value) - min(value))) |>
        ungroup() |>
        mutate(class = factor(class, levels = col_order), key = factor(key, levels = rev(row_order)))
    strips <- data.frame(key = factor(row_order, levels = rev(row_order))) |>
        mutate(patient = sub("\\|.*$", "", key), treatment = sub("^[^|]*\\|", "", key))
    strip_theme <- theme_void() + theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9))
    p_pat <- (
        ggplot(strips, aes(x = 1, y = key, fill = patient)) + geom_tile()
        + scale_fill_manual(values = patient_pal, name = "patient") + strip_theme
    )
    p_trt <- (
        ggplot(strips, aes(x = 1, y = key, fill = treatment)) + geom_tile()
        + scale_fill_manual(values = treatment_pal, name = "treatment") + strip_theme
    )
    p_heat <- (
        ggplot(mat, aes(x = class, y = key, fill = scaled))
        + geom_tile()
        + scale_fill_viridis_c(name = label, na.value = "grey90")
        + labs(x = NULL, y = NULL, title = title)
        + theme_minimal(base_size = 9)
        + theme(
            axis.text.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid = element_blank(), legend.position = "top"
        )
    )
    (p_pat + p_trt + p_heat) + plot_layout(widths = c(0.03, 0.03, 1))
}
add_page(
    plot_clustermap("class_count_matrix.parquet", "clustermap_counts", "column-scaled log1p(count)",
                    "Clustermap of patient x treatment class-count profiles"),
    10, 12
)
add_page(
    plot_clustermap("class_fraction_matrix.parquet", "clustermap_fractions", "column-scaled log1p(fraction)",
                    "Clustermap of patient x treatment class-fraction profiles"),
    10, 12
)

# Task C: which features are treatment-only (class T)?
recurrence <- read_result("treatment_only_feature_recurrence.parquet") |>
    mutate(channel = ifelse(is.na(channel), "none", channel))
top <- head(recurrence, TOP_FEATURES)
top$feature <- factor(top$feature, levels = rev(top$feature))
p_bar <- (
    ggplot(top, aes(x = n_patient_treatments, y = feature, fill = feature_type))
    + geom_col()
    + labs(
        x = "# (patient, treatment) pairs where feature is treatment-only", y = NULL,
        fill = "feature type", title = paste("top", nrow(top), "recurrent treatment-only features")
    )
    + theme_classic(base_size = 10)
    + theme(axis.text.y = element_text(size = 7))
)
grid_df <- recurrence |> count(channel, feature_type, name = "n_features")
p_grid <- (
    ggplot(grid_df, aes(x = feature_type, y = channel, fill = n_features))
    + geom_tile()
    + geom_text(aes(label = n_features), colour = "white", size = 3)
    + scale_fill_viridis_c(name = "# distinct features")
    + labs(x = NULL, y = NULL, title = "distinct treatment-only features: channel x feature type")
    + theme_minimal(base_size = 10)
    + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
)
add_page((p_bar + p_grid) + plot_layout(widths = c(1.4, 1)), 18, max(6, 0.28 * nrow(top)))

per_treatment <- read_result("treatment_only_feature_by_treatment.parquet") |>
    pivot_longer(-feature, names_to = "treatment", values_to = "n_patients") |>
    mutate(feature = factor(feature, levels = rev(top$feature)))
p_by_trt <- (
    ggplot(per_treatment, aes(x = treatment, y = feature, fill = n_patients))
    + geom_tile()
    + scale_fill_viridis_c(name = "# patients (feature is treatment-only)")
    + labs(x = NULL, y = NULL, title = "top treatment-only features across treatments")
    + theme_minimal(base_size = 10)
    + theme(panel.grid = element_blank(), axis.text.y = element_text(size = 7), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
)
add_page(p_by_trt, max(8, 0.4 * length(unique(per_treatment$treatment)) + 5), max(6, 0.28 * nrow(top)))

# Task D: distribution of the top shared treatment-only features
shared <- read_result("top_shared_features_coefficients.parquet")
top_shared <- head(recurrence$feature, N_SHARED)
shared <- shared |> mutate(feature = factor(feature, levels = rev(top_shared)))
hit_pal <- c("TRUE" = "#d62728", "FALSE" = "#7f7f7f")
p_box <- (
    ggplot(shared, aes(x = coefficient, y = feature))
    + geom_boxplot(fill = "lightgrey", outlier.shape = NA)
    + geom_jitter(aes(colour = treatment_only), height = 0.2, size = 0.6, alpha = 0.5)
    + geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4)
    + scale_colour_manual(values = hit_pal, name = "treatment-only hit")
    + labs(x = "treatment coefficient (all patient x treatment models)", y = NULL)
    + theme_classic(base_size = 10)
)
frac <- shared |> group_by(feature) |> summarise(fraction = mean(treatment_only), .groups = "drop")
p_frac <- (
    ggplot(frac, aes(x = fraction, y = feature))
    + geom_col(fill = "#d62728")
    + labs(x = "fraction of models that are treatment-only hits", y = NULL)
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
    n_hits <- shared |> filter(treatment_only) |> count(feature, .data[[by]], name = "n_hits")
    med <- med |> left_join(n_hits, by = c("feature", by)) |> mutate(n_hits = ifelse(is.na(n_hits), 0L, n_hits))
    lim <- quantile(abs(med$median_coef), 0.98, na.rm = TRUE)
    p_med <- (
        ggplot(med, aes(x = .data[[by]], y = feature, fill = median_coef))
        + geom_tile()
        + scale_fill_gradient2(low = "#3b4cc0", mid = "white", high = "#b40426", midpoint = 0, limits = c(-lim, lim), oob = squish, name = "median treatment coefficient")
        + labs(
            x = NULL, y = NULL,
            title = paste0("top ", N_SHARED, " shared treatment-only features across ", by, "s", if (by == "patient") " (numbers = treatment-only hits)" else "")
        )
        + theme_minimal(base_size = 10)
        + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    )
    if (by == "patient") p_med <- p_med + geom_text(aes(label = n_hits), size = 2.5)
    add_page(p_med, if (by == "patient") 8 else max(9, 0.4 * length(treatments) + 5), 5)
}

save_plots_pdf(pages, tag_pdf, width = page_w, height = page_h)
cat(length(pages), "pages ->", tag_pdf, "\n")
