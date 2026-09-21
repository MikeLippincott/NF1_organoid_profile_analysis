list_of_packages <- c("ggplot2", "dplyr", "tidyr", "patchwork", "arrow")
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

results_path <- file.path(root_dir, "4.linear_modeling/results/variate_importance")
figures_path <- file.path(root_dir, "4.linear_modeling/figures/variate_importance")
dir.create(figures_path, recursive = TRUE, showWarnings = FALSE)

venn_df <- arrow::read_parquet(file.path(results_path, "variate_hit_venn_regions_all_scopes.parquet"))
upset_df <- arrow::read_parquet(file.path(results_path, "variate_hit_upset_counts_all_scopes.parquet"))
sizes_df <- arrow::read_parquet(file.path(results_path, "variate_hit_set_sizes_all_scopes.parquet"))

# every scope: the columns that define one Venn + UpSet pair
scopes <- list(
    all_models = character(0),
    per_patient_treatment = c("patient", "treatment"),
    per_patient = c("patient"),
    per_treatment = c("treatment"),
    per_treatment_tumor_type = c("treatment", "tumor_type")
)
model_sets <- c("original", "technical")
UPSET_TOP_N <- 25

term_colors <- c(
    treatment = "#d95f02", cell_count = "#1b9e77", organoid_count = "#7570b3",
    cell_per_organoid_count = "#e7298a", manhattan_distance_from_center = "#66a61e",
    cell_x_position = "#e6ab02", cell_y_position = "#a6761d",
    cell_z_position = "#666666", cell_z_depth = "#1f78b4"
)

# four-set Venn layout (plot coordinates): ellipse x, y, width, height, angle
ellipses <- data.frame(
    set = 1:4,
    x = c(0.350, 0.450, 0.544, 0.644),
    y = c(0.400, 0.500, 0.500, 0.400),
    w = 0.72, h = 0.45,
    angle = c(140, 140, 40, 40)
)
# where the count of each region (bits in term order) is written
region_xy <- data.frame(
    region_key = c(
        "1000", "0100", "0010", "0001", "1100", "1010", "1001", "0110",
        "0101", "0011", "1110", "1101", "1011", "0111", "1111"
    ),
    x = c(0.14, 0.32, 0.68, 0.85, 0.23, 0.29, 0.50, 0.50, 0.71, 0.77, 0.35, 0.61, 0.39, 0.65, 0.50),
    y = c(0.42, 0.72, 0.72, 0.42, 0.59, 0.30, 0.17, 0.66, 0.30, 0.59, 0.50, 0.24, 0.24, 0.50, 0.38)
)
set_label_xy <- data.frame(x = c(0.02, 0.22, 0.78, 0.98), y = c(0.86, 0.97, 0.97, 0.86), hjust = c(0, 0.5, 0.5, 1))

ellipse_polygon <- function(x, y, w, h, angle, set, n = 200) {
    t <- seq(0, 2 * pi, length.out = n)
    theta <- angle * pi / 180
    data.frame(
        set = set,
        px = x + (w / 2) * cos(t) * cos(theta) - (h / 2) * sin(t) * sin(theta),
        py = y + (w / 2) * cos(t) * sin(theta) + (h / 2) * sin(t) * cos(theta)
    )
}

plot_venn <- function(regions, sizes, title) {
    # four-set Venn of feature membership
    terms <- sizes$term
    stopifnot(length(terms) == 4)
    polys <- do.call(rbind, lapply(seq_len(4), function(i) do.call(ellipse_polygon, as.list(ellipses[i, ]))))
    polys$term <- factor(terms[polys$set], levels = terms)
    counts <- region_xy |>
        left_join(regions[, c("region_key", "n_features")], by = "region_key") |>
        mutate(n_features = ifelse(is.na(n_features), 0, n_features))
    labels <- cbind(set_label_xy, term = terms, n = sizes$set_size)
    labels$label <- paste0(labels$term, "\n(", labels$n, " features)")
    (
        ggplot()
        + geom_polygon(data = polys, aes(px, py, group = term, fill = term), alpha = 0.3)
        + geom_path(data = polys, aes(px, py, group = term, colour = term), linewidth = 0.6)
        + geom_text(data = counts, aes(x, y, label = n_features), size = 5)
        + geom_text(
            data = labels, aes(x, y, label = label, colour = term, hjust = hjust),
            fontface = "bold", size = 4, lineheight = 0.9
        )
        + scale_fill_manual(values = term_colors[terms])
        + scale_colour_manual(values = term_colors[terms])
        + coord_fixed(xlim = c(0, 1), ylim = c(0, 1))
        + labs(title = paste0(
            "Features that are a hit for each variate (", title, ")\n",
            sizes$n_none[1], " of ", sizes$n_total[1],
            " features are a hit for none of these ", length(terms), " variates"
        ))
        + theme_void()
        + theme(legend.position = "none", plot.title = element_text(size = 11))
    )
}

plot_upset <- function(combos, sizes, title, top_n = UPSET_TOP_N) {
    # UpSet plot: bars = features in exactly that term combination, dots = which terms
    terms <- sizes$term
    n_set <- length(terms)
    shown <- head(combos, top_n)
    n_combo <- nrow(shown)
    shown$x <- seq_len(n_combo)
    note <- if (n_combo < nrow(combos)) paste0("top ", n_combo, " of ", nrow(combos), " combinations; ") else ""

    p_bar <- (
        ggplot(shown, aes(x = x))
        + geom_rect(
            aes(xmin = x - 0.4, xmax = x + 0.4, ymin = 0.8, ymax = n_features, fill = treatment_specific)
        )
        + geom_text(aes(y = n_features * 1.08, label = n_features), vjust = 0, size = 3.5)
        + scale_fill_manual(values = c("TRUE" = term_colors[["treatment"]], "FALSE" = "#555555"))
        + scale_y_log10(expand = c(0, 0))
        + coord_cartesian(xlim = c(0.4, n_combo + 0.6), ylim = c(0.8, max(shown$n_features) * 3))
        + labs(
            x = NULL, y = "features in exactly\nthis combination (log)",
            title = paste0(
                "Feature membership across variate groups, ", title, "\n(", note, "orange = treatment only)"
            )
        )
        + theme_classic(base_size = 11)
        + theme(
            legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.line.x = element_blank(), plot.title = element_text(size = 11)
        )
    )

    membership <- shown |>
        select(x, all_of(paste0("in_", terms))) |>
        pivot_longer(-x, names_to = "term", values_to = "on") |>
        mutate(term = sub("^in_", "", term), row = match(term, terms))
    stripes <- data.frame(row = seq_len(n_set), shade = seq_len(n_set) %% 2 == 1)
    lines_df <- membership |> filter(on) |> group_by(x) |> filter(n() > 1) |> ungroup()
    p_mat <- (
        ggplot()
        + geom_rect(
            data = filter(stripes, shade),
            aes(xmin = 0.4, xmax = n_combo + 0.6, ymin = row - 0.5, ymax = row + 0.5), fill = "#f2f2f2"
        )
        + geom_point(data = membership, aes(x, row), colour = "#dddddd", size = 4)
        + geom_line(data = lines_df, aes(x, row, group = x), linewidth = 0.9)
        + geom_point(data = filter(membership, on), aes(x, row, fill = term), colour = "black", shape = 21, size = 4)
        + scale_fill_manual(values = term_colors[terms])
        + scale_y_reverse(limits = c(n_set + 0.5, 0.5), expand = c(0, 0))
        + scale_x_continuous(limits = c(0.4, n_combo + 0.6), expand = c(0, 0))
        + theme_void()
        + theme(legend.position = "none")
    )

    sizes$row <- seq_len(n_set)
    p_size <- (
        ggplot(sizes, aes(y = row))
        + geom_rect(aes(xmin = 0, xmax = set_size, ymin = row - 0.4, ymax = row + 0.4), fill = "#555555")
        + geom_text(aes(x = set_size, label = paste0(set_size, " ")), hjust = 1, size = 3)
        + scale_x_reverse(limits = c(max(sizes$set_size) * 1.35, 0))
        + scale_y_reverse(limits = c(n_set + 0.5, 0.5), expand = c(0, 0))
        + labs(x = "features in group", y = NULL)
        + theme_classic(base_size = 10)
        + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
    )
    p_lab <- (
        ggplot(sizes, aes(y = row))
        + geom_text(aes(x = 0.02, label = term, colour = term), hjust = 0, fontface = "bold", size = 3.2)
        + scale_colour_manual(values = term_colors[terms])
        + scale_x_continuous(limits = c(0, 1), expand = c(0, 0))
        + scale_y_reverse(limits = c(n_set + 0.5, 0.5), expand = c(0, 0))
        + theme_void()
        + theme(legend.position = "none")
    )

    (
        (plot_spacer() + plot_spacer() + p_bar + p_size + p_lab + p_mat)
        + plot_layout(ncol = 3, widths = c(1.3, 3.0, max(n_combo, 4)), heights = c(2.2, n_set * 0.55))
    )
}

group_title <- function(group_cols, row) {
    if (length(group_cols) == 0) return("all models")
    paste(sprintf("%s=%s", group_cols, unlist(row[group_cols])), collapse = ", ")
}

# one multi-page pdf per scope (Venn then UpSet page for every group and model set); a pdf is only rebuilt when missing
for (scope in names(scopes)) {
    pdf_path <- if (scope == "all_models") {
        file.path(figures_path, "variate_importance.pdf")
    } else {
        file.path(figures_path, scope, paste0(scope, ".pdf"))
    }
    if (file.exists(pdf_path)) {
        cat("skipping", scope, "(pdf exists)\n")
        next
    }
    dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
    group_cols <- scopes[[scope]]
    plots <- list()
    widths <- c()
    heights <- c()
    for (model in model_sets) {
        sizes_scope <- sizes_df |> filter(model_set == model, scope == !!scope)
        groups <- sizes_scope |> select(all_of(group_cols)) |> distinct() |> arrange(across(all_of(group_cols)))
        if (length(group_cols) == 0) groups <- data.frame(dummy = 1)
        for (i in seq_len(nrow(groups))) {
            in_group <- function(df) {
                df <- df |> filter(model_set == model, scope == !!scope)
                for (col in group_cols) df <- df[df[[col]] == groups[[col]][i], ]
                df
            }
            title <- group_title(group_cols, groups[i, , drop = FALSE])
            sizes_venn <- in_group(sizes_df) |> filter(plot == "venn")
            sizes_upset <- in_group(sizes_df) |> filter(plot == "upset")
            combos <- in_group(upset_df)
            plots[[length(plots) + 1]] <- plot_venn(in_group(venn_df), sizes_venn, title)
            widths <- c(widths, 10)
            heights <- c(heights, 8)
            if (nrow(combos) > 0) {
                plots[[length(plots) + 1]] <- plot_upset(combos, sizes_upset, title)
                widths <- c(widths, max(10, 0.9 * min(nrow(combos), UPSET_TOP_N) + 6))
                heights <- c(heights, 3 + 0.6 * nrow(sizes_upset))
            }
        }
    }
    save_plots_pdf(plots, pdf_path, width = widths, height = heights)
    cat(scope, ":", length(plots), "pages ->", pdf_path, "\n")
}
