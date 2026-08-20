# Shared plot builders so common chart shapes (density, boxplot, patterned
# bar/box, scatter comparison) aren't reimplemented per-notebook.
# All of them apply theme_manuscript() so styling stays consistent.

save_ggplot <- function(plot, path, width, height, dpi = 600) {
    #' Save a ggplot object and return it (invisible passthrough for chaining/printing).
    ggsave(filename = path, plot = plot, width = width, height = height, dpi = dpi)
    plot
}

plot_density_by_facet <- function(
    data, facet_col, facet_val, x_col, fill_col, y_max,
    x_max = NULL, palette = "RdBu", base_size = 18,
    x_lab = NULL, y_lab = "Density", fill_lab = NULL
) {
    if (is.null(x_lab)) x_lab <- "Normalized Object Counts per Well FOV"
    if (is.null(y_lab)) y_lab <- "Density"
    if (is.null(fill_lab)) fill_lab <- "Profile Type"
    #' Density plot for one facet level, sharing a common y-axis max across facets.
    plot_data <- data[data[[facet_col]] == facet_val, ]
    p <- (
        ggplot(plot_data, aes(x = .data[[x_col]], fill = .data[[fill_col]]))
        + geom_density(alpha = 0.3)
        + ylim(0, y_max)
        + labs(
            x = x_lab,
            y = y_lab,
            fill = paste0("Profile Type: ", facet_val)
        )
        + guides(fill = guide_legend(override.aes = list(alpha = 0.5)))
        + scale_fill_brewer(palette = palette)
        + theme_manuscript(base_size = base_size)
    )
    if (!is.null(x_max)) p <- p + xlim(0, x_max)
    p
}

plot_boxplot <- function(
    data, x_col, y_col, fill_col,
    x_lab, y_lab, fill_lab,
    palette = "RdBu", custom_palette = NULL,
    ylim_max = NULL, x_text = "blank", facet_formula = NULL, base_size = 18
) {
    #' Generic boxplot builder shared across the cell/organoid count summary plots.
    p <- (
        ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]], fill = .data[[fill_col]]))
        + geom_boxplot(outlier.shape = NA, alpha = 0.5)
        + labs(x = x_lab, y = y_lab, fill = fill_lab)
        + theme_manuscript(base_size = base_size, x_text = x_text)
    )

    p <- if (!is.null(custom_palette)) {
        p + scale_fill_manual(values = custom_palette)
    } else {
        p + scale_fill_brewer(palette = palette)
    }

    if (!is.null(ylim_max)) p <- p + ylim(0, ylim_max)
    if (!is.null(facet_formula)) p <- p + facet_wrap(facet_formula, scales = "free_x")
    p
}

plot_bar_pattern <- function(
    data, x_col, y_col, fill_col, pattern_col, fill_palette,
    x_lab, y_lab, fill_lab = "Profile Type", pattern_lab = "Dose",
    pattern_values = c("none", "stripe", "crosshatch", "circle"),
    ylim_max = NULL, x_text = "blank", facet_formula = NULL, base_size = 18
) {
    #' Mean bar plot with a dose-encoded pattern fill, shared across treatment bar plots.
    p <- (
        ggplot(
            data,
            aes(
                x = .data[[x_col]],
                y = .data[[y_col]],
                fill = .data[[fill_col]],
                pattern = factor(.data[[pattern_col]])
            )
        )
        + geom_bar_pattern(
            stat = "summary",
            fun = "mean",
            position = position_dodge(width = 0.9),
            color = "black",
            pattern_fill = "black",
            pattern_color = "black",
            pattern_density = 0.1,
            pattern_spacing = 0.03,
            pattern_alpha = 0.6,
            pattern_key_scale_factor = 2.5
        )
        + scale_fill_manual(values = fill_palette)
        + scale_pattern_manual(values = pattern_values)
        + labs(x = x_lab, y = y_lab, fill = fill_lab, pattern = pattern_lab)
        + guides(fill = guide_legend(override.aes = list(pattern = "none")))
        + theme_manuscript(base_size = base_size, x_text = x_text)
    )

    if (!is.null(ylim_max)) p <- p + ylim(0, ylim_max)
    if (!is.null(facet_formula)) p <- p + facet_wrap(facet_formula, scales = "free_y")
    p
}

plot_boxplot_pattern <- function(
    data, x_col, y_col, fill_col, pattern_col, fill_palette,
    x_lab, y_lab, fill_lab = "Profile Type", pattern_lab = "Dose",
    pattern_values = c("none", "stripe", "crosshatch", "circle"),
    ylim_max = NULL, x_text = "angled", facet_formula = NULL, base_size = 18
) {
    #' Boxplot with a dose-encoded pattern fill, shared across treatment box plots.
    p <- (
        ggplot(
            data,
            aes(
                x = .data[[x_col]],
                y = .data[[y_col]],
                fill = .data[[fill_col]],
                pattern = factor(.data[[pattern_col]])
            )
        )
        + geom_boxplot_pattern(
            outlier.shape = NA,
            alpha = 0.5,
            pattern_fill = "black",
            pattern_color = "black",
            pattern_density = 0.1,
            pattern_spacing = 0.03,
            pattern_alpha = 0.6
        )
        + scale_fill_manual(values = fill_palette)
        + scale_pattern_manual(values = pattern_values)
        + labs(x = x_lab, y = y_lab, fill = fill_lab, pattern = pattern_lab)
        + guides(fill = guide_legend(override.aes = list(pattern = "none")))
        + theme_manuscript(base_size = base_size, x_text = x_text)
    )

    if (!is.null(ylim_max)) p <- p + ylim(0, ylim_max)
    if (!is.null(facet_formula)) p <- p + facet_wrap(facet_formula, scales = "free_y")
    p
}

plot_2d_vs_3d_scatter <- function(
    data, x_col, y_col, color_col, shape_col, color_palette,
    x_lab, y_lab, color_lab = "Treatment", shape_lab = "Dose",
    facet_formula = NULL, base_size = 18
) {
    #' Scatter comparison of a 2D vs 3D count metric with a reference y=x line.
    p <- (
        ggplot(
            data,
            aes(
                x = .data[[x_col]],
                y = .data[[y_col]],
                color = .data[[color_col]],
                shape = factor(.data[[shape_col]])
            )
        )
        + geom_point(size = 4, alpha = 0.7)
        + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")
        + scale_color_manual(values = color_palette)
        + labs(x = x_lab, y = y_lab, color = color_lab, shape = shape_lab)
        + theme_manuscript(base_size = base_size)
    )

    if (!is.null(facet_formula)) p <- p + facet_wrap(facet_formula, scales = "free", ncol = 3)
    p
}
