list_of_packages <- c("ggplot2", "dplyr", "tidyr", "patchwork", "arrow", "ggrastr", "RColorBrewer", "png", "grid")
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

module_dir <- file.path(root_dir, "4.linear_modeling")
results_path <- file.path(module_dir, "results")
figures_path <- file.path(module_dir, "figures", "headline_results_figure")
dir.create(figures_path, recursive = TRUE, showWarnings = FALSE)

profiles <- c("organoid", "sc", "organoid_agg", "sc_agg")
main_profiles <- c("organoid", "sc")
profile_palette <- setNames(c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), profiles)

# # Summary figure: the most load-bearing results across the linear-modeling pipeline
#
# One panel per question, all drawn from tables/figures already saved by the calculation
# notebooks (5.calculate_variate_importance) -- nothing is recomputed here.
#
# A. Where does the variance go (treatment vs technical covariates)?
# B. How many/how strong are the treatment effects (volcano)?
# C. Tumor-type UpSet (title-free panel from 6.plot_variate_importance.r)
# D. Treatment-only variate co-occurrence (title-free panel from 6.plot_variate_importance.r)
#
# Panels C and D embed pre-rendered, title-free PNGs saved by 6.plot_variate_importance.r's
# "Headline-figure panels" section (the standalone PDFs in figures/variate_importance/ keep
# their titles; only these dedicated panel PNGs omit them).

# ---- A. variance partition -------------------------------------------------
variate_levels <- c(
    "treatment", "cell_count", "organoid_count", "cell_per_organoid_count",
    "manhattan_distance_from_center", "cell_x_position", "cell_y_position",
    "cell_z_position", "cell_z_depth", "shared/unattributed", "residual"
)
paired <- brewer.pal(12, "Paired")
variate_palette <- setNames(
    c(paired[c(2, 4, 6, 8, 10, 12, 1, 3, 5)], "#7f7f7f", "#d9d9d9"),
    variate_levels
)
partition <- arrow::read_parquet(
    file.path(results_path, "explore_linear_models", "variance_partition_mean.parquet")
) |>
    pivot_longer(-profile, names_to = "term", values_to = "pct") |>
    mutate(profile = factor(profile, levels = rev(profiles)), term = factor(term, levels = variate_levels))
pA <- (
    ggplot(partition, aes(x = pct, y = profile, fill = term))
    + geom_col(position = position_stack(reverse = TRUE))
    + scale_fill_manual(values = variate_palette, drop = FALSE)
    + labs(x = "mean % of total variance", y = NULL, fill = NULL)
    + theme_manuscript(base_size = 14)
    + theme(legend.position = "bottom", legend.text = element_text(size = 13), legend.key.size = unit(6, "mm"))
)

# ---- B. volcano (treatment effect size vs significance) -------------------
volcano <- arrow::read_parquet(
    file.path(results_path, "explore_linear_models", "plot_data", "4_volcano_by_profile__volcano.parquet")
) |>
    filter(profile %in% main_profiles) |>
    mutate(profile = factor(profile, levels = main_profiles))
pB <- (
    ggplot(volcano, aes(x = coefficient, y = neglog_fdr))
    + rasterise(geom_point(aes(colour = rsquared), size = 1, alpha = 0.3), dpi = 600)
    + geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3)
    + scale_colour_viridis_c(option = "A", name = "R2", limits = c(0, 1))
    + facet_wrap(~profile, ncol = 2, scales = "free")
    + labs(x = "treatment coefficient", y = "-log10 FDR")
    + theme_manuscript(base_size = 14)
)

# ---- C. treatment-only variate co-occurrence (title-free panel from 6.plot_variate_importance.r) --
cooccurrence_png <- file.path(
    module_dir, "figures", "variate_importance", "headline_panels",
    "cooccurrence_organoid_original.png"
)
if (!file.exists(cooccurrence_png)) {
    stop(
        "missing ", cooccurrence_png, " -- run scripts/6.plot_variate_importance.r first ",
        "(it saves this PNG, which panel C renders)"
    )
}
pC <- wrap_elements(full = grid::rasterGrob(png::readPNG(cooccurrence_png), interpolate = TRUE))

# ---- D. tumor-type UpSet (title-free panel from 6.plot_variate_importance.r: original model, tumor_type=pNF) --
tumor_type_png <- file.path(
    module_dir, "figures", "variate_importance", "headline_panels",
    "tumor_type_pNF_upset.png"
)
if (!file.exists(tumor_type_png)) {
    stop(
        "missing ", tumor_type_png, " -- run scripts/6.plot_variate_importance.r first ",
        "(it saves this PNG, which panel D renders)"
    )
}
pD <- wrap_elements(full = grid::rasterGrob(png::readPNG(tumor_type_png), interpolate = TRUE))

# ---- assemble ---------------------------------------------------------------
# rows 2-3 give the embedded UpSet (C) and heatmap (D) twice the height of A and B
layout <- "
AABB
CCDD
CCDD
"
headline_results_figure <- (
    wrap_plots(A = pA, B = pB, C = pD, D = pC, design = layout)
    + plot_annotation(tag_levels = "A")
)
save_ggplot(headline_results_figure, file.path(figures_path, "headline_results_figure.png"), width = 16, height = 16, dpi = 600)
