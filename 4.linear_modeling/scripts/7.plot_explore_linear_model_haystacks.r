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

plot_data_path <- file.path(root_dir, "4.linear_modeling/results/explore_linear_models/plot_data")
figures_path <- file.path(root_dir, "4.linear_modeling/figures/explore_linear_models")
dir.create(figures_path, recursive = TRUE, showWarnings = FALSE)
pdf_path <- file.path(figures_path, "explore_linear_models.pdf")

# tables saved by 4.explore_linear_model_haystacks for each figure: plot_data/<figure>__<key>.parquet
read_plot <- function(figure, key) {
    arrow::read_parquet(file.path(plot_data_path, paste0(figure, "__", key, ".parquet")))
}

# "hit" definition threshold, identical to the calculation notebook
FDR_MAX <- 0.05
profiles <- c("organoid", "sc", "organoid_agg", "sc_agg")
main_profiles <- c("organoid", "sc")
profile_palette <- setNames(c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), profiles)

# every figure is a page of one pdf, in the order it is added
pages <- list()
page_w <- c()
page_h <- c()
add_page <- function(plot, width, height) {
    pages[[length(pages) + 1]] <<- plot
    page_w <<- c(page_w, width)
    page_h <<- c(page_h, height)
}

heat_fill <- function(name, ...) scale_fill_viridis_c(option = "F", direction = -1, name = name, ...)
heat_theme <- theme_minimal(base_size = 11) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# 4. volcano plots: significance on the y-axis, effect size on the x-axis, hits coloured by R2
volcano <- read_plot("4_volcano_by_profile", "volcano") |> mutate(profile = factor(profile, levels = profiles))
p_volcano <- (
    ggplot(volcano, aes(x = coefficient, y = neglog_fdr))
    + rasterise(geom_point(data = filter(volcano, !hit), colour = "lightgrey", size = 0.3), dpi = 600)
    + rasterise(geom_point(data = filter(volcano, hit), aes(colour = rsquared), size = 0.5), dpi = 600)
    + geom_text_repel(data = filter(volcano, !is.na(label)), aes(label = label), size = 1.8, max.overlaps = Inf, segment.size = 0.2)
    + geom_hline(yintercept = -log10(FDR_MAX), linetype = "dashed", linewidth = 0.4)
    + scale_colour_viridis_c(option = "A", name = "R2")
    + coord_cartesian(xlim = c(-5, 5))
    + facet_wrap(~profile, ncol = 2)
    + labs(x = "treatment coefficient", y = "-log10 FDR", title = "Volcano plots (colour = R2 of hits)")
    + theme_bw(base_size = 12)
)
add_page(p_volcano, 18, 14)

# 5. variance partitioning: where does the variance go?
partition <- read_plot("5_variance_partition_by_profile", "partition") |>
    pivot_longer(-profile, names_to = "term", values_to = "pct") |>
    mutate(profile = factor(profile, levels = rev(unique(profile))))
p_partition <- (
    ggplot(partition, aes(x = pct, y = profile, fill = term))
    + geom_col()
    + scale_fill_manual(values = hue_pal()(length(unique(partition$term))))
    + labs(x = "mean % of total variance", y = NULL, fill = NULL, title = "Where does the variance go?")
    + theme_classic(base_size = 12)
)
add_page(p_partition, 12, 4)

share <- read_plot("5_treatment_variance_share_heatmap", "share") |> mutate(profile = factor(profile, levels = main_profiles))
p_share <- (
    ggplot(share, aes(x = treatment, y = patient, fill = median_pct_var))
    + geom_tile()
    + heat_fill("median treatment % var")
    + facet_wrap(~profile, ncol = 2)
    + labs(x = NULL, y = NULL, title = "How much does treatment explain?")
    + heat_theme
)
add_page(p_share, 22, 7)

# 6. hit landscape: which patients, drugs, MOAs, doses and tumor types produce the most hits?
landscape <- read_plot("6_hit_rate_by_group", "landscape")
levels_order <- c("patient", "treatment", "drug", "therapeutic_category", "tumor_type", "dose")
landscape <- landscape |> mutate(level = factor(level, levels = levels_order), profile = factor(profile, levels = profiles))
# groups ordered by the organoid hit rate within each level (groups the organoid profile lacks go last)
group_order <- landscape |>
    filter(profile == "organoid") |>
    arrange(level, desc(hit_rate)) |>
    mutate(key = paste(level, group, sep = "___"))
landscape <- landscape |> mutate(key = paste(level, group, sep = "___"))
key_levels <- c(rev(group_order$key), setdiff(unique(landscape$key), group_order$key))
landscape$key <- factor(landscape$key, levels = key_levels)
p_landscape <- (
    ggplot(landscape, aes(x = hit_rate, y = key, fill = profile))
    + geom_col(position = position_dodge(width = 0.8), width = 0.75)
    + scale_y_discrete(labels = function(x) sub("^.*___", "", x))
    + scale_fill_manual(values = profile_palette)
    + facet_wrap(~level, scales = "free_y", ncol = 3)
    + labs(x = "hit rate", y = NULL, fill = "profile", title = "Hit rate by group")
    + theme_classic(base_size = 11)
)
add_page(p_landscape, 22, 11)

hit_rate_pt <- read_plot("6_hit_rate_patient_by_treatment", "hit_rate") |> mutate(profile = factor(profile, levels = main_profiles))
p_hit_rate_pt <- (
    ggplot(hit_rate_pt, aes(x = treatment, y = patient, fill = hit_rate))
    + geom_tile()
    + heat_fill("hit rate")
    + facet_wrap(~profile, ncol = 2)
    + labs(x = NULL, y = NULL, title = "Hit rate, patient x treatment")
    + heat_theme
)
add_page(p_hit_rate_pt, 22, 7)

direction <- read_plot("6_hit_direction_by_treatment", "direction") |>
    pivot_longer(c(down, up), names_to = "direction", values_to = "n_hits") |>
    mutate(profile = factor(profile, levels = main_profiles), direction = factor(direction, levels = c("down", "up")))
p_direction <- (
    ggplot(direction, aes(x = n_hits, y = treatment, fill = direction))
    + geom_col()
    + geom_vline(xintercept = 0, linewidth = 0.3)
    + scale_fill_manual(values = c(down = "#4575b4", up = "#d73027"))
    + facet_wrap(~profile, ncol = 2)
    + labs(x = "number of hits (down < 0 < up)", y = NULL, title = "Hits up vs down vs DMSO")
    + theme_classic(base_size = 11)
)
add_page(p_direction, 20, 7)

# 7. feature space enrichment: Fisher exact test of hit vs non-hit for every compartment / channel / feature type
enrichment <- read_plot("7_feature_family_enrichment", "enrichment") |>
    filter(family %in% c("compartment_label", "channel_label", "Feature_type")) |>
    mutate(
        profile = factor(profile, levels = main_profiles),
        family = factor(family, levels = c("compartment_label", "channel_label", "Feature_type")),
        colour = ifelse(fdr < 0.05, ifelse(log2_or > 0, "enriched", "depleted"), "n.s."),
        key = paste(profile, family, level, sep = "___")
    ) |>
    group_by(profile, family) |>
    arrange(log2_or, .by_group = TRUE) |>
    ungroup()
enrichment$key <- factor(enrichment$key, levels = enrichment$key)
p_enrichment <- (
    ggplot(enrichment, aes(x = log2_or, y = key, fill = colour))
    + geom_col()
    + geom_vline(xintercept = 0, linewidth = 0.3)
    + scale_y_discrete(labels = function(x) sub("^.*___", "", x))
    + scale_fill_manual(values = c(enriched = "#d73027", depleted = "#4575b4", "n.s." = "grey"), name = "FDR < 0.05")
    + facet_grid(profile ~ family, scales = "free_y")
    + labs(x = "log2 odds ratio of being a hit", y = NULL, title = "Feature family enrichment")
    + theme_classic(base_size = 11)
)
add_page(p_enrichment, 22, 12)

channel_heat <- read_plot("7_channel_by_feature_type_heatmaps", "heatmap") |> mutate(profile = factor(profile, levels = main_profiles))
p_channel <- (
    ggplot(channel_heat, aes(x = column, y = channel, fill = hit_rate))
    + geom_tile()
    + geom_text(aes(label = sprintf("%.3f", hit_rate)), size = 2.5, colour = "grey40")
    + heat_fill("hit rate")
    + facet_grid(profile ~ kind, scales = "free_x")
    + labs(x = NULL, y = NULL, title = "Hit rate: channel x feature type / compartment")
    + heat_theme
)
add_page(p_channel, 20, 12)

drug_channel <- read_plot("7_drug_by_channel_hit_rate", "heatmap") |> mutate(profile = factor(profile, levels = main_profiles))
p_drug_channel <- (
    ggplot(drug_channel, aes(x = channel_label, y = drug, fill = hit_rate))
    + geom_tile()
    + heat_fill("hit rate")
    + facet_wrap(~profile, ncol = 2)
    + labs(x = NULL, y = NULL, title = "Drug x channel")
    + heat_theme
)
add_page(p_drug_channel, 22, 8)

# 8. dose response: do two doses of the same drug agree?
concordance <- read_plot("8_dose_concordance", "concordance") |> mutate(profile = factor(profile, levels = profiles))
p_concordance <- (
    ggplot(concordance, aes(x = drug, y = spearman, fill = profile))
    + geom_col(position = position_dodge(width = 0.8), width = 0.75)
    + scale_fill_manual(values = profile_palette)
    + labs(x = NULL, y = "Spearman", title = "Coefficient concordance between two doses of the same drug")
    + theme_classic(base_size = 11)
    + theme(axis.text.x = element_text(angle = 45, hjust = 1))
)
add_page(p_concordance, 10, 5)

scatter_path <- file.path(plot_data_path, "8_dose_scatter_organoid__scatter.parquet")
if (file.exists(scatter_path)) {
    dose_scatter <- arrow::read_parquet(scatter_path) |>
        mutate(
            coef_low = pmin(pmax(coef_low, -3), 3), coef_high = pmin(pmax(coef_high, -3), 3),
            panel = paste0(drug, "\n", dose_low, " vs ", dose_high)
        )
    p_scatter <- (
        ggplot(dose_scatter, aes(x = coef_low, y = coef_high))
        + geom_bin2d(bins = 50)
        + geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed", linewidth = 0.4)
        + scale_fill_viridis_c(trans = "log10", name = "count")
        + facet_wrap(~panel, nrow = 1)
        + labs(x = "coef @ low dose", y = "coef @ high dose", title = "Organoid: low vs high dose")
        + theme_bw(base_size = 11)
    )
    add_page(p_scatter, 5 * max(length(unique(dose_scatter$drug)), 1), 5)
}

save_plots_pdf(pages, pdf_path, width = page_w, height = page_h)
cat(length(pages), "pages ->", pdf_path, "\n")
