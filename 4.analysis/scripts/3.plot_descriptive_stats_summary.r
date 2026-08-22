list_of_packages <- c("ggplot2", "dplyr", "readr", "RColorBrewer")
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

results_dir <- file.path(root_dir, "4.analysis", "results", "descriptive_stats")
figures_dir <- file.path(root_dir, "4.analysis", "figures", "descriptive_stats")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

summary_df <- read_csv(
    file.path(results_dir, "feature_significance_summary.csv"),
    show_col_types = FALSE,
    col_types = cols(modality = col_character(), projection = col_character())
)
summary_df$pct_significant <- summary_df$n_significant / summary_df$n_tested * 100

local_feature_type_palette <- c(
    feature_type_palette,
    "AreaShape" = feature_type_palette[["AreaSizeShape"]],
    "RadialDistribution" = brewer.pal(8, "Paired")[4],
    "Neighbors" = brewer.pal(8, "Paired")[6],
    "Correlation" = brewer.pal(8, "Paired")[7],
    "Location" = "grey50"
)

plot_theme <- theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 11)
)

# --- 3D: per-patient % significant features by category ---
df_3d <- summary_df %>% filter(modality == "3D")
p_3d <- ggplot(df_3d, aes(x = Metadata_patient, y = pct_significant, fill = category)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = local_feature_type_palette, na.value = "grey70") +
    labs(
        title = "3D: % features significant across treatments (Kruskal-Wallis, FDR<0.05)",
        x = "Patient", y = "% significant features", fill = "Feature category"
    ) +
    plot_theme
ggsave(filename = file.path(figures_dir, "3D_significant_features_per_patient.png"),
       plot = p_3d, width = 10, height = 6, dpi = 600, units = "in")

# --- 3D pooled (all patients combined) ---
df_3d_pooled <- df_3d %>%
    group_by(category) %>%
    summarise(n_significant = sum(n_significant), n_tested = sum(n_tested), .groups = "drop") %>%
    mutate(pct_significant = n_significant / n_tested * 100)
p_3d_pooled <- ggplot(df_3d_pooled, aes(x = category, y = pct_significant, fill = category)) +
    geom_col() +
    scale_fill_manual(values = local_feature_type_palette, na.value = "grey70") +
    labs(title = "3D pooled (all patients): % features significant", x = "Feature category",
         y = "% significant features") +
    plot_theme + theme(legend.position = "none")
ggsave(filename = file.path(figures_dir, "3D_significant_features_pooled.png"),
       plot = p_3d_pooled, width = 8, height = 6, dpi = 600, units = "in")

# --- 2D: per-patient % significant features by category, faceted by projection ---
df_2d <- summary_df %>% filter(modality == "2D")
p_2d <- ggplot(df_2d, aes(x = Metadata_patient, y = pct_significant, fill = category)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = local_feature_type_palette, na.value = "grey70") +
    facet_wrap(~projection, ncol = 1) +
    labs(
        title = "2D: % features significant across treatments, by projection method",
        x = "Patient", y = "% significant features", fill = "Feature category"
    ) +
    plot_theme
ggsave(filename = file.path(figures_dir, "2D_significant_features_per_patient.png"),
       plot = p_2d, width = 10, height = 14, dpi = 600, units = "in")

# --- 2D pooled (all patients combined), faceted by projection ---
df_2d_pooled <- df_2d %>%
    group_by(projection, category) %>%
    summarise(n_significant = sum(n_significant), n_tested = sum(n_tested), .groups = "drop") %>%
    mutate(pct_significant = n_significant / n_tested * 100)
p_2d_pooled <- ggplot(df_2d_pooled, aes(x = category, y = pct_significant, fill = category)) +
    geom_col() +
    scale_fill_manual(values = local_feature_type_palette, na.value = "grey70") +
    facet_wrap(~projection) +
    labs(title = "2D pooled (all patients): % features significant, by projection method",
         x = "Feature category", y = "% significant features") +
    plot_theme + theme(legend.position = "none", axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
ggsave(filename = file.path(figures_dir, "2D_significant_features_pooled.png"),
       plot = p_2d_pooled, width = 12, height = 5, dpi = 600, units = "in")

cat("Wrote 4 figures to", figures_dir, "\n")
