list_of_packages <- c("ggplot2", "dplyr", "arrow")

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
cat("Git root directory:", root_dir, "\n")

results_dir <- file.path(root_dir, "2.2d_vs_3d_analysis/results/sparse_cca")
figures_dir <- file.path(root_dir, "2.2d_vs_3d_analysis/figures/sparse_cca")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

canonical_scores <- arrow::read_parquet(file.path(results_dir, "canonical_scores.parquet"))
cor_summary <- arrow::read_parquet(file.path(results_dir, "canonical_correlations.parquet"))

cc_cols <- grep("^CC[0-9]+_2D$", colnames(canonical_scores), value = TRUE)
n_components <- length(cc_cols)

cor_bar_plot <- ggplot(cor_summary, aes(x = component, y = canonical_correlation)) +
    geom_col(fill = "#882E8B") +
    labs(
        title = "Sparse CCA: Canonical Correlations (2D vs 3D)",
        x = "Component",
        y = "Canonical correlation (r)"
    ) +
    theme_bw() +
    theme(text = element_text(size = 14)) +
    scale_y_continuous(limits = c(0, 1))

print(cor_bar_plot)
ggsave(
    file.path(figures_dir, "sparse_cca_canonical_correlations.png"),
    cor_bar_plot, width = 6, height = 6, dpi = 300, units = "in"
)

label_map <- setNames(
    paste0(cor_summary$component, "\n(r = ", round(cor_summary$canonical_correlation, 2), ")"),
    cor_summary$component
)

score_long <- do.call(rbind, lapply(seq_len(n_components), function(k) {
    comp <- paste0("CC", k)
    data.frame(
        Metadata_patient = canonical_scores$Metadata_patient,
        component = comp,
        component_label = label_map[[comp]],
        score_2d = canonical_scores[[paste0(comp, "_2D")]],
        score_3d = canonical_scores[[paste0(comp, "_3D")]]
    )
}))

score_scatter_plot <- ggplot(score_long, aes(x = score_2d, y = score_3d, color = Metadata_patient)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) +
    facet_wrap(~component_label, scales = "free") +
    labs(
        title = "Sparse CCA: 2D vs 3D Canonical Scores",
        x = "2D canonical score",
        y = "3D canonical score",
        color = "Patient"
    ) +
    theme_bw() +
    theme(text = element_text(size = 14))

print(score_scatter_plot)
ggsave(
    file.path(figures_dir, "sparse_cca_canonical_score_scatter.png"),
    score_scatter_plot, width = 12, height = 5, dpi = 300, units = "in"
)
