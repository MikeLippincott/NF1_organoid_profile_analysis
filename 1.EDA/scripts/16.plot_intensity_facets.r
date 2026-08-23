list_of_packages <- c("ggplot2", "dplyr", "arrow", "RColorBrewer")
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

results_dir <- file.path(root_dir, "1.EDA", "results", "intensity")
figures_dir <- file.path(root_dir, "1.EDA", "figures", "intensity")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

plot_theme <- theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 8),
    legend.position = "none"
)

make_stat_figures <- function(df, modality_label, out_subdir, x_col_patient, x_col_treatment) {
    dir.create(file.path(figures_dir, out_subdir), recursive = TRUE, showWarnings = FALSE)
    stats <- unique(df$stat)
    for (s in stats) {
        d <- df %>% filter(stat == s)

        # per-patient: distribution of raw values per patient x channel x compartment
        p_patient <- (
            ggplot(d, aes(x = .data[[x_col_patient]], y = value, fill = .data[[x_col_patient]]))
            + geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.4)
            + facet_grid(channel ~ compartment, scales = "free_y")
            + labs(title = paste0(modality_label, ": ", s, " by patient, faceted by compartment x channel"),
                   x = "Patient", y = paste0(s, " (z-scored)"))
            + plot_theme
        )
        fname <- gsub("[^A-Za-z0-9]+", "_", s)
        ggsave(filename = file.path(figures_dir, out_subdir, paste0(fname, "_per_patient.png")),
               plot = p_patient, width = 11, height = 8, dpi = 600, units = "in")

        # pooled: distribution of raw values per treatment x channel x compartment (all patients combined)
        d_treatment <- d
        d_treatment[[x_col_treatment]] <- factor(d_treatment[[x_col_treatment]],
                                                    levels = intersect(custom_treatment_order, unique(d_treatment[[x_col_treatment]])))
        p_pooled <- (
            ggplot(d_treatment, aes(x = .data[[x_col_treatment]], y = value, fill = .data[[x_col_treatment]]))
            + geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.4)
            + scale_fill_manual(values = custom_treatment_palette, na.value = "grey70")
            + facet_grid(channel ~ compartment, scales = "free_y")
            + labs(title = paste0(modality_label, " pooled: ", s, " by treatment, faceted by compartment x channel"),
                   x = "Treatment", y = paste0(s, " (z-scored)"))
            + plot_theme
        )
        ggsave(filename = file.path(figures_dir, out_subdir, paste0(fname, "_pooled.png")),
               plot = p_pooled, width = 12, height = 8, dpi = 600, units = "in")
    }
    cat("Wrote", 2 * length(stats), "figures to", file.path(figures_dir, out_subdir), "\n")
}

# --- 3D ---
intensity_3d <- read_parquet(file.path(results_dir, "intensity_values_3D.parquet"))
make_stat_figures(intensity_3d, "3D", "3D", "Metadata_patient", "Metadata_treatment")

# --- 2D: one subdirectory per projection method ---
intensity_2d <- read_parquet(file.path(results_dir, "intensity_values_2D.parquet"))
for (proj in unique(intensity_2d$projection)) {
    d_proj <- intensity_2d %>% filter(projection == proj)
    make_stat_figures(d_proj, paste0("2D (", proj, ")"), file.path("2D", proj), "Metadata_patient", "Metadata_treatment")
}

cat("Done.\n")
