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
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(size = 8),
    legend.position = "none"
)

# Biological question: how does intensity, per channel and per patient,
# differ between MEK1/2 inhibitors and DMSO control? Restricted to 3D and
# this treatment class, and split one compartment per plot (rather than
# dodging all compartments into one crowded plot) so each panel stays
# readable -- the full treatment panel and the 2D projection-method
# comparisons didn't serve this question and were dropped.
#
# `value` is z-scored PER PATIENT upstream, so Metadata_patient is always a
# facet (never an x-axis category) with free y-scales -- every panel is
# valid on its own terms. Outlier points on the boxplots are suppressed (the
# violin already shows the tails) to keep the file size down.

mek_treatments <- names(treatment_moa_map)[treatment_moa_map == "MEK1/2 inhibitor"]
keep_treatments <- c("DMSO", mek_treatments)

intensity_3d <- read_parquet(file.path(results_dir, "intensity_values_3D.parquet")) %>%
    filter(Metadata_treatment %in% keep_treatments)
intensity_3d$Metadata_treatment <- factor(intensity_3d$Metadata_treatment,
                                            levels = intersect(custom_treatment_order, keep_treatments))

pdf(file.path(figures_dir, "3D_intensity_MEK_vs_DMSO.pdf"), width = 16, height = 8, onefile = TRUE)
for (s in unique(intensity_3d$stat)) {
    for (cmp in unique(intensity_3d$compartment)) {
        d <- intensity_3d %>% filter(stat == s, compartment == cmp)
        # A handful of rows carry extreme z-scores (near-zero within-patient
        # variance blows up the z-score for the rare non-zero raw value --
        # e.g. ER MedianIntensity hits 1e21). Left alone, one such point
        # forces free_y to stretch to cover it, flattening every other
        # channel/patient panel to a line at 0. Clip per patient x channel
        # to the 1st/99th percentile so panel scales reflect the bulk of
        # the distribution; this only affects the plot, not the data file.
        d <- d %>%
            group_by(Metadata_patient, channel) %>%
            mutate(
                value = pmin(pmax(value, quantile(value, 0.01, na.rm = TRUE)),
                             quantile(value, 0.99, na.rm = TRUE))
            ) %>%
            ungroup()
        p <- (
            ggplot(d, aes(x = Metadata_treatment, y = value, fill = Metadata_treatment))
            + geom_violin(alpha = 0.6, trim = TRUE, scale = "width")
            + geom_boxplot(width = 0.15, alpha = 0.85, outlier.shape = NA)
            + scale_fill_manual(values = custom_treatment_palette, na.value = "grey70")
            + facet_grid(channel ~ Metadata_patient, scales = "free_y")
            + labs(
                title = paste0("3D (", cmp, "): ", s, ", MEK inhibitors vs. DMSO, by patient x channel"),
                x = "Treatment", y = paste0(s, " (z-scored within patient)")
            )
            + plot_theme
        )
        print(p)
    }
}
dev.off()

cat("Wrote 1 PDF to", figures_dir, "\n")
