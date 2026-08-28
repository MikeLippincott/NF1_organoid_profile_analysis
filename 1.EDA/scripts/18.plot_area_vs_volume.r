list_of_packages <- c("ggplot2", "dplyr", "arrow", "RColorBrewer", "ggrastr", "scales")
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

results_dir <- file.path(root_dir, "1.EDA", "results", "area_vs_volume")
figures_dir <- file.path(root_dir, "1.EDA", "figures", "area_vs_volume")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

plot_theme <- theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, size = 13),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    legend.position = "none"
)

# Biological question: does 2D area track 3D volume the way simple organoid
# geometry would predict? Restricted to max_projection -- the standard 2D
# representation -- rather than comparing all 3 projection methods
# side-by-side, and to one plot per kind (organoid/cell) in a single PDF.
#
# Plotted in RAW (unnormalized) units, read directly from each modality's
# pre-normalization stage (script 17) -- not the z-scored values used
# elsewhere in this repo. Area came out of its upstream pipeline already
# z-scored per patient while volume didn't (previously mislabeled
# "z-scored" in this script when it wasn't); rather than patch that
# mismatch, we sidestep it entirely by plotting both in their native scale.
# Both area (px^2) and volume (voxels) span ~3-4 orders of magnitude, so
# both axes are log10.
#
# There's no shared organoid ID between the 2D and 3D pipelines (separate
# segmentations per modality), so individual area/volume records can't be
# paired per-object. We instead randomly pair a capped sample of each
# patient's area values with a capped sample of that patient's volume
# values. This shows the joint occupied range of the two distributions
# under an independence assumption -- NOT a per-organoid correlation -- and
# is capped per patient to keep the 2D density estimate tractable.
set.seed(0)
pair_sample_n <- 500

load_raw <- function(kind) {
    area <- read_parquet(file.path(results_dir, paste0("area_2D_", kind, "_raw.parquet")))
    volume <- read_parquet(file.path(results_dir, paste0("volume_3D_", kind, "_raw.parquet")))
    list(area = area, volume = volume)
}

pair_by_patient <- function(area_df, volume_df) {
    shared_patients <- intersect(unique(area_df$Metadata_patient_tumor),
                                  unique(volume_df$Metadata_patient_tumor))
    bind_rows(lapply(shared_patients, function(p) {
        a <- area_df$area[area_df$Metadata_patient_tumor == p]
        v <- volume_df$volume[volume_df$Metadata_patient_tumor == p]
        a_sample <- sample(a, min(length(a), pair_sample_n))
        v_sample <- sample(v, min(length(v), pair_sample_n))
        expand.grid(area = a_sample, volume = v_sample)
    }))
}

pdf(file.path(figures_dir, "area_vs_volume_density.pdf"), width = 8, height = 7, onefile = TRUE)
for (kind in c("organoid", "cell")) {
    label <- if (kind == "organoid") "Organoid" else "Single-cell"
    d <- load_raw(kind)
    paired <- pair_by_patient(d$area, d$volume)

    p <- (
        ggplot(paired, aes(x = area, y = volume))
        + rasterise(geom_point(alpha = 0.02, size = 0.3, color = "grey40"), dpi = 300)
        + geom_density_2d(color = "#3D7DCC", linewidth = 0.4)
        + scale_x_log10(labels = label_number())
        + scale_y_log10(labels = label_number())
        + labs(
            title = paste0(label, ": area (2D, max projection) vs. volume (3D), raw units"),
            x = "Area (px^2, log10 scale)", y = "Volume (voxels, log10 scale)"
        )
        + plot_theme
    )
    print(p)
}
dev.off()

cat("Wrote 1 PDF to", figures_dir, "\n")
