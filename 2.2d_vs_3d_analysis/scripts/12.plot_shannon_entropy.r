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

entropy_dir <- file.path(root_dir, "2.2d_vs_3d_analysis", "results", "entropy")
figures_dir <- file.path(root_dir, "2.2d_vs_3d_analysis", "figures", "entropy")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Load data

all_files <- list.files(entropy_dir, pattern = "_entropy\\.parquet$", full.names = TRUE)
pooled_files <- all_files[!grepl("_per_patient_entropy\\.parquet$", all_files)]
per_patient_files <- all_files[grepl("_per_patient_entropy\\.parquet$", all_files)]

entropy_pooled <- bind_rows(lapply(pooled_files, arrow::read_parquet))
entropy_per_patient <- bind_rows(lapply(per_patient_files, arrow::read_parquet))

entropy_pooled$group_label <- factor(entropy_pooled$group_label)
entropy_per_patient$group_label <- factor(entropy_per_patient$group_label)

# Plot styling
entropy_theme <- theme(
    plot.title   = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 12, angle = 20, hjust = 1),
    axis.text.y  = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12)
)

stage_titles <- c(normal = "Normal", fs = "Feature-Selected", agg = "Aggregated")
resolution_titles <- c(organoid = "Organoid", sc = "Single-Cell")

make_entropy_plot <- function(df, title) {
    ggplot(df, aes(x = group_label, y = entropy, fill = group_label)) +
        geom_violin(alpha = 0.6, trim = FALSE) +
        geom_boxplot(width = 0.15, alpha = 0.85, outlier.size = 0.5, outlier.alpha = 0.4) +
        labs(
            title = title,
            x     = "Profile type",
            y     = "Shannon Entropy (bits)",
            fill  = "Profile type"
        ) +
        theme_bw() +
        entropy_theme +
        theme(legend.position = "none")
}

width <- 7
height <- 6
dpi <- 300

for (res in names(resolution_titles)) {
    for (stage in names(stage_titles)) {
        df_subset <- entropy_pooled %>% filter(resolution == res, stage == !!stage)
        if (nrow(df_subset) == 0) {
            cat("SKIPPED (no data):", res, stage, "\n")
            next
        }

        title <- paste0(
            "Shannon Entropy: ", resolution_titles[[res]],
            " Profiles (", stage_titles[[stage]], ")"
        )
        p <- make_entropy_plot(df_subset, title)
        print(p)
        ggsave(
            filename = file.path(figures_dir, paste0(res, "_", stage, "_entropy.png")),
            plot     = p,
            width    = width,
            height   = height,
            dpi      = dpi,
            units    = "in"
        )
    }
}

# Per patient plots

make_per_patient_entropy_plot <- function(df, title) {
    ggplot(df, aes(x = patient, y = entropy, fill = patient)) +
        geom_boxplot(alpha = 0.7, outlier.size = 0.4, outlier.alpha = 0.3) +
        facet_wrap(~group_label, scales = "free_x") +
        labs(
            title = title,
            x     = "Patient",
            y     = "Shannon Entropy (bits)",
            fill  = "Patient"
        ) +
        theme_bw() +
        entropy_theme +
        theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1)
        )
}

pp_width  <- 12
pp_height <- 8

for (res in names(resolution_titles)) {
    for (stage in names(stage_titles)) {
        df_subset <- entropy_per_patient %>% filter(resolution == res, stage == !!stage)
        if (nrow(df_subset) == 0) {
            cat("SKIPPED (no data):", res, stage, "\n")
            next
        }

        title <- paste0(
            "Per-Patient Shannon Entropy: ", resolution_titles[[res]],
            " Profiles (", stage_titles[[stage]], ")"
        )
        p <- make_per_patient_entropy_plot(df_subset, title)
        print(p)
        ggsave(
            filename = file.path(figures_dir, paste0(res, "_", stage, "_per_patient_entropy.png")),
            plot     = p,
            width    = pp_width,
            height   = pp_height,
            dpi      = dpi,
            units    = "in"
        )
    }
}
