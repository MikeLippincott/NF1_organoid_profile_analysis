list_of_packages <- c("ggplot2", "dplyr", "tidyr", "circlize","arrow", "ggrepel")
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

# Get the current working directory and find Git root
find_git_root <- function() {
    # Get current working directory
    cwd <- getwd()

    # Check if current directory has .git
    if (dir.exists(file.path(cwd, ".git"))) {
        return(cwd)
    }

    # If not, search parent directories
    current_path <- cwd
    while (dirname(current_path) != current_path) {  # While not at root
        parent_path <- dirname(current_path)
        if (dir.exists(file.path(parent_path, ".git"))) {
            return(parent_path)
        }
        current_path <- parent_path
    }

    # If no Git root found, stop with error
    stop("No Git root directory found.")
}

# Find the Git root directory
root_dir <- find_git_root()
cat("Git root directory:", root_dir, "\n")



# Figure path
figures_path <- file.path("../figures/mAP")
if (!dir.exists(figures_path)) {
  dir.create(figures_path, recursive = TRUE)
}

#Load Data
mAP_2d_organoid_intra <- arrow::read_parquet(file.path(
    root_dir, "1.EDA/results/mAP/2d_organoid_intra_patient_mAP_by_dose.parquet"
))
mAP_3d_organoid_intra <- arrow::read_parquet(file.path(
    root_dir, "1.EDA/results/mAP/3d_organoid_intra_patient_mAP_by_dose.parquet"
))
mAP_2d_sc_intra <- arrow::read_parquet(file.path(
    root_dir, "1.EDA/results/mAP/2d_sc_intra_patient_mAP_by_dose.parquet"
))
mAP_3d_sc_intra <- arrow::read_parquet(file.path(
    root_dir, "1.EDA/results/mAP/3d_sc_intra_patient_mAP_by_dose.parquet"
))
mAP_2d_organoid_inter <- arrow::read_parquet(file.path(
    root_dir, "1.EDA/results/mAP/2d_organoid_inter_patient_mAP_by_dose.parquet"
))
mAP_3d_organoid_inter <- arrow::read_parquet(file.path(
    root_dir, "1.EDA/results/mAP/3d_organoid_inter_patient_mAP_by_dose.parquet"
))
mAP_2d_sc_inter <- arrow::read_parquet(file.path(
    root_dir, "1.EDA/results/mAP/2d_sc_inter_patient_mAP_by_dose.parquet"
))
mAP_3d_sc_inter <- arrow::read_parquet(file.path(
    root_dir, "1.EDA/results/mAP/3d_sc_inter_patient_mAP_by_dose.parquet"
))

width <- 12
height <- 6
options(repr.plot.width = width, repr.plot.height = height)
organoid_intra_plot <- (
    ggplot(data = organoid_fs_intra_patient_mAP_df,
    aes(
        x = mean_average_precision,
        y = `-log10(p-value)`,
        color = treatment
    ))
    + geom_point(size = 2, alpha = 0.7)
    + labs(
        x = "Mean Average Precision",
        y = "-log10(p-value)",
        title = "Intra-patient mAP for Organoid FS"
    )
    + geom_hline(yintercept = 1.3, linetype = "dashed", color = "red")  # Adjust threshold as needed
    + theme_bw()
    + theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 16)
    )
    + xlim(0, 1)
    + facet_wrap(~patient, ncol = 4,
    # scales = "free"
    )
        + geom_text_repel(
        aes(
            label = treatment
        ),
        size = 4,
        nudge_y = 0.1,
        show.legend = FALSE,
        segment.color = "grey50",
        segment.size = 0.2,
        box.padding = 0.5,
        point.padding = 0.5,
        max.overlaps = 10
    )
)
ggsave(
    filename = organoid_fs_intra_patient_mAP_plot_path,
    plot = organoid_intra_plot,
    width = width,
    height = height,
    dpi = 600
)
organoid_intra_plot

width <- 8
height <- 6
options(repr.plot.width = width, repr.plot.height = height)
organoid_inter_plot <- (
    ggplot(data = organoid_fs_inter_patient_mAP_df,
    aes(
        x = mean_average_precision,
        y = `-log10(p-value)`,
        color = treatment
    ))
    + geom_point(size = 4, alpha = 0.7)
    + labs(
        x = "Mean Average Precision",
        y = "-log10(p-value)",
        title = "Inter-patient mAP for Organoid FS"
    )
    + geom_hline(yintercept = 1.3, linetype = "dashed", color = "red")  # Adjust threshold as needed
    + theme_bw()
    + theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 16)
    )
        + geom_text_repel(
        aes(
            label = treatment
        ),
        size = 4,
        nudge_y = 0.1,
        show.legend = FALSE,
        segment.color = "grey50",
        segment.size = 0.2,
        box.padding = 0.5,
        point.padding = 0.5,
        max.overlaps = 10
    )
    + xlim(0, 1)
)
ggsave(
    filename = organoid_fs_inter_patient_mAP_plot_path,
    plot = organoid_inter_plot,
    width = width,
    height = height,
    dpi = 600
)
organoid_inter_plot

width <- 10
height <- 6
options(repr.plot.width = width, repr.plot.height = height)
sc_intra_plot <- (
    ggplot(data = sc_fs_intra_patient_mAP_df,
    aes(
        x = mean_average_precision,
        y = `-log10(p-value)`,
        color = treatment
    ))
    + geom_point(size = 2, alpha = 0.7)
    + labs(
        x = "Mean Average Precision",
        y = "-log10(p-value)",
        title = "Intra-patient mAP for single cell FS"
    )
    + geom_hline(yintercept = 1.3, linetype = "dashed", color = "red")  # Adjust threshold as needed
    + theme_bw()
    + theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 16)
    )
    + xlim(0, 1)
    + facet_wrap(~patient, ncol = 4,
    # scales = "free"
    )
        + geom_text_repel(
        aes(
            label = treatment
        ),
        size = 4,
        nudge_y = 0.1,
        show.legend = FALSE,
        segment.color = "grey50",
        segment.size = 0.2,
        box.padding = 0.5,
        point.padding = 0.5,
        max.overlaps = 8
    )
)
ggsave(
    filename = sc_fs_intra_patient_mAP_plot_path,
    plot = sc_intra_plot,
    width = width,
    height = height,
    dpi = 600
)
sc_intra_plot

width <- 12
height <- 8
options(repr.plot.width = width, repr.plot.height = height)
sc_inter_plot <- (
    ggplot(data = sc_fs_inter_patient_mAP_df,
    aes(
        x = mean_average_precision,
        y = `-log10(p-value)`,
        color = treatment
    ))
    + geom_point(size = 4, alpha = 0.7)
    + labs(
        x = "Mean Average Precision",
        y = "-log10(p-value)",
        title = "Inter-patient mAP for single cell FS"
    )
    + geom_hline(yintercept = 1.3, linetype = "dashed", color = "red")  # Adjust threshold as needed
    + theme_bw()
    + theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 16)
    )
    + xlim(0, 1)
    + geom_text_repel(
        aes(
            label = treatment
        ),
        size = 4,
        nudge_y = 0.1,
        show.legend = FALSE,
        segment.color = "grey50",
        segment.size = 0.2,
        box.padding = 0.5,
        point.padding = 0.5,
        max.overlaps = 10
    )

)
ggsave(
    filename = sc_fs_inter_patient_mAP_plot_path,
    plot = sc_inter_plot,
    width = width,
    height = height,
    dpi = 600
)
sc_inter_plot

organoid_fs_intra_patient_mAP_df <- select(organoid_fs_intra_patient_mAP_df,
    c(
        treatment,mean_average_precision, patient
    )
)
colnames(organoid_fs_intra_patient_mAP_df)[2] <- "intra_patient_mAP"


organoid_fs_inter_patient_mAP_df <- select(organoid_fs_inter_patient_mAP_df,
    c(
        treatment, mean_average_precision
    )
)
# rename mean_average_precision to inter_patient_mAP
colnames(organoid_fs_inter_patient_mAP_df)[2] <- "inter_patient_mAP"
# rename mean_average_precision to intra_patient_mAP
# combine intra and inter patient mAP for organoid
# merge the two dataframes on treatment and shuffle
organoid_sc_fs_mAP_df <- merge(
    organoid_fs_intra_patient_mAP_df,
    organoid_fs_inter_patient_mAP_df,
    by = c("treatment"),
    all = TRUE
)

organoid_sc_fs_mAP_df$intra_to_inter_ratio <- (
    organoid_sc_fs_mAP_df$intra_patient_mAP /
    organoid_sc_fs_mAP_df$inter_patient_mAP
)
head(organoid_sc_fs_mAP_df)

width <- 8
height <- 6
options(repr.plot.width = width, repr.plot.height = height)
ratio_plot <- (
    ggplot(
        data = organoid_sc_fs_mAP_df,
        aes(
            x = treatment,
            y = intra_to_inter_ratio,
            fill = treatment
        )
    )
    + geom_bar(stat = "identity", position = "dodge")
    + labs(
        x = "Treatment",
        y = "Intra to Inter Patient mAP Ratio",
        title = "Intra to Inter Patient mAP Ratio for Organoid and Single Cell FS"
    )
    + theme_bw()
    + theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
    )
    + geom_hline(
        yintercept = 1,
        linetype = "dashed",
        color = "red"
    )  # Adjust threshold as needed

)
ggsave(
    filename = file.path(figures_path, "intra_to_inter_patient_mAP_ratio_bar.png"),
    plot = ratio_plot,
    width = width,
    height = height,
    dpi = 600
)
ratio_plot

width <- 10
height <- 10
options(repr.plot.width = width, repr.plot.height = height)
plot <- (
    ggplot(
        data = organoid_sc_fs_mAP_df,
        aes(
            x = inter_patient_mAP,
            y = intra_patient_mAP,
            color = treatment,
            label = treatment
                    )
    )
    + geom_point(size = 2, alpha = 0.7)
    # + geom_boxplot(
    #     aes(fill = treatment),
    #     alpha = 0.4,
    #     size = 5,
    #     outlier.shape = NA,
    #     jitter = TRUE,
    #     position = position_dodge(width = 0.75)
    # )
    + labs(
        x = "Inter-patient mAP",
        y = "Intra-patient mAP",
        title = "Inter vs Intra Patient mAP for Organoid FS"
    )
    # + ylim(0, 1)
    # + xlim(0, 1)
    + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")  # Diagonal line
    + theme_bw()
    + theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 16)
    )
        + geom_text_repel(
        aes(
            label = treatment
        ),
        size = 4,
        nudge_y = 0.1,
        show.legend = FALSE,
        segment.color = "grey50",
        segment.size = 0.2,
        box.padding = 0.5,
        point.padding = 0.5,
        max.overlaps = 10
    )
    + xlim(0, 1)
)
ggsave(
    filename = file.path(figures_path, "organoid_fs_mAP_inter_vs_intra.png"),
    plot = plot,
    width = 10,
    height = 6,
    dpi = 600
)
plot

# Merge 2D and 3D for plotting
# Intra-patient
mAP_organoid_intra_merged <- inner_join(
    mAP_2d_organoid_intra %>% select(Metadata_treatment, Metadata_patient, mean_average_precision) %>% rename(mAP_2D = mean_average_precision),
    mAP_3d_organoid_intra %>% select(Metadata_treatment, Metadata_patient, mean_average_precision) %>% rename(mAP_3D = mean_average_precision),
    by = c("Metadata_treatment", "Metadata_patient")
)

mAP_sc_intra_merged <- inner_join(
    mAP_2d_sc_intra %>% select(Metadata_treatment, Metadata_patient, mean_average_precision) %>% rename(mAP_2D = mean_average_precision),
    mAP_3d_sc_intra %>% select(Metadata_treatment, Metadata_patient, mean_average_precision) %>% rename(mAP_3D = mean_average_precision),
    by = c("Metadata_treatment", "Metadata_patient")
)

# Inter-patient
mAP_organoid_inter_merged <- inner_join(
    mAP_2d_organoid_inter %>% select(Metadata_treatment_dose, mean_average_precision) %>% rename(mAP_2D = mean_average_precision),
    mAP_3d_organoid_inter %>% select(Metadata_treatment_dose, mean_average_precision) %>% rename(mAP_3D = mean_average_precision),
    by = "Metadata_treatment_dose"
)

mAP_sc_inter_merged <- inner_join(
    mAP_2d_sc_inter %>% select(Metadata_treatment_dose, mean_average_precision) %>% rename(mAP_2D = mean_average_precision),
    mAP_3d_sc_inter %>% select(Metadata_treatment_dose, mean_average_precision) %>% rename(mAP_3D = mean_average_precision),
    by = "Metadata_treatment_dose"
)

# Shared plot theme — defined once so all plots are consistent
plot_theme <- theme(
    plot.title   = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text  = element_text(size = 12)
)

# Organoid: 2D vs 3D Intra-patient mAP scatter plot faceted by patient
width <- 12
height <- 8
options(repr.plot.width = width, repr.plot.height = height)

organoid_2d_vs_3d_plot <- (
    ggplot(
        data = mAP_organoid_intra_merged,
        aes(
            x = mAP_2D,
            y = mAP_3D,
            color = Metadata_treatment
        )
    )
    + geom_point(size = 3, alpha = 0.7)
    + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")
    + labs(
        x = "2D mAP Score",
        y = "3D mAP Score",
        title = "2D vs 3D Intra-patient mAP — Organoid FS"
    )
    + theme_bw()
    + plot_theme
    + xlim(0, 1)
    + ylim(0, 1)
    + facet_wrap(~Metadata_patient, ncol = 4)
    + geom_text_repel(
        aes(label = Metadata_treatment),
        size = 3,
        show.legend = FALSE,
        segment.color = "grey50",
        segment.size = 0.2,
        box.padding = 0.5,
        point.padding = 0.5,
        max.overlaps = 10
    )
)
ggsave(filename = file.path(figures_path, "2d_vs_3d_organoid_intra_mAP.png"),
       plot = organoid_2d_vs_3d_plot, width = width, height = height, dpi = 600)
organoid_2d_vs_3d_plot

# Single cell: 2D vs 3D Intra-patient mAP scatter plot faceted by patient
width <- 12
height <- 8
options(repr.plot.width = width, repr.plot.height = height)

sc_2d_vs_3d_plot <- (
    ggplot(
        data = mAP_sc_intra_merged,
        aes(
            x = mAP_2D,
            y = mAP_3D,
            color = Metadata_treatment
        )
    )
    + geom_point(size = 3, alpha = 0.7)
    + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")
    + labs(
        x = "2D mAP Score",
        y = "3D mAP Score",
        title = "2D vs 3D Intra-patient mAP — Single Cell FS"
    )
    + theme_bw()
    + plot_theme
    + xlim(0, 1)
    + ylim(0, 1)
    + facet_wrap(~Metadata_patient, ncol = 4)
    + geom_text_repel(
        aes(label = Metadata_treatment),
        size = 3,
        show.legend = FALSE,
        segment.color = "grey50",
        segment.size = 0.2,
        box.padding = 0.5,
        point.padding = 0.5,
        max.overlaps = 10
    )
)
ggsave(filename = file.path(figures_path, "2d_vs_3d_sc_intra_mAP.png"),
       plot = sc_2d_vs_3d_plot, width = width, height = height, dpi = 600)
sc_2d_vs_3d_plot

# Organoid: 2D vs 3D Inter-patient mAP scatter plot
width <- 10
height <- 8
options(repr.plot.width = width, repr.plot.height = height)

plot_org_inter_mAP <- (
    ggplot(
        data = mAP_organoid_inter_merged,
        aes(x = mAP_2D, y = mAP_3D, color = Metadata_treatment_dose)
    )
    + geom_point(size = 4, alpha = 0.7)
    + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")
    + labs(
        x = "2D mAP Score",
        y = "3D mAP Score",
        title = "2D vs 3D Inter-patient mAP — Organoid Agg"
    )
    + theme_bw()
    + plot_theme
    + xlim(0, 1)
    + ylim(0, 1)
    + geom_text_repel(
        aes(label = Metadata_treatment_dose),
        size = 4, show.legend = FALSE,
        segment.color = "grey50", segment.size = 0.2,
        box.padding = 0.5, point.padding = 0.5, max.overlaps = 10
    )
)
ggsave(filename = file.path(figures_path, "2d_vs_3d_organoid_inter_mAP.png"),
       plot = plot_org_inter_mAP, width = width, height = height, dpi = 600)
plot_org_inter_mAP

# Single cell: 2D vs 3D Inter-patient mAP scatter plot
width <- 10
height <- 8
options(repr.plot.width = width, repr.plot.height = height)

plot_sc_inter_mAP <- (
    ggplot(
        data = mAP_sc_inter_merged,
        aes(x = mAP_2D, y = mAP_3D, color = Metadata_treatment_dose)
    )
    + geom_point(size = 4, alpha = 0.7)
    + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")
    + labs(
        x = "2D mAP Score",
        y = "3D mAP Score",
        title = "2D vs 3D Inter-patient mAP — Single Cell Agg"
    )
    + theme_bw()
    + plot_theme
    + xlim(0, 1)
    + ylim(0, 1)
    + geom_text_repel(
        aes(label = Metadata_treatment_dose),
        size = 4, show.legend = FALSE,
        segment.color = "grey50", segment.size = 0.2,
        box.padding = 0.5, point.padding = 0.5, max.overlaps = 10
    )
)
ggsave(filename = file.path(figures_path, "2d_vs_3d_sc_inter_mAP.png"),
       plot = plot_sc_inter_mAP, width = width, height = height, dpi = 600)
plot_sc_inter_mAP

#Load Data
dist_2d_organoid_intra <- arrow::read_parquet(file.path(root_dir, "1.EDA/results/distance_metrics/2d_organoid_intra_patient_cosine_distance.parquet"))
dist_3d_organoid_intra <- arrow::read_parquet(file.path(root_dir, "1.EDA/results/distance_metrics/3d_organoid_intra_patient_cosine_distance.parquet"))
dist_2d_organoid_inter <- arrow::read_parquet(file.path(root_dir, "1.EDA/results/distance_metrics/2d_organoid_inter_patient_cosine_distance.parquet"))
dist_3d_organoid_inter <- arrow::read_parquet(file.path(root_dir, "1.EDA/results/distance_metrics/3d_organoid_inter_patient_cosine_distance.parquet"))
dist_2d_sc_intra <- arrow::read_parquet(file.path(root_dir, "1.EDA/results/distance_metrics/2d_sc_intra_patient_cosine_distance.parquet"))
dist_3d_sc_intra <- arrow::read_parquet(file.path(root_dir, "1.EDA/results/distance_metrics/3d_sc_intra_patient_cosine_distance.parquet"))
dist_2d_sc_inter <- arrow::read_parquet(file.path(root_dir, "1.EDA/results/distance_metrics/2d_sc_inter_patient_cosine_distance.parquet"))
dist_3d_sc_inter <- arrow::read_parquet(file.path(root_dir, "1.EDA/results/distance_metrics/3d_sc_inter_patient_cosine_distance.parquet"))

# Distance Metrics output path
dist_figures_path <- file.path(root_dir, "1.EDA/figures/distance_metrics")
if (!dir.exists(dist_figures_path)) {
    dir.create(dist_figures_path, recursive = TRUE)
}

# Distance Metrics: Merge 2D and 3D for plotting
# Intra-patient
dist_organoid_intra_merged <- inner_join(
    dist_2d_organoid_intra %>% select(Metadata_treatment_dose, Metadata_patient, cosine_distance_mean) %>% rename(cosine_2D = cosine_distance_mean),
    dist_3d_organoid_intra %>% select(Metadata_treatment_dose, Metadata_patient, cosine_distance_mean) %>% rename(cosine_3D = cosine_distance_mean),
    by = c("Metadata_treatment_dose", "Metadata_patient")
)

dist_sc_intra_merged <- inner_join(
    dist_2d_sc_intra %>% select(Metadata_treatment_dose, Metadata_patient, cosine_distance_mean) %>% rename(cosine_2D = cosine_distance_mean),
    dist_3d_sc_intra %>% select(Metadata_treatment_dose, Metadata_patient, cosine_distance_mean) %>% rename(cosine_3D = cosine_distance_mean),
    by = c("Metadata_treatment_dose", "Metadata_patient")
)

# Inter-patient
dist_organoid_inter_merged <- inner_join(
    dist_2d_organoid_inter %>% select(Metadata_treatment_dose, cosine_distance_mean) %>% rename(cosine_2D = cosine_distance_mean),
    dist_3d_organoid_inter %>% select(Metadata_treatment_dose, cosine_distance_mean) %>% rename(cosine_3D = cosine_distance_mean),
    by = "Metadata_treatment_dose"
)

dist_sc_inter_merged <- inner_join(
    dist_2d_sc_inter %>% select(Metadata_treatment_dose, cosine_distance_mean) %>% rename(cosine_2D = cosine_distance_mean),
    dist_3d_sc_inter %>% select(Metadata_treatment_dose, cosine_distance_mean) %>% rename(cosine_3D = cosine_distance_mean),
    by = "Metadata_treatment_dose"
)

# Distance Metrics: 2D vs 3D Intra-patient — Organoid (faceted by patient)
width <- 12
height <- 8
options(repr.plot.width = width, repr.plot.height = height)

plot_org_intra <- (
    ggplot(
        data = dist_organoid_intra_merged,
        aes(x = cosine_2D, y = cosine_3D, color = Metadata_treatment_dose)
    )
    + geom_point(size = 3, alpha = 0.7)
    + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")
    + labs(
        x = "2D Cosine Distance",
        y = "3D Cosine Distance",
        title = "2D vs 3D Intra-patient Cosine Distance — Organoid Agg"
    )
    + theme_bw()
    + plot_theme
    + facet_wrap(~Metadata_patient, ncol = 4)
    + geom_text_repel(
        aes(label = Metadata_treatment_dose),
        size = 3, show.legend = FALSE,
        segment.color = "grey50", segment.size = 0.2,
        box.padding = 0.5, point.padding = 0.5, max.overlaps = 10
    )
)
ggsave(filename = file.path(dist_figures_path, "2d_vs_3d_organoid_intra_cosine_distance.png"),
       plot = plot_org_intra, width = width, height = height, dpi = 600)
plot_org_intra

# Distance Metrics: 2D vs 3D Intra-patient — Single Cell (faceted by patient)
width <- 12
height <- 8
options(repr.plot.width = width, repr.plot.height = height)

plot_sc_intra <- (
    ggplot(
        data = dist_sc_intra_merged,
        aes(x = cosine_2D, y = cosine_3D, color = Metadata_treatment_dose)
    )
    + geom_point(size = 3, alpha = 0.7)
    + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")
    + labs(
        x = "2D Cosine Distance",
        y = "3D Cosine Distance",
        title = "2D vs 3D Intra-patient Cosine Distance — Single Cell Agg"
    )
    + theme_bw()
    + plot_theme
    + facet_wrap(~Metadata_patient, ncol = 4)
    + geom_text_repel(
        aes(label = Metadata_treatment_dose),
        size = 3, show.legend = FALSE,
        segment.color = "grey50", segment.size = 0.2,
        box.padding = 0.5, point.padding = 0.5, max.overlaps = 10
    )
)
ggsave(filename = file.path(dist_figures_path, "2d_vs_3d_sc_intra_cosine_distance.png"),
       plot = plot_sc_intra, width = width, height = height, dpi = 600)
plot_sc_intra

# Distance Metrics: 2D vs 3D Inter-patient — Organoid
width <- 10
height <- 8
options(repr.plot.width = width, repr.plot.height = height)

plot_org_inter <- (
    ggplot(
        data = dist_organoid_inter_merged,
        aes(x = cosine_2D, y = cosine_3D, color = Metadata_treatment_dose)
    )
    + geom_point(size = 4, alpha = 0.7)
    + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")
    + labs(
        x = "2D Cosine Distance",
        y = "3D Cosine Distance",
        title = "2D vs 3D Inter-patient Cosine Distance — Organoid Agg"
    )
    + theme_bw()
    + plot_theme
    + geom_text_repel(
        aes(label = Metadata_treatment_dose),
        size = 4, show.legend = FALSE,
        segment.color = "grey50", segment.size = 0.2,
        box.padding = 0.5, point.padding = 0.5, max.overlaps = 10
    )
)
ggsave(filename = file.path(dist_figures_path, "2d_vs_3d_organoid_inter_cosine_distance.png"),
       plot = plot_org_inter, width = width, height = height, dpi = 600)
plot_org_inter

# Distance Metrics: 2D vs 3D Inter-patient — Single Cell
width <- 10
height <- 8
options(repr.plot.width = width, repr.plot.height = height)

plot_sc_inter <- (
    ggplot(
        data = dist_sc_inter_merged,
        aes(x = cosine_2D, y = cosine_3D, color = Metadata_treatment_dose)
    )
    + geom_point(size = 4, alpha = 0.7)
    + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")
    + labs(
        x = "2D Cosine Distance",
        y = "3D Cosine Distance",
        title = "2D vs 3D Inter-patient Cosine Distance — Single Cell Agg"
    )
    + theme_bw()
    + plot_theme
    + geom_text_repel(
        aes(label = Metadata_treatment_dose),
        size = 4, show.legend = FALSE,
        segment.color = "grey50", segment.size = 0.2,
        box.padding = 0.5, point.padding = 0.5, max.overlaps = 10
    )
)
ggsave(filename = file.path(dist_figures_path, "2d_vs_3d_sc_inter_cosine_distance.png"),
       plot = plot_sc_inter, width = width, height = height, dpi = 600)
plot_sc_inter
