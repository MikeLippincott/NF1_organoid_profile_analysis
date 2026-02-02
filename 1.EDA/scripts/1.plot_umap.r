list_of_packages <- c("ggplot2", "dplyr", "tidyr", "circlize")
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
#

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

figures_path <- file.path(root_dir,"1.EDA/figures/umaps")
if (!dir.exists(figures_path)) {
  dir.create(figures_path, recursive = TRUE)
}

umap_theme <- theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size = 12)
    )

#SETTING UP INPUTS FOR GRAPHS

#organoid_umap_results <- arrow::read_parquet(file.path(root_dir,"1.EDA/results/umap/organoid_fs_umap.parquet"))
#head(organoid_umap_results)

max_projection_2D_sc_umap_results <- arrow::read_parquet(file.path(root_dir,"1.EDA/results/umap/2D/max_projection/sc_fs_umap.parquet"))
head(max_projection_2D_sc_umap_results)
#Combining STAURO and Staurosporine into one label
max_projection_2D_sc_umap_results <- max_projection_2D_sc_umap_results %>%
    mutate(Metadata_treatment = ifelse(Metadata_treatment == "STAURO", "Staurosporine", Metadata_treatment))

sc_3D_umap_results <- arrow::read_parquet(file.path(root_dir,"1.EDA/results/umap/3D/sc_fs_umap.parquet"))
head(sc_3D_umap_results)
sc_3D_umap_results <- sc_3D_umap_results %>%
    mutate(Metadata_treatment = ifelse(Metadata_treatment == "STAURO", "Staurosporine", Metadata_treatment))


middle_slice_2D_sc_umap_results <- arrow::read_parquet(file.path(root_dir,"1.EDA/results/umap/2D/middle_slice/sc_fs_umap.parquet"))
head(middle_slice_2D_sc_umap_results)
middle_slice_2D_sc_umap_results <- middle_slice_2D_sc_umap_results %>%
    mutate(Metadata_treatment = ifelse(Metadata_treatment == "STAURO", "Staurosporine", Metadata_treatment))

# set custom colors for each MOA
custom_MOA_palette <- c(
    'Control' = "#5a5c5d",
    'MEK1/2 inhibitor' = "#882E8B",


    'HDAC inhibitor' = "#1E6B61",
    'PI3K and HDAC inhibitor' = "#2E6B8B",
    'PI3K inhibitor'="#0092E0",

    'receptor tyrosine kinase inhibitor' = "#576A20",
    'tyrosine kinase inhibitor' = "#646722",

    'mTOR inhibitor' = "#ACE089",
    'IGF-1R inhibitor' = "#ACE040",

    'HSP90 inhibitor'="#33206A",
    'Apoptosis induction'="#272267",
    'Na+/K+ pump inhibitor' = "#A16C28",
    'histamine H1 receptor antagonist' = "#3A8F00",
    'DNA binding' = "#174F17",
    'BRD4 inhibitor' = "#ff0000"

)


# Set custom colors for each treatment
custom_treatment_palette <- c(
    'DMSO' = "#5a5c5d",              # Control - gray
    'Staurosporine' = "#7D2780",     # Dark purple

    'Fimepinostat' = "#1E6B61",      # Teal (HDAC inhibitor)
    'Copanlisib' = "#0092E0",        # Blue (PI3K inhibitor)

    'Imatinib' = "#576A20",          # Olive green
    'Nilotinib' = "#646722",         # Yellow-green
    'Cabozantinib' = "#758B2D",      # Light olive

    'Everolimus' = "#ACE089",        # Light green (mTOR inhibitor)
    'Rapamycin' = "#90D070",         # Medium green (mTOR inhibitor)
    'Linsitinib' = "#ACE040",        # Yellow-green (IGF-1R inhibitor)

    'Onalespib' = "#33206A",         # Dark purple (HSP90 inhibitor)
    'Digoxin' = "#A16C28",           # Orange-brown
    'Ketotifen' = "#3A8F00",         # Green

    'Binimetinib' = "#ff0000",       # Red (MEK inhibitor)
    'Mirdametinib' = "#cc0000",      # Dark red (MEK inhibitor)
    'Trametinib' = "#ff3333",        # Light red (MEK inhibitor)
    'Selumetinib' = "#ff6666"        # Lighter red (MEK inhibitor)
)

# head(organoid_umap_results)

# width <- 10
# height <- 5
# options(repr.plot.width = width, repr.plot.height = height)
# umap_organoid_plot <- (
#     ggplot(organoid_umap_results, aes(x = UMAP1, y = UMAP2, color = Target, size = single_cell_count))
#     + geom_point(alpha = 0.5)
#     + scale_color_manual(values = custom_MOA_palette)
#     + labs(title = "All patients: Organoid FS Profiles", x = "UMAP 0", y = "UMAP 1")
#     + theme_bw()
#     + umap_theme

#     + guides(
#         size = guide_legend(
#             title = "Single Cell Count",
#             text = element_text(size = 16, hjust = 0.5, position = "top"),
#             nrow = 1,
#             ),
#         color = guide_legend(
#             title = "Target",
#             text = element_text(size = 16, hjust = 0.5),
#             override.aes = list(alpha = 1,size = 5),
#             ncol = 1
#         )
#     )
#     + facet_wrap(~patient, nrow = 2)
# )
# organoid_features_path <- file.path(figures_path, "all_patients_umap_organoid_features_facet_by_patient.png")
# ggsave(umap_organoid_plot, file = organoid_features_path, width = width, height = height, dpi = 300)
# umap_organoid_plot

# patients <- unique(organoid_umap_results$patient)
# hex_codes <- c(
# "#86C436",
# "#BFD468",
# "#36C4BB",
# "#68D4B4",
# "#7336C4",
# "#7E68D4",
# "#C4363F",
# "#D46888"
# )
# patient_color_palette <- setNames(hex_codes[1:length(patients)], patients)

# width <- 10
# height <- 5
# options(repr.plot.width = width, repr.plot.height = height)
# umap_organoid_plot <- (
#     ggplot(organoid_umap_results, aes(x = UMAP1, y = UMAP2, color = patient, size = single_cell_count))
#     + geom_point(alpha = 0.7)
#     + scale_color_manual(values = patient_color_palette)
#     + labs(title = "All patients: Organoid FS Profiles", x = "UMAP 0", y = "UMAP 1")
#     + theme_bw()
#     + umap_theme
#     + guides(
#         size = guide_legend(
#             title = "Single Cell Count",
#             text = element_text(size = 16, hjust = 0.5, position = "top"),
#             nrow = 1,
#             ),
#         color = guide_legend(
#             title = "Patient",
#             text = element_text(size = 16, hjust = 0.5),
#             override.aes = list(alpha = 1,size = 5),
#             ncol = 1
#         )
#     )
#     + facet_wrap(~Target, nrow = 4)
# )
# organoid_features_path <- file.path(figures_path, "all_patients_umap_organoid_features_facet_by_target.png")

# ggsave(umap_organoid_plot, file = organoid_features_path, width = width, height = height, dpi = 300)
# umap_organoid_plot

# get just the DMSO treatments
# organoid_umap_results_dmsos <- organoid_umap_results %>% filter(Target == "Control")

# width <- 10
# height <- 5
# options(repr.plot.width = width, repr.plot.height = height)
# umap_organoid_plot <- (
#     ggplot(organoid_umap_results_dmsos, aes(x = UMAP1, y = UMAP2, color = patient, size = single_cell_count))
#     + geom_point(alpha = 0.7)
#     + scale_color_manual(values = patient_color_palette)
#     + labs(title = "All patients: Organoid FS Profiles", x = "UMAP 0", y = "UMAP 1")
#     + theme_bw()
#     + umap_theme

#     + guides(
#         size = guide_legend(
#             title = "Single Cell Count",
#             text = element_text(size = 16, hjust = 0.5, position = "top"),
#             nrow = 1,
#             ),
#         color = guide_legend(
#             title = "Patients",
#             text = element_text(size = 16, hjust = 0.5),
#             override.aes = list(alpha = 1,size = 5),
#             ncol = 1
#         )
#     )
# )
# organoid_features_path <- file.path(figures_path, "all_patients_controls_only_umap_organoid_features.png")

# ggsave(umap_organoid_plot, file = organoid_features_path, width = width, height = height, dpi = 300)
# umap_organoid_plot

# Define all possible patients upfront
all_patients <- c("NF0014", "NF0016", "NF0018", "NF0021", "NF0030",  "NF0040", "SARCO219", "SARCO361")

hex_codes <- c(
    "#86C436",
    "#BFD468",
    "#36C4BB",
    "#68D4B4",
    "#7336C4",
    "#7E68D4",
    "#C4363F",
    "#D46888"
)

# Create master palette
master_patient_palette <- setNames(hex_codes[1:length(all_patients)], all_patients)

patient_color_palette <- master_patient_palette[names(master_patient_palette) %in% unique(max_projection_2D_sc_umap_results$Metadata_patient)]

width <- 10
height <- 5
options(repr.plot.width = width, repr.plot.height = height)
umap_sc_plot <- (
    ggplot(max_projection_2D_sc_umap_results, aes(x = UMAP1, y = UMAP2, color = Metadata_treatment))
    + geom_point(alpha = 0.5, size = 1)
    + scale_color_manual(values = custom_treatment_palette)
    + labs(title = "All patients: 2D MIP Single Cell FS Profiles", x = "UMAP 0", y = "UMAP 1")
    + theme_bw()
    + umap_theme

    + guides(
        color = guide_legend(
            title = "Treatment",
            override.aes = list(alpha = 1,size = 5),
            ncol = 1
        )
    )
    + facet_wrap(~Metadata_patient, nrow = 2)
)
sc_2D_MIP_facet_by_features_path <- file.path(figures_path, "2D", "max_projection", "all_patients_umap_2D_max_intensity_sc_features_facet_by_patient.png")
ggsave(umap_sc_plot, file = sc_2D_MIP_facet_by_features_path, width = width, height = height, dpi = 300)
umap_sc_plot

width <- 10
height <- 5
options(repr.plot.width = width, repr.plot.height = height)
umap_sc_plot_by_patient <- (
    ggplot(max_projection_2D_sc_umap_results, aes(x = UMAP1, y = UMAP2, color = Metadata_patient))
    + geom_point(alpha = 0.4, size = 1)
    + scale_color_manual(values = patient_color_palette)  # You'll need to define this palette
    + labs(title = "All patients: 2D MIP Single Cell FS Profiles", x = "UMAP 0", y = "UMAP 1")
    + theme_bw()
    + umap_theme

    + guides(
        color = guide_legend(
            title = "Patient",
            override.aes = list(alpha = 1, size = 5),
            ncol = 1
        )
    )
)
sc_features_2D_MIP_by_patient_path <- file.path(figures_path, "2D", "max_projection", "all_patients_umap_2D_max_intensity_sc_features_by_patient.png")
ggsave(umap_sc_plot_by_patient, file = sc_features_2D_MIP_by_patient_path, width = width, height = height, dpi = 300)
umap_sc_plot_by_patient

width <- 10
height <- 5
options(repr.plot.width = width, repr.plot.height = height)
umap_sc_plot <- (
    ggplot(max_projection_2D_sc_umap_results, aes(x = UMAP1, y = UMAP2, color = Metadata_treatment))
    + geom_point(alpha = 0.4, size = 1)
    + scale_color_manual(values = custom_treatment_palette)
    + labs(title = "All patients: 2D MIP Single Cell FS Profiles", x = "UMAP 0", y = "UMAP 1")
    + theme_bw()
    + umap_theme

    + guides(
        color = guide_legend(
            title = "Treatment",
            override.aes = list(alpha = 1, size = 5),
            ncol = 1
        )
    )
)
sc_features_2D_MIP_by_treatment_path <- file.path(figures_path,"2D", "max_projection", "all_patients_umap_2D_max_intensity_sc_features_by_treatment.png")
ggsave(umap_sc_plot, file = sc_features_2D_MIP_by_treatment_path, width = width, height = height, dpi = 300)
umap_sc_plot

width <- 10
height <- 8
options(repr.plot.width = width, repr.plot.height = height)
umap_sc_plot_facet_treatment <- (
    ggplot(max_projection_2D_sc_umap_results, aes(x = UMAP1, y = UMAP2, color = Metadata_patient))
    + geom_point(alpha = 0.4, size = 1)
    + scale_color_manual(values = patient_color_palette)
    + labs(title = "All patients: 2D MIP Single Cell FS Profiles", x = "UMAP 0", y = "UMAP 1")
    + theme_bw()
    + umap_theme

    + guides(
        color = guide_legend(
            override.aes = list(alpha = 1, size = 5),
            ncol = 1
        )
    )
    + facet_wrap(~Metadata_treatment, nrow = 4)  # Changed from Target to Metadata_treatment
)
sc_features_2D_MIP_facet_by_treatment_path <- file.path(figures_path,"2D", "max_projection", "all_patients_umap_2D_max_intensity_sc_features_facet_by_treatment.png")
ggsave(umap_sc_plot_facet_treatment, file = sc_features_2D_MIP_facet_by_treatment_path, width = width, height = height, dpi = 300)
umap_sc_plot_facet_treatment

# Filter for DMSO only
sc_umap_results_dmsos <- max_projection_2D_sc_umap_results %>%
    filter(Metadata_treatment == "DMSO")

width <- 10
height <- 5
options(repr.plot.width = width, repr.plot.height = height)
umap_sc_plot_dmso <- (
    ggplot(sc_umap_results_dmsos, aes(x = UMAP1, y = UMAP2, color = Metadata_patient))
    + geom_point(alpha = 0.7, size = 1)
    + scale_color_manual(values = patient_color_palette)
    + labs(title = "All patients: 2D MIP Single Cell FS Profiles (DMSO Controls)", x = "UMAP 0", y = "UMAP 1")
    + theme_bw()
    + umap_theme
    + theme(
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, hjust = 0.5)
    )
    + guides(
        color = guide_legend(
            title = "Patient",
            override.aes = list(alpha = 1, size = 5),
            ncol = 1
        )
    )
)
sc_features_2D_MIP_dmso_path <- file.path(figures_path, "2D", "max_projection","all_patient_umap_2D_max_intensity_sc_fs_dmso_only.png")
ggsave(umap_sc_plot_dmso, file = sc_features_2D_MIP_dmso_path, width = width, height = height, dpi = 300)
umap_sc_plot_dmso

width <- 10
height <- 5
options(repr.plot.width = width, repr.plot.height = height)
umap_sc_plot <- (
    ggplot(middle_slice_2D_sc_umap_results, aes(x = UMAP1, y = UMAP2, color = Metadata_treatment))
    + geom_point(alpha = 0.4, size = 1)
    + scale_color_manual(values = custom_treatment_palette)
    + labs(title = "All patients: 2D Middle Slice Single Cell FS Profiles", x = "UMAP 0", y = "UMAP 1")
    + theme_bw()
    + umap_theme

    + guides(
        color = guide_legend(
            title = "Treatment",
            override.aes = list(alpha = 1,size = 5),
            ncol = 1
        )
    )
    + facet_wrap(~Metadata_patient, nrow = 2)
)
sc_features_2D_MS_facet_by_patient_path <- file.path(figures_path, "2D", "middle_slice","all_patients_umap_2D_middle_slice_sc_features_facet_by_patient.png")
ggsave(umap_sc_plot, file = sc_features_2D_MS_facet_by_patient_path, width = width, height = height, dpi = 300)
umap_sc_plot

width <- 10
height <- 5
options(repr.plot.width = width, repr.plot.height = height)
umap_sc_plot_by_patient <- (
    ggplot(middle_slice_2D_sc_umap_results, aes(x = UMAP1, y = UMAP2, color = Metadata_patient))
    + geom_point(alpha = 0.4, size = 1)
    + scale_color_manual(values = patient_color_palette)  # You'll need to define this palette
    + labs(title = "All patients: 2D Middle Slice Single Cell FS Profiles", x = "UMAP 0", y = "UMAP 1")
    + theme_bw()
    + umap_theme

    + guides(
        color = guide_legend(
            title = "Patient",
            override.aes = list(alpha = 1, size = 5),
            ncol = 1
        )
    )
)
sc_features_2D_MS_by_patient_path <- file.path(figures_path, "2D", "middle_slice", "all_patients_umap_2D_middle_slice_sc_features_by_patient.png")
ggsave(umap_sc_plot_by_patient, file = sc_features_2D_MS_by_patient_path, width = width, height = height, dpi = 300)
umap_sc_plot_by_patient

width <- 10
height <- 5
options(repr.plot.width = width, repr.plot.height = height)
umap_sc_plot_by_treatment <- (
    ggplot(middle_slice_2D_sc_umap_results, aes(x = UMAP1, y = UMAP2, color = Metadata_treatment))
    + geom_point(alpha = 0.4, size = 1)
    + scale_color_manual(values = custom_treatment_palette)  # You'll need to define this palette
    + labs(title = "All patients: 2D Middle Slice Single Cell FS Profiles", x = "UMAP 0", y = "UMAP 1")
    + theme_bw()
    + umap_theme

    + guides(
        color = guide_legend(
            title = "Treatment",
            override.aes = list(alpha = 1, size = 5),
            ncol = 1
        )
    )
)
sc_features_2D_MS_by_treatment_path <- file.path(figures_path, "2D", "middle_slice", "all_patients_umap_2D_middle_slice_sc_features_by_treatment.png")
ggsave(umap_sc_plot_by_treatment, file = sc_features_2D_MS_by_treatment_path, width = width, height = height, dpi = 300)
umap_sc_plot_by_treatment

width <- 10
height <- 8
options(repr.plot.width = width, repr.plot.height = height)
umap_sc_plot_MS_facet_treatment <- (
    ggplot(middle_slice_2D_sc_umap_results, aes(x = UMAP1, y = UMAP2, color = Metadata_patient))
    + geom_point(alpha = 0.4, size = 1)
    + scale_color_manual(values = patient_color_palette)
    + labs(title = "All patients: 2D Middle Slice Single Cell FS Profiles", x = "UMAP 0", y = "UMAP 1")
    + theme_bw()
    + umap_theme

    + guides(
        color = guide_legend(
            override.aes = list(alpha = 1, size = 5),
            ncol = 1
        )
    )
    + facet_wrap(~Metadata_treatment, nrow = 4)
)
sc_features_2D_MS_facet_by_treatment_path <- file.path(figures_path, "2D", "middle_slice", "all_patients_umap_2D_middle_slice_sc_features_facet_by_treatment.png")
ggsave(umap_sc_plot_MS_facet_treatment, file = sc_features_2D_MS_facet_by_treatment_path, width = width, height = height, dpi = 300)
umap_sc_plot_MS_facet_treatment

# Filter for DMSO only
sc_umap_results_dmsos <- middle_slice_2D_sc_umap_results %>%
    filter(Metadata_treatment == "DMSO")

width <- 10
height <- 5
options(repr.plot.width = width, repr.plot.height = height)
umap_sc_plot_2D_MS_dmso <- (
    ggplot(sc_umap_results_dmsos, aes(x = UMAP1, y = UMAP2, color = Metadata_patient))
    + geom_point(alpha = 0.7, size = 1)
    + scale_color_manual(values = patient_color_palette)
    + labs(title = "All patients: 2D MIP Single Cell FS Profiles (DMSO Controls)", x = "UMAP 0", y = "UMAP 1")
    + theme_bw()
    + umap_theme
    + theme(
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, hjust = 0.5)
    )
    + guides(
        color = guide_legend(
            title = "Patient",
            override.aes = list(alpha = 1, size = 5),
            ncol = 1
        )
    )
)
sc_features_2D_MS_dmso_path <- file.path(figures_path, "2D", "middle_slice", "all_patients_umap_2D_middle_slice_sc_fs_dmso_only.png")
ggsave(umap_sc_plot_2D_MS_dmso, file = sc_features_2D_MS_dmso_path, width = width, height = height, dpi = 300)
umap_sc_plot_2D_MS_dmso

width <- 10
height <- 5
options(repr.plot.width = width, repr.plot.height = height)
umap_sc_plot <- (
    ggplot(sc_3D_umap_results, aes(x = UMAP1, y = UMAP2, color = Metadata_treatment))
    + geom_point(alpha = 0.5, size = 1)
    + scale_color_manual(values = custom_treatment_palette)
    + labs(title = "All patients: 3D Single Cell FS Profiles", x = "UMAP 0", y = "UMAP 1")
    + theme_bw()
    + umap_theme

    + guides(
        color = guide_legend(
            title = "Treatment",
            override.aes = list(alpha = 1,size = 5),
            ncol = 1
        )
    )
    + facet_wrap(~Metadata_patient, nrow = 2)
)
sc_3D_features_by_treatment_path <- file.path(figures_path, "3D", "all_patients_umap_3D_sc_features_facet_by_patient.png")
ggsave(umap_sc_plot, file = sc_3D_features_by_treatment_path, width = width, height = height, dpi = 300)
umap_sc_plot

width <- 10
height <- 5
options(repr.plot.width = width, repr.plot.height = height)
umap_sc_plot_by_patient <- (
    ggplot(sc_3D_umap_results, aes(x = UMAP1, y = UMAP2, color = Metadata_patient))
    + geom_point(alpha = 0.7, size = 1)
    + scale_color_manual(values = patient_color_palette)  # You'll need to define this palette
    + labs(title = "All patients: 3D Single Cell FS Profiles", x = "UMAP 0", y = "UMAP 1")
    + theme_bw()
    + umap_theme

    + guides(
        color = guide_legend(
            title = "Patient",
            override.aes = list(alpha = 1, size = 5),
            ncol = 1
        )
    )
)
sc_features_3D_by_patient_path <- file.path(figures_path, "3D", "all_patients_umap_3D_features_by_patient.png")
ggsave(umap_sc_plot_by_patient, file = sc_features_3D_by_patient_path, width = width, height = height, dpi = 300)
umap_sc_plot_by_patient

width <- 10
height <- 5
options(repr.plot.width = width, repr.plot.height = height)
umap_sc_plot_by_treatment <- (
    ggplot(sc_3D_umap_results, aes(x = UMAP1, y = UMAP2, color = Metadata_treatment))
    + geom_point(alpha = 0.7, size = 1)
    + scale_color_manual(values = custom_treatment_palette)
    + labs(title = "All patients: 3D Single Cell FS Profiles", x = "UMAP 0", y = "UMAP 1")
    + theme_bw()
    + umap_theme

    + guides(
        color = guide_legend(
            title = "Treatment",
            override.aes = list(alpha = 1, size = 5),
            ncol = 1
        )
    )
)
sc_3D_features_by_treatment_path <- file.path(figures_path, "3D", "all_patients_umap_3D_features_by_treatment.png")
ggsave(umap_sc_plot_by_treatment, file = sc_3D_features_by_treatment_path, width = width, height = height, dpi = 300)
umap_sc_plot_by_treatment

width <- 10
height <- 8
options(repr.plot.width = width, repr.plot.height = height)
umap_sc_plot_facet_treatment <- (
    ggplot(sc_3D_umap_results, aes(x = UMAP1, y = UMAP2, color = Metadata_patient))
    + geom_point(alpha = 0.7, size = 1)
    + scale_color_manual(values = patient_color_palette)
    + labs(title = "All patients: 3D Single Cell FS Profiles", x = "UMAP 0", y = "UMAP 1")
    + theme_bw()
    + umap_theme

    + guides(
        color = guide_legend(
            override.aes = list(alpha = 1, size = 5),
            ncol = 1
        )
    )
    + facet_wrap(~Metadata_treatment, nrow = 4)
)
sc_3D_features_facet_treatment_path <- file.path(figures_path, "3D", "all_patients_umap_3D_sc_features_facet_by_treatment.png")
ggsave(umap_sc_plot_facet_treatment, file = sc_3D_features_facet_treatment_path, width = width, height = height, dpi = 300)
umap_sc_plot_facet_treatment

# Filter for DMSO only
sc_umap_results_dmsos <- sc_3D_umap_results %>%
    filter(Metadata_treatment == "DMSO")

width <- 10
height <- 5
options(repr.plot.width = width, repr.plot.height = height)
umap_sc_plot_dmso <- (
    ggplot(sc_umap_results_dmsos, aes(x = UMAP1, y = UMAP2, color = Metadata_patient))
    + geom_point(alpha = 0.7, size = 1)
    + scale_color_manual(values = patient_color_palette)
    + labs(title = "All patients: 3DSingle Cell FS Profiles (DMSO Controls)", x = "UMAP 0", y = "UMAP 1")
    + theme_bw()
    + umap_theme
    + theme(
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, hjust = 0.5)
    )
    + guides(
        color = guide_legend(
            title = "Patient",
            override.aes = list(alpha = 1, size = 5),
            ncol = 1
        )
    )
)
sc_features_3D_dmso_path <- file.path(figures_path, "3D", "all_patients_sc_3D_fs_by_patient_dmso_only.png")
ggsave(umap_sc_plot_dmso, file = sc_features_3D_dmso_path, width = width, height = height, dpi = 300)
umap_sc_plot_dmso

individual_patients <- c("NF0014_T1", "NF0016_T1", "NF0018_T6", "NF0021_T1", "NF0030_T1",  "NF0040_T1", "SARCO219_T2", "SARCO361_T1")
hex_codes <- c(
    "#86C436",
    "#BFD468",
    "#36C4BB",
    "#68D4B4",
    "#7336C4",
    "#7E68D4",
    "#C4363F",
    "#D46888"
)

for (patient in individual_patients) {
    # patient_umap_file_path <- file.path(root_dir, paste0("5.EDA/results/patient_results/",patient,"_organoid_fs_umap.parquet"))
    # umap_results_patient <- arrow::read_parquet(patient_umap_file_path)
    # umap_results_patient

    # width <- 10
    # height <- 5
    # options(repr.plot.width = width, repr.plot.height = height)
    # umap_organoid_plot <- (
    #     ggplot(umap_results_patient, aes(x = UMAP1, y = UMAP2, color = Target, size = single_cell_count))
    #     + geom_point(alpha = 0.7)
    #     + scale_color_manual(values = custom_MOA_palette)
    #     + labs(title = paste(patient," - Organoid FS Profiles"), x = "UMAP 0", y = "UMAP 1")
    #     + theme_bw()
    #     + umap_theme

    #     + guides(
    #         size = guide_legend(
    #             title = "Single Cell Count",
    #             text = element_text(size = 16, hjust = 0.5),
    #             # make two columns for the legend
    #             nrow = 1,
    #             # title position to top

    #             # move to bottom
    #             # position = "bottom"

    #             ),
    #         color = guide_legend(
    #             title = "Target",
    #             text = element_text(size = 16, hjust = 0.5, position = "top"),
    #             override.aes = list(alpha = 1,size = 5),
    #             # nrow = 2,
    #             ncol = 1
    #             # position = "bottom"
    #         )
    #     )

    # )
    # print(umap_organoid_plot)

    patient_umap_file_path <- file.path(root_dir, paste0("1.EDA/results/umap/patient_results/2D/max_projection/",patient,"_sc_fs_umap.parquet"))
    umap_results_patient <- arrow::read_parquet(patient_umap_file_path)

    width <- 10
    height <- 5
    options(repr.plot.width = width, repr.plot.height = height)

    umap_sc_2D_MIP_color_by_treatment_plot <- (
        ggplot(umap_results_patient, aes(x = UMAP1, y = UMAP2, color = Metadata_treatment))
        + geom_point(alpha = 0.6)
        + scale_color_manual(values = custom_treatment_palette)  # adjust if you need a different palette
        + labs(title = paste0(patient," - Single-cells MIP FS Profiles"), x = "UMAP 0", y = "UMAP 1")
        + theme_bw()
        + umap_theme
        + guides(
            color = guide_legend(
                title = "Treatment",
                override.aes = list(alpha = 1, size = 5),
                ncol = 1
            )
        )
    )

    print(umap_sc_2D_MIP_color_by_treatment_plot)

    patient_umap_file_path <- file.path(root_dir, paste0("1.EDA/results/umap/patient_results/2D/middle_slice/",patient,"_sc_fs_umap.parquet"))
    umap_results_patient <- arrow::read_parquet(patient_umap_file_path)

    width <- 10
    height <- 5
    options(repr.plot.width = width, repr.plot.height = height)

    umap_sc_2D_middle_slice_color_by_treatment_plot <- (
        ggplot(umap_results_patient, aes(x = UMAP1, y = UMAP2, color = Metadata_treatment))
        + geom_point(alpha = 0.6)
        + scale_color_manual(values = custom_treatment_palette)  # adjust if you need a different palette
        + labs(title = paste0(patient," - Single-cells Middle Slice FS Profiles"), x = "UMAP 0", y = "UMAP 1")
        + theme_bw()
        + umap_theme
        + guides(
            color = guide_legend(
                title = "Treatment",
                override.aes = list(alpha = 1, size = 5),
                ncol = 1
            )
        )
    )

    print(umap_sc_2D_middle_slice_color_by_treatment_plot)

    single_cell_MIP_features_path <- file.path(figures_path, "2D", "max_projection", paste0(patient, "_sc_fs_umap_color_by_treatment.png"))
    ggsave(umap_sc_2D_MIP_color_by_treatment_plot, file = single_cell_MIP_features_path, width = width, height = height, dpi = 400)
    single_cell_MS_features_path <- file.path(figures_path, "2D", "middle_slice", paste0(patient, "_sc_fs_umap_color_by_treatment.png"))
    ggsave(umap_sc_2D_middle_slice_color_by_treatment_plot, file = single_cell_MS_features_path, width = width, height = height, dpi = 400)
}
