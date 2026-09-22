list_of_packages <- c("ggplot2", "dplyr", "tidyr", "patchwork", "arrow", "RColorBrewer", "ComplexHeatmap", "circlize", "grid")
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

results_path <- file.path(root_dir, "4.linear_modeling/results/variate_importance")
figures_path <- file.path(root_dir, "4.linear_modeling/figures/variate_importance")
dir.create(figures_path, recursive = TRUE, showWarnings = FALSE)

venn_df <- arrow::read_parquet(file.path(results_path, "variate_hit_venn_regions_all_scopes.parquet"))
upset_df <- arrow::read_parquet(file.path(results_path, "variate_hit_upset_counts_all_scopes.parquet"))
sizes_df <- arrow::read_parquet(file.path(results_path, "variate_hit_set_sizes_all_scopes.parquet"))

# every scope: the columns that define one Venn + UpSet pair
scopes <- list(
    all_models = character(0),
    per_patient_treatment = c("patient", "treatment"),
    per_patient = c("patient"),
    per_treatment = c("treatment"),
    per_treatment_tumor_type = c("treatment", "tumor_type"),
    per_tumor_type = c("tumor_type")
)
model_sets <- c("original", "technical")
UPSET_TOP_N <- 25
# UpSet column widths (0.5 in per unit): set-size bars, term labels (fits "cell_per_organoid_count"),
# and the bar-row height (fits the two-line y-axis title)
UPSET_SIZE_COL <- 4.2
UPSET_LABEL_COL <- 6.4
UPSET_BAR_ROW <- 3.4
upset_dims <- function(n_combo, n_set) {
    c(
        width = 0.5 * (UPSET_SIZE_COL + UPSET_LABEL_COL + max(min(n_combo, UPSET_TOP_N), 4) + 1.5) + 0.3,
        height = 3.8 + 0.42 * n_set
    )
}

term_colors <- c(
    treatment = "#d95f02", cell_count = "#1b9e77", organoid_count = "#7570b3",
    cell_per_organoid_count = "#e7298a", manhattan_distance_from_center = "#66a61e",
    cell_x_position = "#e6ab02", cell_y_position = "#a6761d",
    cell_z_position = "#666666", cell_z_depth = "#1f78b4"
)

# four-set Venn layout (plot coordinates): ellipse x, y, width, height, angle
ellipses <- data.frame(
    set = 1:4,
    x = c(0.350, 0.450, 0.544, 0.644),
    y = c(0.400, 0.500, 0.500, 0.400),
    w = 0.72, h = 0.45,
    angle = c(140, 140, 40, 40)
)
# where the count of each region (bits in term order) is written
region_xy <- data.frame(
    region_key = c(
        "1000", "0100", "0010", "0001", "1100", "1010", "1001", "0110",
        "0101", "0011", "1110", "1101", "1011", "0111", "1111"
    ),
    x = c(0.14, 0.32, 0.68, 0.85, 0.23, 0.29, 0.50, 0.50, 0.71, 0.77, 0.35, 0.61, 0.39, 0.65, 0.50),
    y = c(0.42, 0.72, 0.72, 0.42, 0.59, 0.30, 0.17, 0.66, 0.30, 0.59, 0.50, 0.24, 0.24, 0.50, 0.38)
)
set_label_xy <- data.frame(x = c(0.02, 0.22, 0.78, 0.98), y = c(0.86, 0.97, 0.97, 0.86), hjust = c(0, 0.5, 0.5, 1))

ellipse_polygon <- function(x, y, w, h, angle, set, n = 200) {
    t <- seq(0, 2 * pi, length.out = n)
    theta <- angle * pi / 180
    data.frame(
        set = set,
        px = x + (w / 2) * cos(t) * cos(theta) - (h / 2) * sin(t) * sin(theta),
        py = y + (w / 2) * cos(t) * sin(theta) + (h / 2) * sin(t) * cos(theta)
    )
}

plot_venn <- function(regions, sizes, title) {
    # four-set Venn of feature membership
    terms <- sizes$term
    stopifnot(length(terms) == 4)
    polys <- do.call(rbind, lapply(seq_len(4), function(i) do.call(ellipse_polygon, as.list(ellipses[i, ]))))
    polys$term <- factor(terms[polys$set], levels = terms)
    counts <- region_xy |>
        left_join(regions[, c("region_key", "n_features")], by = "region_key") |>
        mutate(n_features = ifelse(is.na(n_features), 0, n_features))
    labels <- cbind(set_label_xy, term = terms, n = sizes$set_size)
    labels$label <- paste0(labels$term, "\n(", labels$n, " features)")
    (
        ggplot()
        + geom_polygon(data = polys, aes(px, py, group = term, fill = term), alpha = 0.3)
        + geom_path(data = polys, aes(px, py, group = term, colour = term), linewidth = 0.6)
        + geom_text(data = counts, aes(x, y, label = n_features), size = 6.5)
        + geom_text(
            data = labels, aes(x, y, label = label, colour = term, hjust = hjust),
            fontface = "bold", size = 5, lineheight = 0.9
        )
        + scale_fill_manual(values = term_colors[terms])
        + scale_colour_manual(values = term_colors[terms])
        + coord_fixed(xlim = c(0, 1), ylim = c(0, 1))
        + labs(title = paste0(
            "Features that are significant for each variate (", title, ")\n",
            sizes$n_none[1], " of ", sizes$n_total[1],
            " features are significant for none of these ", length(terms), " variates"
        ))
        + theme_void()
        + theme(legend.position = "none", plot.title = element_text(size = 14))
    )
}

plot_upset <- function(combos, sizes, title, top_n = UPSET_TOP_N, show_title = TRUE) {
    # UpSet plot: bars = features in exactly that term combination, dots = which terms.
    # show_title = FALSE drops the overall plot_annotation title (used for the headline-figure
    # panel, which embeds this plot without a title even though the standalone PDF keeps one).
    terms <- sizes$term
    n_set <- length(terms)
    shown <- head(combos, top_n)
    n_combo <- nrow(shown)
    shown$x <- seq_len(n_combo)
    note <- if (n_combo < nrow(combos)) paste0("top ", n_combo, " of ", nrow(combos), " combinations; ") else ""

    p_bar <- (
        ggplot(shown, aes(x = x))
        + geom_rect(
            aes(xmin = x - 0.45, xmax = x + 0.45, ymin = 0.8, ymax = n_features, fill = treatment_specific)
        )
        + geom_text(aes(y = n_features * 1.08, label = n_features), vjust = 0, size = 5)
        + scale_fill_manual(values = c("TRUE" = term_colors[["treatment"]], "FALSE" = "#555555"))
        + scale_y_log10(expand = c(0, 0))
        + coord_cartesian(xlim = c(0.4, n_combo + 0.6), ylim = c(0.8, max(shown$n_features) * 3), expand = FALSE)
        + labs(x = NULL, y = "features in exactly\nthis combination (log)")
        + theme_classic(base_size = 16)
        + theme(
            legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.line.x = element_blank(), plot.margin = margin(2, 2, 0, 2)
        )
    )

    membership <- shown |>
        select(x, all_of(paste0("in_", terms))) |>
        pivot_longer(-x, names_to = "term", values_to = "on") |>
        mutate(term = sub("^in_", "", term), row = match(term, terms))
    stripes <- data.frame(row = seq_len(n_set), shade = seq_len(n_set) %% 2 == 1)
    lines_df <- membership |> filter(on) |> group_by(x) |> filter(n() > 1) |> ungroup()
    p_mat <- (
        ggplot()
        + geom_rect(
            data = filter(stripes, shade),
            aes(xmin = 0.4, xmax = n_combo + 0.6, ymin = row - 0.5, ymax = row + 0.5), fill = "#f2f2f2"
        )
        + geom_point(data = membership, aes(x, row), colour = "#dddddd", size = 5.5)
        + geom_line(data = lines_df, aes(x, row, group = x), linewidth = 0.9)
        + geom_point(data = filter(membership, on), aes(x, row, fill = term), colour = "black", shape = 21, size = 5.5)
        + scale_fill_manual(values = term_colors[terms])
        + scale_y_reverse(limits = c(n_set + 0.5, 0.5), expand = c(0, 0))
        + scale_x_continuous(limits = c(0.4, n_combo + 0.6), expand = c(0, 0))
        + theme_void()
        + theme(legend.position = "none")
    )

    sizes$row <- seq_len(n_set)
    p_size <- (
        ggplot(sizes, aes(y = row))
        + geom_rect(aes(xmin = 0, xmax = set_size, ymin = row - 0.45, ymax = row + 0.45), fill = "#555555")
        + geom_text(aes(x = set_size, label = paste0(set_size, " ")), hjust = 1, size = 4.5)
        + scale_x_reverse(limits = c(max(sizes$set_size) * 1.35, 0))
        + scale_y_reverse(limits = c(n_set + 0.5, 0.5), expand = c(0, 0))
        + labs(x = "features in group", y = NULL)
        + theme_classic(base_size = 15)
        + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), plot.margin = margin(0, 0, 2, 2))
    )
    p_lab <- (
        ggplot(sizes, aes(y = row))
        + geom_text(aes(x = 0.02, label = term, colour = term), hjust = 0, fontface = "bold", size = 5)
        + scale_colour_manual(values = term_colors[terms])
        + scale_x_continuous(limits = c(0, 1), expand = c(0, 0))
        + scale_y_reverse(limits = c(n_set + 0.5, 0.5), expand = c(0, 0))
        + theme_void()
        + theme(legend.position = "none", plot.margin = margin(0, 2, 0, 2))
    )

    combined <- (
        (plot_spacer() + plot_spacer() + p_bar + p_size + p_lab + p_mat)
        + plot_layout(ncol = 3, widths = c(UPSET_SIZE_COL, UPSET_LABEL_COL, max(n_combo, 4) + 1.5), heights = c(UPSET_BAR_ROW, n_set * 0.55))
    )
    if (!show_title) return(combined)
    combined + plot_annotation(
        title = paste0("Feature membership across variate groups, ", title, " (", note, "orange = treatment only)"),
        theme = theme(plot.title = element_text(size = 14), plot.margin = margin(4, 4, 2, 4))
    )
}

group_title <- function(group_cols, row) {
    if (length(group_cols) == 0) return("all models")
    paste(sprintf("%s=%s", group_cols, unlist(row[group_cols])), collapse = ", ")
}

# one multi-page pdf per scope (Venn then UpSet page for every group and model set); a pdf is only rebuilt when missing
for (scope in names(scopes)) {
    pdf_path <- if (scope == "all_models") {
        file.path(figures_path, "variate_importance.pdf")
    } else {
        file.path(figures_path, scope, paste0(scope, ".pdf"))
    }
    if (file.exists(pdf_path)) {
        cat("skipping", scope, "(pdf exists)\n")
        next
    }
    dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
    group_cols <- scopes[[scope]]
    plots <- list()
    widths <- c()
    heights <- c()
    for (model in model_sets) {
        sizes_scope <- sizes_df |> filter(model_set == model, scope == !!scope)
        groups <- sizes_scope |> select(all_of(group_cols)) |> distinct() |> arrange(across(all_of(group_cols)))
        if (length(group_cols) == 0) groups <- data.frame(dummy = 1)
        for (i in seq_len(nrow(groups))) {
            in_group <- function(df) {
                df <- df |> filter(model_set == model, scope == !!scope)
                for (col in group_cols) df <- df[df[[col]] == groups[[col]][i], ]
                df
            }
            title <- group_title(group_cols, groups[i, , drop = FALSE])
            sizes_venn <- in_group(sizes_df) |> filter(plot == "venn")
            sizes_upset <- in_group(sizes_df) |> filter(plot == "upset")
            combos <- in_group(upset_df)
            plots[[length(plots) + 1]] <- plot_venn(in_group(venn_df), sizes_venn, title)
            widths <- c(widths, 7.5)
            heights <- c(heights, 6)
            if (nrow(combos) > 0) {
                plots[[length(plots) + 1]] <- plot_upset(combos, sizes_upset, title)
                dims <- upset_dims(nrow(combos), nrow(sizes_upset))
                widths <- c(widths, dims[["width"]])
                heights <- c(heights, dims[["height"]])
            }
        }
    }
    save_plots_pdf(plots, pdf_path, width = widths, height = heights)
    cat(scope, ":", length(plots), "pages ->", pdf_path, "\n")
}

lm_results_path <- file.path(root_dir, "4.linear_modeling/results/linear_modeling")
cooccurrence_figures_path <- file.path(figures_path, "treatment_only_cooccurrence")
dir.create(cooccurrence_figures_path, recursive = TRUE, showWarnings = FALSE)
cooccurrence_pdf <- file.path(cooccurrence_figures_path, "treatment_only_cooccurrence.pdf")

# profile -> (original model file, technical model file)
lm_files <- list(
    organoid = c(original = "organoid_norm", technical = "organoid_norm_technical_model"),
    sc = c(original = "sc_norm", technical = "sc_norm_technical_model"),
    organoid_agg = c(original = "organoid_agg", technical = "organoid_agg_technical_model"),
    sc_agg = c(original = "single_cell_agg", technical = "single_cell_agg_technical_model")
)
meki_drugs <- names(treatment_moa_map)[treatment_moa_map == "MEK1/2 inhibitor"]

# annotation colours: patient, tumor type, drug (treatment) and dose come from the manuscript theme
patient_palette <- setNames(tab20_palette_for_patients[seq_along(tumor_type_lookup)], names(tumor_type_lookup))
dose_label_palette <- c("1uM" = dose_palette[["1"]], "10uM" = dose_palette[["10"]], "10nM" = "#8da0cb")

read_significant <- function(file_name) {
    df <- arrow::read_parquet(file.path(lm_results_path, paste0(file_name, ".parquet")))
    # the technical models name the columns / treatment term differently
    df <- df |>
        rename(any_of(c(patient = "Metadata_Biology_PatientTumor", treatment = "Metadata_Experiment_Treatment"))) |>
        mutate(term = ifelse(term == "Metadata_Experiment_Treatment", "treatment", term))
    df |>
        filter(pvalue_fdr < 0.05) |>
        group_by(patient, treatment, feature) |>
        # treatment is significant and it is the only variate that is
        summarise(treatment_only = all(term == "treatment"), Feature_type = first(Feature_type), .groups = "drop") |>
        filter(treatment_only)
}

cooccurrence_heatmap <- function(sig, title, show_title = TRUE, legend_title_size = 14, legend_label_size = 12) {
    # show_title = FALSE drops column_title (used for the headline-figure panel, which embeds
    # this heatmap without a title even though the standalone PDF keeps one).
    # legend_*_size set the legend font sizes; the headline panel enlarges them because it is shrunk to fit.
    legend_title_gp <- gpar(fontsize = legend_title_size, fontface = "bold")
    legend_labels_gp <- gpar(fontsize = legend_label_size)
    sig <- sig |> mutate(patient_treatment = paste(patient, treatment, sep = " | "))
    columns <- sig |>
        distinct(patient, treatment, patient_treatment) |>
        mutate(
            tumor_type = tumor_type_lookup[patient],
            drug = sub("_.*$", "", treatment),
            dose = sub("^.*_", "", treatment)
        )
    features <- sig |> distinct(feature, Feature_type) |> mutate(Feature_type = ifelse(is.na(Feature_type), "Other", Feature_type))
    mat <- matrix(0, nrow = nrow(features), ncol = nrow(columns), dimnames = list(features$feature, columns$patient_treatment))
    mat[cbind(match(sig$feature, features$feature), match(sig$patient_treatment, columns$patient_treatment))] <- 1

    feature_type_colors <- c(feature_type_palette, Other = "#999999")
    top_annotation <- HeatmapAnnotation(
        Patient = columns$patient, `Tumor type` = columns$tumor_type,
        Treatment = columns$drug, Dose = columns$dose,
        col = list(Patient = patient_palette, `Tumor type` = tumor_type_palette,
                   Treatment = custom_treatment_palette, Dose = dose_label_palette),
        annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 14),
        annotation_legend_param = list(
            Treatment = list(ncol = 2, title_gp = legend_title_gp, labels_gp = legend_labels_gp),
            Patient = list(title_gp = legend_title_gp, labels_gp = legend_labels_gp),
            `Tumor type` = list(title_gp = legend_title_gp, labels_gp = legend_labels_gp),
            Dose = list(title_gp = legend_title_gp, labels_gp = legend_labels_gp)
        ),
        simple_anno_size = unit(5, "mm")
    )
    left_annotation <- rowAnnotation(
        `Feature type` = features$Feature_type,
        col = list(`Feature type` = feature_type_colors),
        annotation_name_gp = gpar(fontsize = 14), simple_anno_size = unit(5, "mm"),
        annotation_legend_param = list(`Feature type` = list(title_gp = legend_title_gp, labels_gp = legend_labels_gp))
    )
    # the matrix is a yes/no occurrence, so it gets a discrete two-colour legend (clustering uses the 0/1 matrix)
    row_hc <- if (nrow(mat) > 1) hclust(dist(mat, method = "binary"), method = "average") else FALSE
    col_hc <- if (ncol(mat) > 1) hclust(dist(t(mat), method = "binary"), method = "average") else FALSE
    sig_mat <- matrix(ifelse(mat == 1, "yes", "no"), nrow = nrow(mat), dimnames = dimnames(mat))
    Heatmap(
        sig_mat, name = "treatment-only significant",
        col = c(no = "white", yes = "#d95f02"),
        heatmap_legend_param = list(
            title = "treatment\nonly terms", at = c("yes", "no"), border = "black",
            title_gp = legend_title_gp, labels_gp = legend_labels_gp
        ),
        cluster_rows = row_hc, cluster_columns = col_hc,
        top_annotation = top_annotation, left_annotation = left_annotation,
        show_row_names = nrow(mat) <= 100, row_names_gp = gpar(fontsize = 8),
        show_column_names = FALSE, use_raster = TRUE, raster_quality = 2,
        column_title = if (show_title) paste0(title, "\n", nrow(mat), " features x ", ncol(mat), " patient-treatments") else NULL,
        column_title_gp = gpar(fontsize = 16, fontface = "bold"),
        row_title = "feature (clustered)", row_title_gp = gpar(fontsize = 14), column_title_side = "top"
    )
}

occurrence_tables <- list()
pages <- list()
for (profile in names(lm_files)) {
    for (model in model_sets) {
        sig <- read_significant(lm_files[[profile]][[model]])
        occurrence_tables[[length(occurrence_tables) + 1]] <- sig |>
            select(patient, treatment, feature, Feature_type, treatment_only) |>
            mutate(profile = profile, model_set = model)
        label <- paste0("Treatment-only significant: ", profile, ", ", model, " model")
        meki_sig <- sig |> filter(sub("_.*$", "", treatment) %in% meki_drugs)
        # the full matrix, then the same matrix subset to the MEKi treatments
        pages[[length(pages) + 1]] <- heatmap_grid_page(list(cooccurrence_heatmap(sig, label)))
        pages[[length(pages) + 1]] <- heatmap_grid_page(list(cooccurrence_heatmap(meki_sig, paste0(label, " (MEKi treatments only)"))))
        cat(profile, model, ":", nrow(sig), "treatment-only significant,", nrow(meki_sig), "in MEKi treatments\n")
    }
}
save_plots_pdf(pages, cooccurrence_pdf, width = 15, height = 10)
cat(length(pages), "pages ->", cooccurrence_pdf, "\n")
arrow::write_parquet(
    bind_rows(occurrence_tables),
    file.path(results_path, "treatment_only_cooccurrence.parquet")
)

headline_panels_path <- file.path(figures_path, "headline_panels")
dir.create(headline_panels_path, recursive = TRUE, showWarnings = FALSE)

# panel: tumor-type UpSet, tumor_type = pNF, original model
pnf_sizes <- sizes_df |> filter(model_set == "original", scope == "per_tumor_type", tumor_type == "pNF")
pnf_upset_sizes <- pnf_sizes |> filter(plot == "upset")
pnf_combos <- upset_df |> filter(model_set == "original", scope == "per_tumor_type", tumor_type == "pNF")
p_pnf_upset <- plot_upset(pnf_combos, pnf_upset_sizes, show_title = FALSE)
pnf_dims <- upset_dims(nrow(pnf_combos), nrow(pnf_upset_sizes))
ggsave(
    file.path(headline_panels_path, "tumor_type_pNF_upset.png"), p_pnf_upset,
    width = pnf_dims[["width"]], height = pnf_dims[["height"]], dpi = 600, bg = "white"
)

# panel: treatment-only co-occurrence heatmap, organoid, original model
sig_organoid <- read_significant(lm_files[["organoid"]][["original"]])
ht_organoid_notitle <- cooccurrence_heatmap(
    sig_organoid, title = NULL, show_title = FALSE, legend_title_size = 22, legend_label_size = 20
)
png(file.path(headline_panels_path, "cooccurrence_organoid_original.png"), width = 15, height = 10, units = "in", res = 600)
ComplexHeatmap::draw(ht_organoid_notitle)
dev.off()

cat("headline panels ->", headline_panels_path, "\n")
