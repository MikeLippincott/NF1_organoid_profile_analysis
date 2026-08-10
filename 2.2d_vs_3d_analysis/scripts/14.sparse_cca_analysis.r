list_of_packages <- c("arrow", "dplyr", "tidyr", "PMA")

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

path_2d_organoid <- file.path(root_dir, "data/2D_profiles/all_patient_profiles/max_projection/organoid_agg_profiles.parquet")
path_3d_organoid <- file.path(root_dir, "data/3D_profiles/all_patient_profiles/organoid_agg_profiles.parquet")

results_dir <- file.path(root_dir, "2.2d_vs_3d_analysis/results/sparse_cca")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

df_2d <- arrow::read_parquet(path_2d_organoid)
df_3d <- arrow::read_parquet(path_3d_organoid)

merge_cols <- c("Metadata_patient", "Metadata_tumor", "Metadata_Well")

keys_2d <- df_2d %>% select(all_of(merge_cols)) %>% distinct()
keys_3d <- df_3d %>% select(all_of(merge_cols)) %>% distinct()
shared_keys <- inner_join(keys_2d, keys_3d, by = merge_cols)

if (nrow(shared_keys) == 0) {
    stop("No matched patient-tumor-well combinations between 2D and 3D.")
}

df_2d_matched <- df_2d %>%
    inner_join(shared_keys, by = merge_cols) %>%
    arrange(across(all_of(merge_cols)))
df_3d_matched <- df_3d %>%
    inner_join(shared_keys, by = merge_cols) %>%
    arrange(across(all_of(merge_cols)))

stopifnot(identical(
    df_2d_matched %>% select(all_of(merge_cols)),
    df_3d_matched %>% select(all_of(merge_cols))
))

metadata <- df_2d_matched %>% select(all_of(merge_cols))
cat("Matched samples (pooled across all patients):", nrow(metadata), "\n")

# Cleans a profile dataframe into a numeric feature matrix ready for CCA:
# drops metadata/constant/all-NA columns, then mean-imputes remaining NAs.
prep_features <- function(df) {
    meta_cols <- grep("^Metadata_", colnames(df), value = TRUE)
    feat_df <- df %>% select(-all_of(meta_cols))

    numeric_cols <- sapply(feat_df, is.numeric)
    feat_df <- feat_df[, numeric_cols, drop = FALSE]

    # drop constant columns (zero variance, uninformative for CCA)
    constant_cols <- sapply(feat_df, function(x) length(unique(na.omit(x))) <= 1)
    feat_df <- feat_df[, !constant_cols, drop = FALSE]

    # drop all-NA columns, then mean-impute remaining NAs
    # (PMA::CCA requires complete matrices, no NA support)
    all_na_cols <- sapply(feat_df, function(x) all(is.na(x)))
    feat_df <- feat_df[, !all_na_cols, drop = FALSE]
    feat_df <- as.data.frame(lapply(feat_df, function(x) {
        if (any(is.na(x))) x[is.na(x)] <- mean(x, na.rm = TRUE)
        x
    }))

    feat_df
}

set.seed(0)
# 1 component: 3 was only ever exploratory (no rationale behind that number),
# and CCA.permute only validates the first canonical pair with a p-value, so
# there was never a real reason to keep the unvalidated CC2/CC3 around.
n_components <- 1

perm_out <- CCA.permute(
    x = mat_2d,
    z = mat_3d,
    typex = "standard",
    typez = "standard",
    standardize = TRUE,
    nperms = 25 # CCA.permute's own package default
)

cat("Best penalty (2D view):", perm_out$bestpenaltyx, "\n")
cat("Best penalty (3D view):", perm_out$bestpenaltyz, "\n")

permute_summary <- data.frame(
    penaltyx = perm_out$penaltyxs,
    penaltyz = perm_out$penaltyzs,
    cor = perm_out$cors,
    zstat = perm_out$zstats,
    pval = perm_out$pvals,
    n_nonzero_2d = perm_out$nnonzerous,
    n_nonzero_3d = perm_out$nnonzerovs
)
permute_summary$is_best <- (permute_summary$penaltyx == perm_out$bestpenaltyx) &
    (permute_summary$penaltyz == perm_out$bestpenaltyz)

arrow::write_parquet(permute_summary, file.path(results_dir, "cca_permute_summary.parquet"))

cca_out <- CCA(
    x = mat_2d,
    z = mat_3d,
    typex = "standard",
    typez = "standard",
    penaltyx = perm_out$bestpenaltyx,
    penaltyz = perm_out$bestpenaltyz,
    K = n_components,
    standardize = TRUE
)

scores_2d <- mat_2d %*% cca_out$u
scores_3d <- mat_3d %*% cca_out$v
colnames(scores_2d) <- paste0("CC", seq_len(n_components), "_2D")
colnames(scores_3d) <- paste0("CC", seq_len(n_components), "_3D")

canonical_scores <- cbind(metadata, as.data.frame(scores_2d), as.data.frame(scores_3d))
arrow::write_parquet(canonical_scores, file.path(results_dir, "canonical_scores.parquet"))

canonical_cors <- sapply(seq_len(n_components), function(k) {
    cor(scores_2d[, k], scores_3d[, k])
})

best_idx <- which(permute_summary$is_best)[1]
cor_summary <- data.frame(
    component = paste0("CC", seq_len(n_components)),
    canonical_correlation = canonical_cors,
    perm_pvalue = c(permute_summary$pval[best_idx], rep(NA, n_components - 1))
)

arrow::write_parquet(cor_summary, file.path(results_dir, "canonical_correlations.parquet"))
print(cor_summary)

# Save feature loadings (weights) per component, in long format.
loadings_2d <- as.data.frame(cca_out$u)
colnames(loadings_2d) <- paste0("CC", seq_len(n_components))
loadings_2d$feature <- colnames(mat_2d)
loadings_2d_long <- loadings_2d %>%
    tidyr::pivot_longer(cols = starts_with("CC"), names_to = "component", values_to = "weight")
arrow::write_parquet(loadings_2d_long, file.path(results_dir, "loadings_2d.parquet"))

loadings_3d <- as.data.frame(cca_out$v)
colnames(loadings_3d) <- paste0("CC", seq_len(n_components))
loadings_3d$feature <- colnames(mat_3d)
loadings_3d_long <- loadings_3d %>%
    tidyr::pivot_longer(cols = starts_with("CC"), names_to = "component", values_to = "weight")
arrow::write_parquet(loadings_3d_long, file.path(results_dir, "loadings_3d.parquet"))

cat("Sparse CCA results saved to:", results_dir, "\n")
