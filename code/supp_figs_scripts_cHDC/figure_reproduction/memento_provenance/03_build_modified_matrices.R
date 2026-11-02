#!/usr/bin/env Rscript

repo_dir     <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
cluster_dir  <- file.path(repo_dir, "cellranger_cluster-mode_trial")
memento_dir  <- file.path(cluster_dir, "analysis/datasets/memento")
analysis_dir <- file.path(repo_dir, "analysis/epithelial_subcluster24_excl_sensitivity")
out_dir      <- file.path(analysis_dir, "outputs")

published_res_var_csv <- file.path(
  memento_dir, "new_res_var_by_cell_type_ind_1000_cts_new_predicted_id.csv")
published_mean_csv <- file.path(
  memento_dir, "new_mean_expression_by_cell_type_ind_1000_cts_new_predicted_id.csv")

new_epi_res_var_csv <- file.path(out_dir, "res_var_epithelial_excl_cluster24.csv")
new_epi_mean_csv    <- file.path(out_dir, "mean_expression_epithelial_excl_cluster24.csv")

# ANALYSIS_OK[fixed-design]: the 3 pooled cHDC iPSC lines (the entire study
# cohort), not an exclusion filter -- matches the identical constant in
# subcluster_dispersion_robustness/scripts/01_subcluster_and_export.R.
individuals <- c("NA19099", "NA19093", "NA18522")
cell_types <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)
n_expected_cols <- length(individuals) * length(cell_types)  # 18
epi_cols <- paste0(individuals, "_Epithelial_cells")
min_shared_genes <- 100  # sanity floor on the gene-row intersection size

source(file.path(analysis_dir, "scripts", "utils.R"))
log_msg <- make_logger("03_build_matrices")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load ---------------------------------------------------------------
for (f in c(published_res_var_csv, published_mean_csv, new_epi_res_var_csv, new_epi_mean_csv)) {
  if (!file.exists(f)) stop("Missing input file: ", f)
}

# ANALYSIS_OK[within-file-helper]: load_published/load_new_epi/build_pair are
# called later in this same file (res_var_pair/mean_pair below); R026's
# cross-file rule cannot see intra-file callers -- matches the identical
# waiver pattern in subcluster_dispersion_robustness/scripts/01_subcluster_and_export.R.
load_published <- function(path) {
  df <- read.csv(path, row.names = 1, check.names = FALSE)
  if (ncol(df) != n_expected_cols) {
    stop(sprintf("%s has %d columns; expected %d (%d individuals x %d cell types).",
                 basename(path), ncol(df), n_expected_cols,
                 length(individuals), length(cell_types)))
  }
  as.matrix(df)
}

# ANALYSIS_OK[within-file-helper]: called later in this same file (build_pair
# below); R026's cross-file rule cannot see intra-file callers.
load_new_epi <- function(path) {
  df <- read.csv(path, row.names = 1, check.names = FALSE)
  if (!setequal(colnames(df), epi_cols)) {
    stop(sprintf("%s columns are %s; expected exactly %s.",
                 basename(path), paste(colnames(df), collapse = ", "),
                 paste(epi_cols, collapse = ", ")))
  }
  as.matrix(df[, epi_cols, drop = FALSE])
}

# ANALYSIS_OK[within-file-helper]: called later in this same file
# (res_var_pair/mean_pair below); R026's cross-file rule cannot see intra-file
# callers.
build_pair <- function(published_path, new_epi_path, tag) {
  published_mat <- load_published(published_path)
  new_epi_mat <- load_new_epi(new_epi_path)
  log_msg("[%s] published matrix: %d genes x %d cols.", tag,
          nrow(published_mat), ncol(published_mat))
  log_msg("[%s] fresh Epithelial re-estimate: %d genes x %d cols.", tag,
          nrow(new_epi_mat), ncol(new_epi_mat))

  shared_genes <- intersect(rownames(published_mat), rownames(new_epi_mat))
  n_pub_only <- length(setdiff(rownames(published_mat), rownames(new_epi_mat)))
  n_new_only <- length(setdiff(rownames(new_epi_mat), rownames(published_mat)))
  log_msg(paste0(
    "[%s] gene-row intersection: %d shared genes kept; %d published-only genes ",
    "dropped (no usable fresh Epithelial estimate -- likely failed Memento's ",
    "min_perc_group filter on the cluster-24-excluded cell set); %d ",
    "recompute-only genes ignored (not in the published gene universe)."),
    tag, length(shared_genes), n_pub_only, n_new_only)
  if (length(shared_genes) < min_shared_genes) {
    stop(sprintf("[%s] Only %d shared genes survive -- too few for a dream refit. ",
                 tag, length(shared_genes)))
  }

  baseline_mat <- published_mat[shared_genes, , drop = FALSE]
  excl_mat <- baseline_mat
  excl_mat[, epi_cols] <- new_epi_mat[shared_genes, epi_cols]

  # Sanity: the two matrices must differ ONLY in the 3 Epithelial columns.
  non_epi_cols <- setdiff(colnames(baseline_mat), epi_cols)
  if (!isTRUE(all.equal(baseline_mat[, non_epi_cols], excl_mat[, non_epi_cols]))) {
    stop(sprintf("[%s] Non-Epithelial columns differ between baseline and ",
                 "excl_cluster24 -- splice logic bug.", tag))
  }
  if (isTRUE(all.equal(baseline_mat[, epi_cols], excl_mat[, epi_cols]))) {
    stop(sprintf("[%s] Epithelial columns are IDENTICAL between baseline and ",
                 "excl_cluster24 -- the substitution did not take effect.", tag))
  }
  log_msg("[%s] Verified: non-Epithelial columns identical, Epithelial columns differ.", tag)

  list(baseline = baseline_mat, excl_cluster24 = excl_mat)
}

res_var_pair <- build_pair(published_res_var_csv, new_epi_res_var_csv, "DD/res_var")
mean_pair    <- build_pair(published_mean_csv, new_epi_mean_csv, "DE/mean_expression")

write.csv(res_var_pair$baseline, file.path(out_dir, "res_var_baseline.csv"))
write.csv(res_var_pair$excl_cluster24, file.path(out_dir, "res_var_excl_cluster24.csv"))
write.csv(mean_pair$baseline, file.path(out_dir, "mean_expression_baseline.csv"))
write.csv(mean_pair$excl_cluster24, file.path(out_dir, "mean_expression_excl_cluster24.csv"))

log_msg("Wrote res_var_{baseline,excl_cluster24}.csv (%d genes) to %s",
        nrow(res_var_pair$baseline), out_dir)
log_msg("Wrote mean_expression_{baseline,excl_cluster24}.csv (%d genes) to %s",
        nrow(mean_pair$baseline), out_dir)
log_msg("DONE. Next: scripts/04_dream_dd_de_refit_and_compare.R")
