#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
fig_repro_dir <- file.path(repo_dir, "figure_reproduction")
data_dir <- file.path(fig_repro_dir, "data")

dispersion_csv <- file.path(data_dir, "res_var_excl_cluster24.csv")

out_dir <- file.path(fig_repro_dir, "outputs")
ranked_dir <- file.path(out_dir, "ranked")
retained_matrix_csv <- file.path(out_dir, "retained_gene_matrix.csv")

ribo_pattern <- "^RP"
n_quantiles <- 5   # source Rmd line 96

# ANALYSIS_OK[experimental-design]: the 3 pooled iPSC lines are the fixed experiment (data/DATA_MANIFEST.md chdc-scrnaseq), not an exclusion list; asserted against matrix columns downstream
individuals <- c("NA19099", "NA19093", "NA18522")
cell_types <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)
n_expected_cols <- length(individuals) * length(cell_types)  # 18

# Membership rule for the consistently-low/high sets: quintile 1 (or 5) in all
# 6 cell types (source uses rowSums(. == q) >= 6).
n_cell_types_required <- length(cell_types)

# ------------------------------- HELPERS -------------------------------------

source(file.path(fig_repro_dir, "scripts", "utils.R"))

# Rank the given genes (lowest dispersion -> lowest rank) within one cell type,
# using the per-CT median dispersion across the 3 individuals (source lines
# 83-88). Returns a named numeric rank vector over `genes_to_rank`.
# ANALYSIS_OK[local-helper]: used within this standalone pipeline-stage script (run via run.sh); not a shared-library export
rank_cell_type <- function(res_var_mat, cell_type, genes_to_rank) {
  cols <- ct_columns(cell_type)
  ct_res_var <- res_var_mat[genes_to_rank, cols, drop = FALSE]
  median_res_var <- apply(ct_res_var, 1, median)
  ranks <- rank(median_res_var, ties.method = "average")  # source line 88
  ranks
}

# Cut a rank vector into n_quantiles equal-count bins (source lines 98-103).
# ANALYSIS_OK[local-helper]: used within this standalone pipeline-stage script (run via run.sh); not a shared-library export
ranks_to_quintiles <- function(ranks) {
  cut(ranks,
      breaks = quantile(ranks, probs = seq(0, 1, length.out = n_quantiles + 1),
                        na.rm = TRUE),
      labels = FALSE,
      include.lowest = TRUE)
}

# consistently-low (target_quintile = 1) or consistently-high (= 5) gene set:
# genes assigned target_quintile in ALL cell types (NA quintiles never count).
# ANALYSIS_OK[local-helper]: used within this standalone pipeline-stage script (run via run.sh); not a shared-library export
consistent_genes <- function(quintile_df, target_quintile) {
  hits <- rowSums(quintile_df == target_quintile, na.rm = TRUE)
  rownames(quintile_df)[hits >= n_cell_types_required]
}

# ANALYSIS_OK[local-helper]: used within this standalone pipeline-stage script (run via run.sh); not a shared-library export
jaccard <- function(a, b) {
  if (length(union(a, b)) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

# --------------------------------- MAIN --------------------------------------

dir.create(ranked_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load dispersion matrix and drop ^RP genes (source lines 50-55) ----------
if (!file.exists(dispersion_csv)) stop("Missing dispersion matrix: ", dispersion_csv)
res_var_df <- read.csv(dispersion_csv, row.names = 1, check.names = FALSE)
message(sprintf("[load] dispersion %d genes x %d cols", nrow(res_var_df), ncol(res_var_df)))
if (ncol(res_var_df) != n_expected_cols) {
  stop("Dispersion matrix has ", ncol(res_var_df), " cols; expected ", n_expected_cols)
}
if (any(is.na(res_var_df))) stop("Dispersion matrix contains NA values")
res_var_mat <- as.matrix(res_var_df)

ribo_genes <- rownames(res_var_mat)[grepl(ribo_pattern, rownames(res_var_mat))]
# ANALYSIS_OK[filtering]: drop ^RP ribosomal genes (repo convention, source lines 50-55); dropped + retained counts logged immediately below
res_var_mat <- res_var_mat[!rownames(res_var_mat) %in% ribo_genes, , drop = FALSE]
message(sprintf("[ribo] dropped %d genes; %d non-ribo genes remain",
                length(ribo_genes), nrow(res_var_mat)))
all_genes <- rownames(res_var_mat)

# --- Load per-CT retained indicator from script 01 ---------------------------
if (!file.exists(retained_matrix_csv)) {
  stop("Missing retained_gene_matrix.csv — run 01_exclude_topexpressed_genes.R first")
}
retained_ind <- read.csv(retained_matrix_csv, row.names = 1, check.names = FALSE)
retained_ind <- as.matrix(retained_ind)
if (!setequal(colnames(retained_ind), cell_types)) {
  stop("retained_gene_matrix.csv columns do not match the cell-type set")
}
# ANALYSIS_OK[column-reorder]: reorder columns to canonical cell-type order only; the column set was asserted equal to cell_types at line above (setequal); no genes/rows dropped
retained_ind <- retained_ind[, cell_types, drop = FALSE]
if (!all(rownames(retained_ind) %in% all_genes)) {
  stop("retained_gene_matrix.csv contains genes absent from the non-ribo dispersion matrix")
}

# --- Build ORIGINAL quintile table (no exclusion; all non-ribo genes) --------
# Faithful reproduction of source lines 72-106.
orig_quintiles <- matrix(NA_integer_, nrow = length(all_genes),
                         ncol = length(cell_types),
                         dimnames = list(all_genes, cell_types))
for (cell_type in cell_types) {
  ranks <- rank_cell_type(res_var_mat, cell_type, all_genes)
  orig_quintiles[, cell_type] <- ranks_to_quintiles(ranks)
}
orig_quintiles <- as.data.frame(orig_quintiles)
message("[orig] built full-universe quintile table")

# --- Build RETAINED quintile table (per-CT top-expressed genes excluded) -----
# For each cell type, ranks/quintiles are computed among that cell type's
# retained genes only; genes excluded in a cell type are left NA in its column.
retained_quintiles <- matrix(NA_integer_, nrow = length(all_genes),
                             ncol = length(cell_types),
                             dimnames = list(all_genes, cell_types))
for (cell_type in cell_types) {
  ct_retained <- rownames(retained_ind)[retained_ind[, cell_type] == 1]
  ct_retained <- intersect(all_genes, ct_retained)
  message(sprintf("[rank] %-25s ranking %d retained genes",
                  cell_type, length(ct_retained)))
  ranks <- rank_cell_type(res_var_mat, cell_type, ct_retained)
  retained_quintiles[ct_retained, cell_type] <- ranks_to_quintiles(ranks)
}
retained_quintiles <- as.data.frame(retained_quintiles)

write.csv(orig_quintiles,     file.path(ranked_dir, "var_quintile_dfs_original.csv"))
write.csv(retained_quintiles, file.path(ranked_dir, "var_quintile_dfs_retained.csv"))

# --- Re-derive consistently-low / consistently-high gene sets ----------------
orig_low  <- consistent_genes(orig_quintiles,     target_quintile = 1)
orig_high <- consistent_genes(orig_quintiles,     target_quintile = n_quantiles)
new_low   <- consistent_genes(retained_quintiles, target_quintile = 1)
new_high  <- consistent_genes(retained_quintiles, target_quintile = n_quantiles)

writeLines(new_low,  file.path(ranked_dir, "consistently_low_retained.txt"))
writeLines(new_high, file.path(ranked_dir, "consistently_high_retained.txt"))
writeLines(orig_low,  file.path(ranked_dir, "consistently_low_original.txt"))
writeLines(orig_high, file.path(ranked_dir, "consistently_high_original.txt"))

# --- Compare original vs retained gene sets ----------------------------------
comparison <- data.frame(
  gene_set = c("consistently_low", "consistently_high"),
  n_original = c(length(orig_low), length(orig_high)),
  n_retained = c(length(new_low), length(new_high)),
  n_shared = c(length(intersect(orig_low, new_low)),
               length(intersect(orig_high, new_high))),
  n_lost = c(length(setdiff(orig_low, new_low)),
             length(setdiff(orig_high, new_high))),
  n_gained = c(length(setdiff(new_low, orig_low)),
               length(setdiff(new_high, orig_high))),
  jaccard = c(jaccard(orig_low, new_low), jaccard(orig_high, new_high)),
  stringsAsFactors = FALSE)

write.csv(comparison, file.path(ranked_dir, "consistent_set_comparison.csv"),
          row.names = FALSE)

message("[done] wrote quintile tables, gene sets, and comparison to ", ranked_dir)
print(comparison)
