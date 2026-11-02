#!/usr/bin/env Rscript

suppressPackageStartupMessages({})

# ------------------------------- CONSTANTS -----------------------------------
# PATCHED for figure_reproduction (see PORT_CHANGELOG_qthreshold.md).
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
analysis_dir <- file.path(repo_dir, "figure_reproduction", "outputs", "q_threshold_sensitivity")
script_dir <- file.path(repo_dir, "figure_reproduction", "scripts", "q_threshold_sensitivity")
memento_dir <- file.path(analysis_dir, "memento")
quintile_dir <- file.path(analysis_dir, "quintiles")

q_values <- c("0.05", "0.1", "0.2", "0.25")  # 0.25 is the manuscript's own canonical value
ribo_pattern <- "^RP"
n_quantiles <- 5  # main_figs_1-4_related_supp.Rmd line 96

individuals <- c("NA19099", "NA19093", "NA18522")
cell_types <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)
n_expected_cols <- length(individuals) * length(cell_types)  # 18

source(file.path(script_dir, "utils.R"))

# ------------------------------- HELPERS -------------------------------------

# Rank genes within one cell type by median dispersion across the 3
# individuals (source lines 83-88).
# ANALYSIS_OK[local-helper]: used within this standalone pipeline-stage script (run via run.sh); not a shared-library export
rank_cell_type <- function(res_var_mat, cell_type, genes_to_rank) {
  cols <- ct_columns(cell_type)
  ct_res_var <- res_var_mat[genes_to_rank, cols, drop = FALSE]
  median_res_var <- apply(ct_res_var, 1, median)
  rank(median_res_var, ties.method = "average")
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

# --------------------------------- MAIN --------------------------------------

dir.create(quintile_dir, recursive = TRUE, showWarnings = FALSE)

for (q_str in q_values) {
  res_var_csv <- file.path(memento_dir, paste0("q", q_str), "res_var.csv")
  if (!file.exists(res_var_csv)) {
    stop("Missing res_var matrix for q=", q_str, ": ", res_var_csv,
         " -- run 01_memento_q_sweep.py for this Q_VALUE first")
  }
  res_var_df <- read.csv(res_var_csv, row.names = 1, check.names = FALSE)
  if (ncol(res_var_df) != n_expected_cols) {
    stop("res_var matrix for q=", q_str, " has ", ncol(res_var_df),
         " cols; expected ", n_expected_cols)
  }
  res_var_mat <- as.matrix(res_var_df)
  # ANALYSIS_OK[filtering]: drop ^RP ribosomal genes (repo convention, matches
  # dd_excl_topexpressed script 02); retained non-ribo count logged below
  res_var_mat <- res_var_mat[!grepl(ribo_pattern, rownames(res_var_mat)), , drop = FALSE]
  all_genes <- rownames(res_var_mat)
  message(sprintf("[q=%s] %d non-ribo genes", q_str, length(all_genes)))

  quintile_mat <- matrix(NA_integer_, nrow = length(all_genes), ncol = length(cell_types),
                         dimnames = list(all_genes, cell_types))
  for (cell_type in cell_types) {
    ranks <- rank_cell_type(res_var_mat, cell_type, all_genes)
    quintile_mat[, cell_type] <- ranks_to_quintiles(ranks)
  }
  quintile_df <- as.data.frame(quintile_mat)

  out_dir <- file.path(quintile_dir, paste0("q", q_str))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(quintile_df, file.path(out_dir, "quintile_df.csv"))
  message(sprintf("[q=%s] wrote quintile table (%d genes x %d cell types) to %s",
                  q_str, nrow(quintile_df), ncol(quintile_df), out_dir))
}

message("[done] quintile tables written for all q values: ", paste(q_values, collapse = ", "))
