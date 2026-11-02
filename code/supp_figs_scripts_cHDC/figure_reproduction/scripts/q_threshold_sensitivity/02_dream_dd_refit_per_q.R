#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)              # normalizeQuantiles
  library(variancePartition)  # dream, makeContrastsDream, eBayes, topTable
  library(BiocParallel)       # SnowParam
})

# ------------------------------- CONSTANTS -----------------------------------
# PATCHED for figure_reproduction (see PORT_CHANGELOG_qthreshold.md): outputs/
# script source repointed into this folder.
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
analysis_dir <- file.path(repo_dir, "figure_reproduction", "outputs", "q_threshold_sensitivity")
script_dir <- file.path(repo_dir, "figure_reproduction", "scripts", "q_threshold_sensitivity")
memento_dir <- file.path(analysis_dir, "memento")
refit_dir <- file.path(analysis_dir, "dream_refit")

ribo_pattern <- "^RP"
dd_adjp_thresh <- 0.05
n_dream_cores <- 4

individuals <- c("NA19099", "NA19093", "NA18522")
cell_types <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)
n_expected_cols <- length(individuals) * length(cell_types)  # 18

source(file.path(script_dir, "utils.R"))

if (Sys.getenv("Q_VALUE") == "") stop("Q_VALUE environment variable not set")
q_str <- Sys.getenv("Q_VALUE")
message(sprintf("[q=%s] starting dream DD refit", q_str))

# --------------------------------- MAIN --------------------------------------

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
# dd_excl_topexpressed script 03); retained non-ribo count logged immediately below
res_var_mat <- res_var_mat[!grepl(ribo_pattern, rownames(res_var_mat)), , drop = FALSE]
all_genes <- rownames(res_var_mat)
message(sprintf("[q=%s] %d non-ribo genes in this q's universe", q_str, length(all_genes)))

sample_cols <- colnames(res_var_mat)
col_ind <- sub("_.*$", "", sample_cols)
col_ct  <- sub("^[^_]*_", "", sample_cols)
stopifnot(all(col_ind %in% individuals), all(col_ct %in% cell_types),
          length(sample_cols) == n_expected_cols)
metadata <- data.frame(
  cell_type = factor(col_ct, levels = cell_types),
  cell_line = factor(col_ind, levels = individuals),
  row.names = sample_cols)

# 15 pairwise contrasts named compare_{B}_{A} (matches dd_excl_topexpressed
# script 03 / the published reformat_and_dream_de_dd_analysis.R convention).
pairs <- combn(cell_types, 2, simplify = FALSE)
contrast_matrix <- setNames(
  vapply(pairs, function(x) paste0("cell_type", x[2], " - ", "cell_type", x[1]),
         character(1)),
  vapply(pairs, function(x) paste0("compare_", x[2], "_", x[1]), character(1)))

# ANALYSIS_OK[model-formula]: dream mixed-model design formula (cell_type
# fixed effect only, 0+ absorbs the intercept), matching the published DD
# model reused verbatim from dd_excl_topexpressed script 03.
dream_formula <- ~ 0 + cell_type

mat <- normalizeQuantiles(res_var_mat[all_genes, sample_cols, drop = FALSE])
contrast_dream <- makeContrastsDream(dream_formula, metadata,
                                     contrasts = contrast_matrix)
param <- SnowParam(n_dream_cores, type = "SOCK",
                   progressbar = FALSE, exportglobals = FALSE)
message(sprintf("[q=%s] dream fit: %d genes x %d samples, %d contrasts",
                q_str, length(all_genes), length(sample_cols),
                length(contrast_matrix)))
fit <- dream(mat, dream_formula, metadata, contrast_dream, BPPARAM = param)
fit <- eBayes(fit)

dir.create(file.path(refit_dir, paste0("q", q_str)), recursive = TRUE, showWarnings = FALSE)
out_dir <- file.path(refit_dir, paste0("q", q_str))
saveRDS(fit, file.path(out_dir, "dream_fit.rds"))
writeLines(all_genes, file.path(out_dir, "gene_universe.txt"))

# Per-contrast DD gene sets (adj.P.Val <= dd_adjp_thresh).
dd_sets <- setNames(lapply(names(contrast_matrix), function(cn) {
  tt <- topTable(fit, coef = cn, number = Inf, sort.by = "none")
  rownames(tt)[!is.na(tt$adj.P.Val) & tt$adj.P.Val <= dd_adjp_thresh]
}), names(contrast_matrix))
saveRDS(dd_sets, file.path(out_dir, "dd_gene_sets.rds"))

dd_counts <- data.frame(
  contrast = names(dd_sets),
  n_dd = vapply(dd_sets, length, integer(1)),
  n_universe = length(all_genes),
  stringsAsFactors = FALSE)
write.csv(dd_counts, file.path(out_dir, "dd_counts.csv"), row.names = FALSE)

message(sprintf("[q=%s] done -- wrote dream fit + %d DD gene sets to %s",
                q_str, length(dd_sets), out_dir))
print(dd_counts)
