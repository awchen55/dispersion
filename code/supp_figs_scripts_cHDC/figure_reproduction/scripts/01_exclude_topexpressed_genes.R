#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
fig_repro_dir <- file.path(repo_dir, "figure_reproduction")
data_dir <- file.path(fig_repro_dir, "data")

dispersion_csv <- file.path(data_dir, "res_var_excl_cluster24.csv")
mean_csv <- file.path(data_dir, "mean_expression_excl_cluster24.csv")

out_dir <- file.path(fig_repro_dir, "outputs")
gene_list_dir <- file.path(out_dir, "gene_lists")
filt_mat_dir  <- file.path(out_dir, "filtered_matrices")

# Alexander's exclusion threshold: exclude genes with Memento mean > -6.
mean_excl_threshold <- -6

# Ribosomal-gene prefix filter, matching the source Rmd (grepl("^RP", ...)).
ribo_pattern <- "^RP"

# ANALYSIS_OK[experimental-design]: the 3 pooled iPSC lines are the fixed experiment (data/DATA_MANIFEST.md chdc-scrnaseq), not an exclusion list; asserted against matrix columns downstream
individuals <- c("NA19099", "NA19093", "NA18522")
cell_types <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)
n_expected_cols <- length(individuals) * length(cell_types)  # 18

# ------------------------------- HELPERS -------------------------------------

# Load a Memento matrix (genes x 18) with genes as row names; fail loudly on any
# structural surprise so a malformed input never propagates silently.
# ANALYSIS_OK[local-helper]: used within this standalone pipeline-stage script (run via run.sh); not a shared-library export
load_memento_matrix <- function(path, label) {
  if (!file.exists(path)) stop("Missing input matrix: ", path)
  m <- read.csv(path, row.names = 1, check.names = FALSE)
  message(sprintf("[load] %-11s %d genes x %d cols  (%s)",
                  label, nrow(m), ncol(m), basename(path)))
  if (ncol(m) != n_expected_cols) {
    stop(label, " matrix has ", ncol(m), " columns; expected ", n_expected_cols)
  }
  expected_cols <- as.vector(outer(individuals, cell_types,
                                   function(i, ct) paste0(i, "_", ct)))
  missing_cols <- setdiff(expected_cols, colnames(m))
  if (length(missing_cols) > 0) {
    stop(label, " matrix missing expected columns: ",
         paste(missing_cols, collapse = ", "))
  }
  if (any(is.na(m))) stop(label, " matrix contains NA values; investigate before use")
  as.matrix(m)
}

# Column names ({IND}_{CT}) for one cell type across the three individuals.
source(file.path(fig_repro_dir, "scripts", "utils.R"))

# --------------------------------- MAIN --------------------------------------

dir.create(gene_list_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(filt_mat_dir,  recursive = TRUE, showWarnings = FALSE)

res_var_mat <- load_memento_matrix(dispersion_csv, "dispersion")
mean_mat    <- load_memento_matrix(mean_csv, "mean")

# Both matrices must describe the same genes in the same order; the exclusion
# indexes dispersion genes by their mean value, so a mismatch would be silent.
if (!identical(rownames(res_var_mat), rownames(mean_mat))) {
  stop("Dispersion and mean matrices have non-identical row names/order")
}
if (!identical(colnames(res_var_mat), colnames(mean_mat))) {
  stop("Dispersion and mean matrices have non-identical column names/order")
}

# --- Drop ribosomal (^RP) genes from BOTH matrices (source Rmd lines 52-60) ---
ribo_genes <- rownames(res_var_mat)[grepl(ribo_pattern, rownames(res_var_mat))]
message(sprintf("[ribo] dropping %d genes matching '%s'", length(ribo_genes),
                ribo_pattern))
keep_nonribo <- !rownames(res_var_mat) %in% ribo_genes
# ANALYSIS_OK[filtering]: drop ^RP ribosomal genes (repo convention, source Rmd lines 52-60); dropped count logged at line above, retained count logged below
res_var_mat <- res_var_mat[keep_nonribo, , drop = FALSE]
mean_mat    <- mean_mat[keep_nonribo, , drop = FALSE]
message(sprintf("[ribo] retained %d non-ribosomal genes", nrow(res_var_mat)))
stopifnot(identical(rownames(res_var_mat), rownames(mean_mat)))

all_genes <- rownames(mean_mat)

# --- Per-cell-type exclusion of top-expressed genes --------------------------
# Exclusion rule (Brendan 2026-07-20): exclude a gene from a cell type if ANY of
#   its 3 individual mean values in that cell type is > -6. This is the more
#   aggressive of the interpretations considered (vs. per-CT median > -6): a gene
#   is dropped from a cell type if it is highly expressed in even one individual,
#   because a single poorly-mean-corrected individual is enough to distort that
#   gene's dispersion estimate in that cell type. The per-CT median rule is still
#   computed and reported per cell type for sensitivity comparison, but it does
#   NOT drive the retained-gene matrix or the filtered matrices below.

summary_rows <- list()
retained_indicator <- matrix(
  FALSE, nrow = length(all_genes), ncol = length(cell_types),
  dimnames = list(all_genes, cell_types))

for (cell_type in cell_types) {
  cols <- ct_columns(cell_type)
  ct_mean <- mean_mat[, cols, drop = FALSE]

  # PRIMARY rule: exclude if ANY individual's mean in this CT exceeds threshold.
  excluded_ct <- names(which(apply(ct_mean, 1,
                                   function(x) any(x > mean_excl_threshold))))
  retained_ct <- setdiff(all_genes, excluded_ct)

  # SENSITIVITY rule (reported only): per-CT median across the 3 individuals.
  ct_median_mean <- apply(ct_mean, 1, median)
  excluded_median <- names(ct_median_mean)[ct_median_mean > mean_excl_threshold]

  retained_indicator[retained_ct, cell_type] <- TRUE

  message(sprintf(
    "[excl] %-25s any-ind rule: %4d excluded / %5d retained | median-rule: %4d excluded",
    cell_type, length(excluded_ct), length(retained_ct),
    length(excluded_median)))

  writeLines(retained_ct,
             file.path(gene_list_dir, paste0(cell_type, "_retained.txt")))
  writeLines(excluded_ct,
             file.path(gene_list_dir, paste0(cell_type, "_excluded.txt")))

  # Per-cell-type filtered matrices (3 individual cols, any-individual retained).
  write.csv(res_var_mat[retained_ct, cols, drop = FALSE],
            file.path(filt_mat_dir, paste0(cell_type, "_dispersion_retained.csv")))
  write.csv(mean_mat[retained_ct, cols, drop = FALSE],
            file.path(filt_mat_dir, paste0(cell_type, "_mean_retained.csv")))

  summary_rows[[cell_type]] <- data.frame(
    cell_type = cell_type,
    n_genes_total = length(all_genes),
    n_excluded_any_individual = length(excluded_ct),
    n_retained_any_individual = length(retained_ct),
    n_excluded_median = length(excluded_median),
    stringsAsFactors = FALSE)
}

summary_df <- do.call(rbind, summary_rows)
rownames(summary_df) <- NULL
write.csv(summary_df, file.path(out_dir, "excluded_gene_summary.csv"),
          row.names = FALSE)

# Genes-x-CT boolean retained indicator (TRUE = retained under median rule),
# consumed by scripts 02 and 03.
write.csv(retained_indicator, file.path(out_dir, "retained_gene_matrix.csv"))

message("[done] wrote gene lists, filtered matrices, and summaries to ", out_dir)
print(summary_df)
