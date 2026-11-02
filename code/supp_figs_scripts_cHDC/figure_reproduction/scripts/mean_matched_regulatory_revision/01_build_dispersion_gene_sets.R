#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

args        <- commandArgs(trailingOnly = FALSE)
script_dir  <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
source(file.path(script_dir, "config.R"))
gs_dir <- ensure_dir(file.path(out_dir, "gene_sets"))

# ---- Load dispersion matrix + ribo filter (Rmd 50-55) -----------------------
assert_that(file.exists(dispersion_matrix), paste("missing:", dispersion_matrix))
res_var_df <- read.csv(dispersion_matrix, row.names = 1)
log_step("loaded dispersion matrix (genes)", nrow(res_var_df))
assert_that(ncol(res_var_df) == length(cell_types) * length(individuals),
            sprintf("expected %d columns, got %d",
                    length(cell_types) * length(individuals), ncol(res_var_df)))

ribo_genes_full <- rownames(res_var_df)[grepl(ribo_pattern, rownames(res_var_df))]
# ANALYSIS_OK[filtering]: drop ribosomal (^RP) genes per the manuscript convention;
# the number removed and the remaining gene count are logged immediately below.
res_var_df <- res_var_df %>% filter(!rownames(.) %in% ribo_genes_full)
log_step(sprintf("removed %d ribo (%s) genes; remaining",
                 length(ribo_genes_full), ribo_pattern), nrow(res_var_df))
assert_that(nrow(res_var_df) > 0, "no genes left after ribo filter")

# ---- -6 exclusion indicator (figure_reproduction script 01) ------------------
assert_that(file.exists(retained_matrix), paste("missing:", retained_matrix))
retained_ind <- as.matrix(read.csv(retained_matrix, row.names = 1, check.names = FALSE))
assert_that(setequal(colnames(retained_ind), cell_types),
            "retained_gene_matrix.csv columns do not match cell_types")
assert_that(setequal(rownames(retained_ind), rownames(res_var_df)),
            "retained_gene_matrix.csv genes differ from the ribo-filtered dispersion matrix")

# ---- Per-cell-type median-rank quintile table (Rmd 65-107) ------------------
# Ranked over each cell type's retained genes only; excluded genes are NA.
ct_ranks_df <- as.data.frame(matrix(NA, nrow = nrow(res_var_df),
                                    ncol = length(cell_types)))
colnames(ct_ranks_df) <- cell_types
rownames(ct_ranks_df) <- rownames(res_var_df)
quantile_bins <- ct_ranks_df

for (cell_type in cell_types) {
  cols <- paste0(individuals, "_", cell_type)
  assert_that(all(cols %in% colnames(res_var_df)),
              sprintf("missing dispersion columns for %s", cell_type))
  ct_retained    <- rownames(retained_ind)[retained_ind[, cell_type] == 1]
  log_step(sprintf("%s: ranking retained genes", cell_type), length(ct_retained))
  ct_res_var     <- res_var_df[ct_retained, cols, drop = FALSE]
  median_res_var <- apply(ct_res_var, 1, median)
  ranks          <- rank(median_res_var, ties.method = "average")
  ct_ranks_df[ct_retained, cell_type] <- ranks
  quantile_bins[ct_retained, cell_type] <- cut(
    ranks,
    breaks = quantile(ranks, probs = seq(0, 1, length.out = n_quantiles + 1)),
    labels = FALSE, include.lowest = TRUE)
}

# Restrict to genes retained in all 6 cell types (the dream DE/DD fit universe).
union_retained <- rownames(retained_ind)[rowSums(retained_ind == 1) == length(cell_types)]
# ANALYSIS_OK[filtering]: -6 top-expressed-gene exclusion (repo convention, figure_reproduction script 01); count logged below
var_quintile_dfs <- quantile_bins[union_retained, , drop = FALSE]
log_step("genes retained in all cell types (quintile-table universe)", nrow(var_quintile_dfs))
assert_that(!anyNA(var_quintile_dfs), "NA quintile among all-cell-type-retained genes")
assert_that(all(unlist(var_quintile_dfs) %in% seq_len(n_quantiles)),
            "quintile bins outside 1..n_quantiles")
log_step("built quintile table (genes x cell types)", nrow(var_quintile_dfs))

# ---- Derive the three consistent dispersion gene sets (Rmd 1433-1441 / 2777)-
var_first_quint_genes <- var_quintile_dfs %>%
  dplyr::filter(rowSums(. == lowest_quintile) >= quintile_membership_min) %>% rownames()
var_fifth_quint_genes <- var_quintile_dfs %>%
  dplyr::filter(rowSums(. == highest_quintile) >= quintile_membership_min) %>% rownames()
never_first_fifth_genes <- var_quintile_dfs %>%
  dplyr::filter(rowSums(!. == lowest_quintile) >= quintile_membership_min &
                rowSums(!. == highest_quintile) >= quintile_membership_min) %>% rownames()

log_step("var_first_quint_genes (consistently low)",  length(var_first_quint_genes))
log_step("var_fifth_quint_genes (consistently high)", length(var_fifth_quint_genes))
log_step("never_first_fifth_genes (neither)",         length(never_first_fifth_genes))

# Sanity: low/high sets must be disjoint; both must be non-empty.
assert_that(length(intersect(var_first_quint_genes, var_fifth_quint_genes)) == 0,
            "low and high dispersion gene sets overlap")
assert_that(length(var_first_quint_genes) > 0 && length(var_fifth_quint_genes) > 0,
            "a consistent dispersion gene set is empty — check membership threshold")

# ---- Save -------------------------------------------------------------------
saveRDS(list(var_quintile_dfs        = var_quintile_dfs,
             var_first_quint_genes   = var_first_quint_genes,
             var_fifth_quint_genes   = var_fifth_quint_genes,
             never_first_fifth_genes = never_first_fifth_genes,
             ribo_genes_removed      = ribo_genes_full,
             # Raw per-cell-type dispersion rank (1..N per cell type), kept
             # alongside the quintile-bin table -- added 2026-09-03 for 06's
             # per-cell-type (not consistent-across-cell-types) gene-set
             # definition (Fig3C, Rmd 1536-1551).
             ct_ranks_df              = ct_ranks_df),
        file.path(gs_dir, "dispersion_gene_sets.rds"))
write.csv(var_quintile_dfs %>% rownames_to_column("gene"),
          file.path(gs_dir, "var_quintile_table.csv"), row.names = FALSE)
writeLines(var_first_quint_genes,   file.path(gs_dir, "var_first_quint_genes.txt"))
writeLines(var_fifth_quint_genes,   file.path(gs_dir, "var_fifth_quint_genes.txt"))
writeLines(never_first_fifth_genes, file.path(gs_dir, "never_first_fifth_genes.txt"))
log_step(paste("wrote gene sets to", gs_dir))
