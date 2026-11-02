#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(limma)              # normalizeQuantiles
  library(variancePartition)  # dream, makeContrastsDream, eBayes, topTable
  library(BiocParallel)       # SnowParam
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
fig_repro_dir <- file.path(repo_dir, "figure_reproduction")
data_dir <- file.path(fig_repro_dir, "data")
cluster_dir <- file.path(repo_dir, "cellranger_cluster-mode_trial")
pub_dir <- file.path(cluster_dir, "analysis", "first_pub_datasets")

dispersion_csv <- file.path(data_dir, "res_var_excl_cluster24.csv")

# Erik's eQTL results (external, read-only, unrelated to the Epithelial
# exclusion). eGenes are genes with a significant cis-eQTL in the Control
# condition for a given cell type.
# ANALYSIS_OK[external-path]: read-only path to a frozen collaborator results dir (Erik's map_eqtls/archive); identity pinned by the explicit path
eqtl_dir <- file.path(
  "/project/gilad/emcintire/chemo_project_2023/pipeline/map_eqtls",
  "archive", "no_cvrt_condox_only", "eqtl_results", "Control", "%s")
eqtl_file_tmpl <- "steadyState_%s_Control_AllChroms.txt"  # sprintf with cell type

orig_fisher_csv <- file.path(pub_dir, "egene_quintile_fisher_tests.csv")

out_dir <- file.path(fig_repro_dir, "outputs")
dd_dir <- file.path(out_dir, "dd_rank")
refit_dir <- file.path(out_dir, "dd_refit")
retained_matrix_csv <- file.path(out_dir, "retained_gene_matrix.csv")

ribo_pattern <- "^RP"
n_quantiles <- 5              # source line 1264
tolerance_frac <- 0.15        # source line 1262: tolerance = 11847 * 0.15
egene_pvalue_thresh <- 0.05   # source line 1217 (bh_adj_pvalue <= 0.05)
egene_min_count <- 49         # source line 1222 (cell types with >= 49 eGenes)
dd_adjp_thresh <- 0.05        # DD call: adj.P.Val <= 0.05

# Full dream DD re-fit (source: reformat_and_dream_de_dd_analysis.R lines 87-123)
# ANALYSIS_OK[model-formula]: this is a dream mixed-model design formula, not a
#   weighted composite score; the terms are model effects (cell_type fixed,
#   cell_line random) with no numeric weights, matching the published DD model.
dream_formula <- ~ 0 + cell_type  # 0+ absorbs intercept into cell_type
n_dream_cores <- 4
min_universe_genes <- 100     # sanity floor on the rank / re-fit gene-set size

individuals <- c("NA19099", "NA19093", "NA18522")
cell_types <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)
n_expected_cols <- length(individuals) * length(cell_types)  # 18

# ------------------------------- HELPERS -------------------------------------

source(file.path(fig_repro_dir, "scripts", "utils.R"))

# Load eGenes for one cell type (unique gene symbols passing the BH threshold).
# ANALYSIS_OK[local-helper]: used within this standalone pipeline-stage script (run via run.sh); not a shared-library export
load_egenes <- function(cell_type) {
  path <- file.path(sprintf(eqtl_dir, cell_type), sprintf(eqtl_file_tmpl, cell_type))
  if (!file.exists(path)) stop("Missing eQTL results file: ", path)
  tab <- readr::read_table(path, show_col_types = FALSE)
  req <- c("gene", "bh_adj_pvalue")
  if (!all(req %in% colnames(tab))) {
    stop("eQTL file ", basename(path), " missing columns: ",
         paste(setdiff(req, colnames(tab)), collapse = ", "))
  }
  unique(tab$gene[tab$bh_adj_pvalue <= egene_pvalue_thresh])
}

# --------------------------------- MAIN --------------------------------------

dir.create(dd_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load dispersion matrix, drop ^RP genes ----------------------------------
if (!file.exists(dispersion_csv)) stop("Missing dispersion matrix: ", dispersion_csv)
res_var_df <- read.csv(dispersion_csv, row.names = 1, check.names = FALSE)
if (ncol(res_var_df) != n_expected_cols) {
  stop("Dispersion matrix has ", ncol(res_var_df), " cols; expected ", n_expected_cols)
}
res_var_mat <- as.matrix(res_var_df)
# ANALYSIS_OK[filtering]: drop ^RP ribosomal genes (repo convention); retained non-ribo count logged immediately below
res_var_mat <- res_var_mat[!grepl(ribo_pattern, rownames(res_var_mat)), , drop = FALSE]
all_genes <- rownames(res_var_mat)
message(sprintf("[load] %d non-ribo genes in dispersion matrix", length(all_genes)))

# --- Load per-CT retained indicator from script 01 ---------------------------
if (!file.exists(retained_matrix_csv)) {
  stop("Missing retained_gene_matrix.csv — run 01_exclude_topexpressed_genes.R first")
}
retained_ind <- as.matrix(read.csv(retained_matrix_csv, row.names = 1, check.names = FALSE))
# ANALYSIS_OK[column-reorder]: reorder/select the cell-type columns to canonical order only; no genes/rows dropped, values unchanged
retained_ind <- retained_ind[, cell_types, drop = FALSE]

# =============================================================================
# PART A — eGene enrichment across dispersion quintiles (Fig 3A), on retained
# =============================================================================

# Load eGenes and keep only cell types with >= egene_min_count (source 1222).
egenes_by_ct <- list()
for (cell_type in cell_types) {
  eg <- load_egenes(cell_type)
  message(sprintf("[egene] %-25s %d eGenes (bh_adj_pvalue <= %.2f)",
                  cell_type, length(eg), egene_pvalue_thresh))
  if (length(eg) >= egene_min_count) egenes_by_ct[[cell_type]] <- eg
}
all_egenes <- unique(unlist(egenes_by_ct, use.names = FALSE))
message(sprintf("[egene] %d cell types included; %d unique eGenes total",
                length(egenes_by_ct), length(all_egenes)))

# Global-ranking universe.
# TODO(assumption): the global ranking (source lines 1271-1314) needs a per-CT
#   rank for every gene in every cell type with no NA. Because the exclusion is
#   per cell type, we restrict the ranking universe to genes RETAINED IN ALL 6
#   cell types (intersection). Confirm this is the intended universe for the
#   revised Fig 3A rather than, e.g., per-CT retained ranks reconciled some
#   other way.
ranking_genes <- all_genes[rowSums(retained_ind == 1) == length(cell_types)]
message(sprintf("[rankA] global-ranking universe (retained in all CTs): %d genes",
                length(ranking_genes)))
if (length(ranking_genes) < min_universe_genes) {
  stop("Global-ranking universe too small (", length(ranking_genes),
       ") — check the retained indicator")
}

# Per-CT median-dispersion ranks over the ranking universe (source 1277-1290).
ct_ranks <- matrix(NA_real_, nrow = length(ranking_genes), ncol = length(cell_types),
                   dimnames = list(ranking_genes, cell_types))
for (cell_type in cell_types) {
  ct_res_var <- res_var_mat[ranking_genes, ct_columns(cell_type), drop = FALSE]
  med <- apply(ct_res_var, 1, median)
  ct_ranks[, cell_type] <- rank(med, ties.method = "average")
}

# Keep genes whose per-CT ranks are all within `tolerance` of their median rank
# (source lines 1294-1300). tolerance recomputed against this universe's size.
tolerance <- tolerance_frac * nrow(ct_ranks)
avg_rank <- apply(ct_ranks, 1, median)
within_tol <- apply(abs(ct_ranks - avg_rank) <= tolerance, 1, all)
filtered_ranks <- ct_ranks[within_tol, , drop = FALSE]
message(sprintf("[rankA] tolerance = %.1f positions; %d genes within tolerance",
                tolerance, nrow(filtered_ranks)))

global_rank <- rank(apply(filtered_ranks, 1, median), ties.method = "average")
global_quintile <- cut(
  global_rank,
  breaks = quantile(global_rank, probs = seq(0, 1, length.out = n_quantiles + 1),
                    na.rm = TRUE),
  labels = FALSE, include.lowest = TRUE)
names(global_quintile) <- rownames(filtered_ranks)

# Per-quintile Fisher test for eGene enrichment (source lines 1320-1369).
background_genes <- names(global_quintile)
egenes_bg <- intersect(all_egenes, background_genes)
fisher_rows <- list()
for (i in seq_len(n_quantiles)) {
  quint_genes <- names(global_quintile)[global_quintile == i]
  a <- length(intersect(egenes_bg, quint_genes))
  b <- length(setdiff(egenes_bg, quint_genes))
  c <- length(setdiff(quint_genes, egenes_bg))
  d <- length(setdiff(background_genes, union(egenes_bg, quint_genes)))
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
  fisher_rows[[i]] <- data.frame(
    quintile = i, egene_count = a, p.value = ft$p.value,
    OR = unname(ft$estimate),
    ci_lower = ft$conf.int[1], ci_upper = ft$conf.int[2],
    log2OR = log2(unname(ft$estimate)),
    log2_ci_lower = log2(ft$conf.int[1]), log2_ci_upper = log2(ft$conf.int[2]),
    stringsAsFactors = FALSE)
}
fisher_df <- do.call(rbind, fisher_rows)
fisher_df$padj <- p.adjust(fisher_df$p.value, method = "BH")  # source line 1369
write.csv(fisher_df, file.path(dd_dir, "egene_quintile_fisher_retained.csv"),
          row.names = FALSE)
message("[rankA] wrote retained-gene eGene-quintile Fisher table")
print(fisher_df[, c("quintile", "egene_count", "p.value", "log2OR", "padj")])

# Side-by-side with the manuscript's stored result, if available.
if (file.exists(orig_fisher_csv)) {
  orig_fisher <- read.csv(orig_fisher_csv, row.names = 1)
  cmp <- data.frame(
    quintile = fisher_df$quintile,
    log2OR_original = orig_fisher$log2OR[match(fisher_df$quintile, orig_fisher$quintile)],
    log2OR_retained = fisher_df$log2OR,
    p_original = orig_fisher$p.value[match(fisher_df$quintile, orig_fisher$quintile)],
    p_retained = fisher_df$p.value)
  write.csv(cmp, file.path(dd_dir, "egene_quintile_fisher_comparison.csv"),
            row.names = FALSE)
  message("[rankA] wrote original-vs-retained Fisher comparison")
} else {
  message("[rankA] NOTE: ", basename(orig_fisher_csv),
          " not found; skipping original comparison")
}

# =============================================================================
# PART B — FULL dream DD re-fit (joint, union exclusion): full vs retained
# =============================================================================
dir.create(refit_dir, recursive = TRUE, showWarnings = FALSE)

# Union-retained gene set: genes excluded in NO cell type (retained in all 6).
# This equals Part A's ranking universe; recomputed here with a clear name.
union_retained_genes <- all_genes[rowSums(retained_ind == 1) == length(cell_types)]
message(sprintf(
  "[ddB] union-retained genes (excluded in no cell type): %d of %d non-ribo",
  length(union_retained_genes), length(all_genes)))
if (length(union_retained_genes) < min_universe_genes) {
  stop("Union-retained gene set too small (", length(union_retained_genes), ")")
}

# Build sample metadata by PARSING column names ({IND}_{CT}); do NOT assume
# column order. Cell-type names contain underscores, so split on the FIRST '_'
# only (individual IDs are NA##### with no underscore).
sample_cols <- colnames(res_var_mat)
col_ind <- sub("_.*$", "", sample_cols)
col_ct  <- sub("^[^_]*_", "", sample_cols)
stopifnot(all(col_ind %in% individuals), all(col_ct %in% cell_types),
          length(sample_cols) == n_expected_cols)
metadata <- data.frame(
  cell_type = factor(col_ct, levels = cell_types),
  cell_line = factor(col_ind, levels = individuals),
  row.names = sample_cols)

# 15 pairwise contrasts named compare_{B}_{A} (source lines 93-107).
pairs <- combn(cell_types, 2, simplify = FALSE)
contrast_matrix <- setNames(
  vapply(pairs, function(x) paste0("cell_type", x[2], " - ", "cell_type", x[1]),
         character(1)),
  vapply(pairs, function(x) paste0("compare_", x[2], "_", x[1]), character(1)))

# One dream + eBayes fit on a given gene set; identical pipeline for full/retained
# so the only difference between the two fits is which genes are included.
# ANALYSIS_OK[local-helper]: used within this standalone pipeline-stage script (run via run.sh); not a shared-library export
fit_dream_dd <- function(gene_set, tag) {
  mat <- normalizeQuantiles(res_var_mat[gene_set, sample_cols, drop = FALSE])
  contrast_dream <- makeContrastsDream(dream_formula, metadata,
                                       contrasts = contrast_matrix)
  param <- SnowParam(n_dream_cores, type = "SOCK",
                     progressbar = FALSE, exportglobals = FALSE)
  message(sprintf("[ddB] dream fit '%s': %d genes x %d samples, %d contrasts",
                  tag, length(gene_set), length(sample_cols),
                  length(contrast_matrix)))
  fit <- dream(mat, dream_formula, metadata, contrast_dream, BPPARAM = param)
  eBayes(fit)
}

# Per-contrast DD gene sets (adj.P.Val <= dd_adjp_thresh) from a fitted object.
# ANALYSIS_OK[local-helper]: used within this standalone pipeline-stage script (run via run.sh); not a shared-library export
dd_sets_from_fit <- function(fit) {
  setNames(lapply(names(contrast_matrix), function(cn) {
    tt <- topTable(fit, coef = cn, number = Inf, sort.by = "none")
    rownames(tt)[!is.na(tt$adj.P.Val) & tt$adj.P.Val <= dd_adjp_thresh]
  }), names(contrast_matrix))
}

fit_full     <- fit_dream_dd(all_genes, "full")
fit_retained <- fit_dream_dd(union_retained_genes, "retained")
saveRDS(fit_full,     file.path(refit_dir, "dream_fit_full.rds"))
saveRDS(fit_retained, file.path(refit_dir, "dream_fit_retained.rds"))

dd_full     <- dd_sets_from_fit(fit_full)
dd_retained <- dd_sets_from_fit(fit_retained)

# Compare per contrast. A full-fit DD gene can only persist in the retained fit
# if it was in the retained universe at all, so lost/gained are computed against
# that eligible set.
cmp_rows <- list()
for (cn in names(contrast_matrix)) {
  f <- dd_full[[cn]]
  r <- dd_retained[[cn]]
  f_in_univ <- intersect(f, union_retained_genes)
  cmp_rows[[cn]] <- data.frame(
    contrast = cn,
    n_dd_full = length(f),
    n_dd_full_in_retained_universe = length(f_in_univ),
    n_dd_retained = length(r),
    n_dd_lost = length(setdiff(f_in_univ, r)),
    n_dd_gained = length(setdiff(r, f)),
    stringsAsFactors = FALSE)
}
dd_refit_compare <- do.call(rbind, cmp_rows)
rownames(dd_refit_compare) <- NULL
write.csv(dd_refit_compare, file.path(refit_dir, "dd_refit_full_vs_retained.csv"),
          row.names = FALSE)

message("[done] Part A Fisher tables in ", dd_dir,
        "; Part B dream re-fit results in ", refit_dir)
print(dd_refit_compare)
