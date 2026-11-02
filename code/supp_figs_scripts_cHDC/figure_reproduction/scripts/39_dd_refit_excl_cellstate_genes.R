#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)              # normalizeQuantiles
  library(variancePartition)  # dream, makeContrastsDream, eBayes, topTable
  library(BiocParallel)       # SnowParam
  library(Seurat)             # cc.genes.updated.2019
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
fig_repro_dir <- file.path(repo_dir, "figure_reproduction")
dispersion_csv <- file.path(fig_repro_dir, "data", "res_var_excl_cluster24.csv")
out_dir <- file.path(fig_repro_dir, "outputs")
retained_matrix_csv <- file.path(out_dir, "retained_gene_matrix.csv")
baseline_fit_rds <- file.path(out_dir, "dd_refit", "dream_fit_retained.rds")
excl_dir <- file.path(out_dir, "dd_refit", "excl_cellstate_genes")

# Identical to 03_rerun_dd_rank_analyses.R.
ribo_pattern <- "^RP"
dd_adjp_thresh <- 0.05
# ANALYSIS_OK[model-formula]: dream mixed-model design formula identical to 03 Part B, no numeric weights.
dream_formula <- ~ 0 + cell_type
n_dream_cores <- 4
min_universe_genes <- 100

individuals <- c("NA19099", "NA19093", "NA18522")
cell_types <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)
n_expected_cols <- length(individuals) * length(cell_types)  # 18

# --- Gene sets to exclude ----------------------------------------------------
cc_genes <- Seurat::cc.genes.updated.2019
cell_cycle <- unique(c(cc_genes$s.genes, cc_genes$g2m.genes))
stopifnot(length(cell_cycle) == 97L)

# ANALYSIS_OK[hardcoded-gene-list]: literal HGNC symbols copied from marker_geneset_dispersion_sliding_window/scripts/02 (erik_celltype_markers).
erik_celltype_markers <- c(
  "MYH6", "NKX2-5", "TNNT2", "ACTC1", "RYR2",      # Cardiomyocytes
  "TWIST1", "PRRX2", "FOXC1",                      # Mesenchymal_cells
  "TCF21", "ITGA8", "PDGFRA", "COL1A1", "COL3A1",  # Cardiac_fibroblasts
  "POSTN", "DCN", "LUM", "PDPN",                   # Activated_fibroblasts
  "WFDC2", "GRHL2", "CLDN6",                       # Epithelial_cells
  "HHEX", "ALDH1A1", "ONECUT1"                     # Foregut_endodermal_cells
)
stopifnot(length(erik_celltype_markers) == 23L)
# ANALYSIS_OK[hardcoded-gene-list]: literal HGNC symbols copied from marker_geneset_dispersion_enrichment/scripts/01 (qc_plots_cHDC.Rmd DotPlot features).
celltype_marker_panel <- c(
  "POU5F1", "L1TD1", "SOX2", "CRABP1", "HHEX", "FOXA2", "AFP",
  "KDR", "CDH5", "PECAM1", "NFATC1", "TBX5", "TBX20",
  "WT1", "TBX18", "BMP4", "COL3A1", "LUM", "POSTN",
  "EYA1", "ISL1", "FGF10", "TNNT2", "ACTN2", "RYR2", "CACNA1C", "TNNI3"
)
stopifnot(length(celltype_marker_panel) == 27L)
celltype_markers <- union(erik_celltype_markers, celltype_marker_panel)

# ANALYSIS_OK[hardcoded-gene-list]: literal HGNC symbols copied from marker_geneset_dispersion_enrichment/scripts/01 (Brendan's chamber panel).
chamber_specificity <- c("MYL2", "MYL7", "HEY2", "IRX4", "TBX5", "NR2F2")

gene_sets <- list(
  cell_cycle = cell_cycle,
  celltype_markers = celltype_markers,
  chamber_specificity = chamber_specificity
)
excluded_genes <- unique(unlist(gene_sets, use.names = FALSE))
message(sprintf("[sets] cell_cycle=%d celltype_markers=%d chamber_specificity=%d; union=%d",
                length(cell_cycle), length(celltype_markers),
                length(chamber_specificity), length(excluded_genes)))

# --------------------------------- LOAD --------------------------------------
dir.create(excl_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(dispersion_csv)) stop("Missing dispersion matrix: ", dispersion_csv)
res_var_mat <- as.matrix(read.csv(dispersion_csv, row.names = 1, check.names = FALSE))
if (ncol(res_var_mat) != n_expected_cols) {
  stop("Dispersion matrix has ", ncol(res_var_mat), " cols; expected ", n_expected_cols)
}
# ANALYSIS_OK[filtering]: drop ^RP ribosomal genes, identical to 03; count logged below
res_var_mat <- res_var_mat[!grepl(ribo_pattern, rownames(res_var_mat)), , drop = FALSE]
all_genes <- rownames(res_var_mat)
message(sprintf("[load] %d non-ribo genes in dispersion matrix", length(all_genes)))

if (!file.exists(retained_matrix_csv)) stop("Missing ", retained_matrix_csv, " -- run 01 first")
retained_ind <- as.matrix(read.csv(retained_matrix_csv, row.names = 1, check.names = FALSE))
# ANALYSIS_OK[column-reorder]: select cell-type columns in canonical order; no rows dropped
retained_ind <- retained_ind[, cell_types, drop = FALSE]
if (!all(all_genes %in% rownames(retained_ind))) {
  stop("retained_gene_matrix.csv is missing ", sum(!all_genes %in% rownames(retained_ind)),
       " non-ribo genes from the dispersion matrix")
}
retained_ind <- retained_ind[all_genes, , drop = FALSE]
union_retained_genes <- all_genes[rowSums(retained_ind == 1) == length(cell_types)]

if (!file.exists(baseline_fit_rds)) stop("Missing ", baseline_fit_rds, " -- run 03 first")
fit_retained <- readRDS(baseline_fit_rds)
if (!setequal(rownames(fit_retained), union_retained_genes)) {
  stop("dream_fit_retained.rds genes do not match the union-retained set recomputed ",
       "from retained_gene_matrix.csv (", nrow(fit_retained), " vs ",
       length(union_retained_genes), ") -- 01/03 outputs are out of sync")
}
message(sprintf("[load] baseline union-retained universe: %d genes", length(union_retained_genes)))

# --- Excluded-gene table -----------------------------------------------------
excl_tab <- data.frame(
  gene = excluded_genes,
  in_cell_cycle = excluded_genes %in% cell_cycle,
  in_celltype_markers = excluded_genes %in% celltype_markers,
  in_chamber_specificity = excluded_genes %in% chamber_specificity,
  in_dispersion_matrix = excluded_genes %in% all_genes,
  in_union_retained_universe = excluded_genes %in% union_retained_genes,
  stringsAsFactors = FALSE)
write.csv(excl_tab, file.path(excl_dir, "excluded_genes.csv"), row.names = FALSE)
for (s in names(gene_sets)) {
  message(sprintf("[sets] %-20s %3d genes, %3d in union-retained universe",
                  s, length(gene_sets[[s]]),
                  sum(gene_sets[[s]] %in% union_retained_genes)))
}

# ANALYSIS_OK[filtering]: the requested exclusion itself; removed count logged and every gene listed in excluded_genes.csv
excl_universe <- setdiff(union_retained_genes, excluded_genes)
n_removed <- length(union_retained_genes) - length(excl_universe)
message(sprintf("[ddX] removed %d genes from the universe: %d -> %d",
                n_removed, length(union_retained_genes), length(excl_universe)))
if (length(excl_universe) < min_universe_genes) stop("Excluded universe too small")

# ------------------------------- DREAM FIT -----------------------------------
# Metadata and contrasts built exactly as in 03 Part B.
sample_cols <- colnames(res_var_mat)
col_ind <- sub("_.*$", "", sample_cols)
col_ct  <- sub("^[^_]*_", "", sample_cols)
stopifnot(all(col_ind %in% individuals), all(col_ct %in% cell_types))
metadata <- data.frame(
  cell_type = factor(col_ct, levels = cell_types),
  cell_line = factor(col_ind, levels = individuals),
  row.names = sample_cols)

pairs <- combn(cell_types, 2, simplify = FALSE)
contrast_matrix <- setNames(
  vapply(pairs, function(x) paste0("cell_type", x[2], " - ", "cell_type", x[1]),
         character(1)),
  vapply(pairs, function(x) paste0("compare_", x[2], "_", x[1]), character(1)))
stopifnot(all(names(contrast_matrix) %in% colnames(coef(fit_retained))))

mat <- normalizeQuantiles(res_var_mat[excl_universe, sample_cols, drop = FALSE])
contrast_dream <- makeContrastsDream(dream_formula, metadata, contrasts = contrast_matrix)
param <- SnowParam(n_dream_cores, type = "SOCK", progressbar = FALSE, exportglobals = FALSE)
message(sprintf("[ddX] dream fit: %d genes x %d samples, %d contrasts",
                length(excl_universe), length(sample_cols), length(contrast_matrix)))
fit_excl <- eBayes(dream(mat, dream_formula, metadata, contrast_dream, BPPARAM = param))
saveRDS(fit_excl, file.path(excl_dir, "dream_fit_excl_cellstate.rds"))

# ------------------------------- COMPARE -------------------------------------
# Lost/gained are computed against the genes eligible in BOTH fits (the
# excluded universe); baseline DD calls among the removed genes are counted
# separately, per gene set.
cmp_rows <- list()
for (cn in names(contrast_matrix)) {
  tt_r <- topTable(fit_retained, coef = cn, number = Inf, sort.by = "none")
  tt_x <- topTable(fit_excl, coef = cn, number = Inf, sort.by = "none")
  dd_r <- rownames(tt_r)[!is.na(tt_r$adj.P.Val) & tt_r$adj.P.Val <= dd_adjp_thresh]
  dd_x <- rownames(tt_x)[!is.na(tt_x$adj.P.Val) & tt_x$adj.P.Val <= dd_adjp_thresh]
  dd_r_kept <- intersect(dd_r, excl_universe)
  dd_r_removed <- setdiff(dd_r, excl_universe)
  cmp_rows[[cn]] <- data.frame(
    contrast = cn,
    n_dd_retained = length(dd_r),
    n_dd_retained_removed = length(dd_r_removed),
    n_dd_retained_removed_cell_cycle = sum(dd_r_removed %in% cell_cycle),
    n_dd_retained_removed_celltype_markers = sum(dd_r_removed %in% celltype_markers),
    n_dd_retained_removed_chamber = sum(dd_r_removed %in% chamber_specificity),
    n_dd_retained_in_excl_universe = length(dd_r_kept),
    n_dd_excl = length(dd_x),
    n_dd_lost = length(setdiff(dd_r_kept, dd_x)),
    n_dd_gained = length(setdiff(dd_x, dd_r_kept)),
    spearman_logFC = cor(tt_r[excl_universe, "logFC"], tt_x[excl_universe, "logFC"],
                         method = "spearman"),
    spearman_t = cor(tt_r[excl_universe, "t"], tt_x[excl_universe, "t"], method = "spearman"),
    stringsAsFactors = FALSE)
}
cmp <- do.call(rbind, cmp_rows)
rownames(cmp) <- NULL
write.csv(cmp, file.path(excl_dir, "dd_refit_retained_vs_excl_cellstate.csv"), row.names = FALSE)
print(cmp)

# ---- reportable values (same format as 14) ----------------------------------
reportable <- data.frame(
  key = c("dd_excl_cellstate_n_genes_removed",
          "dd_excl_cellstate_n_genes_universe",
          "dd_excl_cellstate_total_dd_retained",
          "dd_excl_cellstate_total_dd_excl",
          "dd_excl_cellstate_total_dd_retained_removed",
          "dd_excl_cellstate_total_dd_lost",
          "dd_excl_cellstate_total_dd_gained",
          "dd_excl_cellstate_min_spearman_logFC"),
  value = c(n_removed, length(excl_universe),
            sum(cmp$n_dd_retained), sum(cmp$n_dd_excl),
            sum(cmp$n_dd_retained_removed), sum(cmp$n_dd_lost),
            sum(cmp$n_dd_gained), round(min(cmp$spearman_logFC), 6)),
  provenance = "figure_reproduction/scripts/39_dd_refit_excl_cellstate_genes.R")
write.csv(reportable, file.path(out_dir, "reportable_values_39.csv"), row.names = FALSE)

message("[done] outputs in ", excl_dir)
