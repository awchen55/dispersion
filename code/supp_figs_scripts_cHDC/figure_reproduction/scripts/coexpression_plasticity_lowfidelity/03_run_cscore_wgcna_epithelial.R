#!/usr/bin/env Rscript

stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(Seurat)
  library(CSCORE)
  library(WGCNA)
  library(dynamicTreeCut)
})

IN_DIR  <- "figure_reproduction/outputs/coexpression_plasticity_lowfidelity/cscore_inputs"
OUT_DIR <- "figure_reproduction/outputs/coexpression_plasticity_lowfidelity/cscore_modules"

CELL_TYPES  <- c("Epithelial_cells")
INDIVIDUALS <- c("NA19099", "NA19093", "NA18522")

P_ADJUST_THRESH   <- 0.05
ADJACENCY_POWER   <- 1
DEEP_SPLIT        <- 2
MIN_CLUSTER_SIZE  <- 10

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) message("[03_run_cscore_wgcna_epithelial] ", ...)

if (!nzchar(Sys.getenv("CELL_TYPE")) || !nzchar(Sys.getenv("INDIVIDUAL"))) {
  stop("CELL_TYPE and INDIVIDUAL environment variables must both be set -- aborting (no default).")
}
CELL_TYPE  <- Sys.getenv("CELL_TYPE")
INDIVIDUAL <- Sys.getenv("INDIVIDUAL")
if (!(CELL_TYPE %in% CELL_TYPES)) {
  stop("CELL_TYPE=", CELL_TYPE, " not one of the expected cell types: ",
       paste(CELL_TYPES, collapse = ", "))
}
if (!(INDIVIDUAL %in% INDIVIDUALS)) {
  stop("INDIVIDUAL=", INDIVIDUAL, " not one of the 3 expected individuals: ",
       paste(INDIVIDUALS, collapse = ", "))
}
sample_label <- paste0(CELL_TYPE, "_", INDIVIDUAL)
log_msg("Sample: ", sample_label)

in_rds <- file.path(IN_DIR, paste0(sample_label, ".rds"))
if (!file.exists(in_rds)) {
  stop("Missing input: ", in_rds, ". Run 02_prepare_cscore_inputs_epithelial.R first.")
}

seurat_obj <- readRDS(in_rds)
stopifnot("expected a Seurat object" = inherits(seurat_obj, "Seurat"))
genes_selected <- rownames(seurat_obj)
n_genes <- length(genes_selected)
log_msg("Loaded ", ncol(seurat_obj), " cells x ", n_genes, " genes.")

if (!("nCount_RNA" %in% colnames(seurat_obj[[]]))) {
  stop("nCount_RNA column missing on the per-sample object -- script 02 should ",
       "have set this to the full-transcriptome UMI total per cell.")
}
if (any(seurat_obj$nCount_RNA <= 0)) {
  stop("Found cell(s) with non-positive nCount_RNA -- cannot run CSCORE.")
}

t0 <- Sys.time()
CSCORE_result <- CSCORE(seurat_obj, genes = genes_selected)
log_msg("CSCORE() done in ", round(difftime(Sys.time(), t0, units = "mins"), 2), " min.")

stopifnot(
  "CSCORE_result missing $est"     = "est" %in% names(CSCORE_result),
  "CSCORE_result missing $p_value" = "p_value" %in% names(CSCORE_result),
  "est matrix has wrong dimensions"     = all(dim(CSCORE_result$est) == n_genes),
  "p_value matrix has wrong dimensions" = all(dim(CSCORE_result$p_value) == n_genes)
)

est <- CSCORE::post_process_est(CSCORE_result$est)

p_value <- CSCORE_result$p_value
p_matrix_BH <- matrix(0, n_genes, n_genes)
p_matrix_BH[upper.tri(p_matrix_BH)] <- stats::p.adjust(
  p_value[upper.tri(p_value)], method = "BH")
p_matrix_BH <- p_matrix_BH + t(p_matrix_BH)

coexp <- est
coexp[p_matrix_BH > P_ADJUST_THRESH] <- 0
n_sig_pairs <- sum(p_matrix_BH[upper.tri(p_matrix_BH)] <= P_ADJUST_THRESH)
n_total_pairs <- sum(upper.tri(p_matrix_BH))
log_msg("Significant co-expression pairs (BH-adjusted p <= ", P_ADJUST_THRESH, "): ",
        n_sig_pairs, " of ", n_total_pairs,
        " (", round(100 * n_sig_pairs / n_total_pairs, 2), "%).")

t1 <- Sys.time()
adj <- WGCNA::adjacency.fromSimilarity(abs(coexp), power = ADJACENCY_POWER)
TOM <- WGCNA::TOMsimilarity(adj)
dissTOM <- 1 - TOM
rownames(dissTOM) <- colnames(dissTOM) <- genes_selected
hclust_dist <- stats::hclust(stats::as.dist(dissTOM), method = "average")
memb <- dynamicTreeCut::cutreeDynamic(
  dendro = hclust_dist,
  distM = dissTOM,
  deepSplit = DEEP_SPLIT,
  pamRespectsDendro = FALSE,
  minClusterSize = MIN_CLUSTER_SIZE
)
names(memb) <- genes_selected
log_msg("WGCNA/dynamicTreeCut done in ", round(difftime(Sys.time(), t1, units = "mins"), 2), " min.")

stopifnot("module membership vector length mismatch" = length(memb) == n_genes)

n_modules <- length(unique(memb[memb != 0]))
n_in_module <- sum(memb != 0)
log_msg(n_modules, " module(s) detected (excluding label 0/unassigned); ",
        n_in_module, " of ", n_genes, " genes (",
        round(100 * n_in_module / n_genes, 1), "%) assigned to a module.")

module_df <- data.frame(
  gene = genes_selected,
  module = as.integer(memb),
  in_module = memb != 0,
  cell_type = CELL_TYPE,
  individual = INDIVIDUAL,
  row.names = NULL
)
out_csv <- file.path(OUT_DIR, paste0(sample_label, "_module_membership.csv"))
utils::write.csv(module_df, out_csv, row.names = FALSE)
log_msg("Wrote: ", out_csv)

out_coexp_rds <- file.path(OUT_DIR, paste0(sample_label, "_coexp_thresholded.rds"))
saveRDS(coexp, out_coexp_rds)
log_msg("Wrote: ", out_coexp_rds)

log_msg("DONE for ", sample_label, ".")
