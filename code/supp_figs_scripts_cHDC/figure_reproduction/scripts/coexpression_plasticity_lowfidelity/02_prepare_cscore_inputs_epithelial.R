#!/usr/bin/env Rscript

stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# ---------------------------------------------------------------------------
# NAMED CONSTANTS
# ---------------------------------------------------------------------------
RAW_COUNTS_SEURAT_RDS <- "figure_reproduction/outputs/coexpression_plasticity_lowfidelity/epithelial_excl_cluster24_seurat_subset.rds"

CELLTYPE_COL   <- "new.predicted.id"
INDIVIDUAL_COL <- "vireo.individual"
UMI_COL        <- "nCount_RNA"

CELL_TYPES <- c("Epithelial_cells")
INDIVIDUALS <- c("NA19099", "NA19093", "NA18522")

RIBO_PATTERN <- "^RP"
MITO_PATTERN <- "^MT-"

N_GENES_SELECTED <- 5000L

OUT_DIR <- "figure_reproduction/outputs/coexpression_plasticity_lowfidelity/cscore_inputs"
LOG_FILE <- file.path(OUT_DIR, "prepare_inputs_log.txt")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

log_lines <- character(0)
log_msg <- function(...) {
  line <- paste0(...)
  message(line)
  log_lines[[length(log_lines) + 1L]] <<- line
}

assert_raw_counts <- function(mat, label) {
  if (!(nrow(mat) > 0 && ncol(mat) > 0)) {
    stop(label, ": counts matrix is empty")
  }
  set.seed(1L)
  samp_idx <- sample(length(mat@x), size = min(2e6L, length(mat@x)))
  vals <- mat@x[samp_idx]
  if (!all(vals == round(vals))) {
    stop(label, ": counts contain non-integer values -- not raw counts")
  }
  if (!all(vals >= 0)) {
    stop(label, ": counts contain negative values -- not raw counts")
  }
}

stopifnot("pre-subsetted Epithelial excl-cluster24 Seurat RDS not found -- run script 00 first" =
            file.exists(RAW_COUNTS_SEURAT_RDS))

# ---------------------------------------------------------------------------
# LOAD (much smaller than the original ~23GB full object -- already
# restricted to Epithelial_cells, sub-cluster 24 excluded, by script 00)
# ---------------------------------------------------------------------------
log_msg("Loading pre-subsetted Seurat object: ", RAW_COUNTS_SEURAT_RDS)
obj <- readRDS(RAW_COUNTS_SEURAT_RDS)
stopifnot("expected a Seurat object" = inherits(obj, "Seurat"))

md <- obj[[]]
stopifnot(
  "cell type column missing"  = CELLTYPE_COL %in% colnames(md),
  "individual column missing" = INDIVIDUAL_COL %in% colnames(md)
)
log_msg("Total cells: ", ncol(obj), " | total genes: ", nrow(obj))
stopifnot(
  "expected this pre-subsetted object to contain only Epithelial_cells" =
    setequal(unique(as.character(md[[CELLTYPE_COL]])), CELL_TYPES)
)

has_umi_col <- UMI_COL %in% colnames(md)
if (!has_umi_col) {
  log_msg("  '", UMI_COL, "' column not found in metadata -- will compute ",
          "per-cell UMI totals directly from the full (pre-filter) counts ",
          "matrix instead (matches CS-CORE's own seq_depth fallback logic).")
}

counts_all <- GetAssayData(obj, assay = "RNA", layer = "counts")
assert_raw_counts(counts_all, "full object")
log_msg("Raw-count assertion passed on full object.")

# ---------------------------------------------------------------------------
# PER-SAMPLE (Epithelial_cells x individual): subset, gene-filter/select, save
# ---------------------------------------------------------------------------
for (ct in CELL_TYPES) {
  for (ind in INDIVIDUALS) {
    sample_label <- paste0(ct, "_", ind)
    log_msg("\n=== ", sample_label, " ===")

    cell_mask <- (as.character(md[[CELLTYPE_COL]]) == ct) &
      (as.character(md[[INDIVIDUAL_COL]]) == ind)
    cells_this <- rownames(md)[cell_mask]
    log_msg("  Cells: ", length(cells_this))
    if (length(cells_this) == 0L) {
      log_msg("  SKIP: no cells for this (cell type, individual) combination.")
      next
    }

    if (has_umi_col) {
      full_umi <- md[cells_this, UMI_COL]
    } else {
      full_umi <- Matrix::colSums(counts_all[, cells_this, drop = FALSE])
    }
    names(full_umi) <- cells_this
    stopifnot("UMI-per-cell vector length mismatch" = length(full_umi) == length(cells_this))
    if (any(full_umi <= 0)) {
      stop(sample_label, ": found cell(s) with zero/negative total UMI -- ",
           "cannot compute mean-expression ranking.")
    }
    log_msg("  Full-transcriptome UMI/cell: median=", stats::median(full_umi),
            ", range=[", min(full_umi), ", ", max(full_umi), "]")

    cts_mat <- counts_all[, cells_this, drop = FALSE]

    ribo <- grepl(RIBO_PATTERN, rownames(cts_mat))
    mito <- grepl(MITO_PATTERN, rownames(cts_mat))
    log_msg("  Dropping ", sum(ribo), " ribosomal + ", sum(mito),
            " mitochondrial genes (of ", nrow(cts_mat), ").")
    keep_gene <- !ribo & !mito
    cts_filtered <- cts_mat[keep_gene, , drop = FALSE]
    log_msg("  Genes after ribo/mito filter: ", nrow(cts_filtered))

    norm_mat <- cts_filtered %*% Matrix::Diagonal(x = 1 / full_umi)
    mean_exp <- Matrix::rowMeans(norm_mat)

    n_select <- min(N_GENES_SELECTED, length(mean_exp))
    if (n_select < N_GENES_SELECTED) {
      log_msg("  WARNING: only ", length(mean_exp), " genes available after ",
              "ribo/mito filtering -- using all of them (< ",
              N_GENES_SELECTED, " requested).")
    }
    genes_selected <- names(sort(mean_exp, decreasing = TRUE))[seq_len(n_select)]
    log_msg("  Selected top ", n_select, " genes by mean expression ",
            "(selected range: [", signif(min(mean_exp[genes_selected]), 4), ", ",
            signif(max(mean_exp[genes_selected]), 4), "]).")

    cts_final <- cts_filtered[genes_selected, , drop = FALSE]
    assert_raw_counts(cts_final, sample_label)

    seurat_sample <- CreateSeuratObject(counts = cts_final, project = sample_label)
    stopifnot(
      "cell count mismatch after CreateSeuratObject" = ncol(seurat_sample) == length(cells_this),
      "gene count mismatch after CreateSeuratObject" = nrow(seurat_sample) == n_select
    )
    new_ncount <- full_umi[colnames(seurat_sample)]
    if (anyNA(new_ncount)) {
      stop(sample_label, ": nCount_RNA overwrite produced NA -- full_umi is ",
           "not correctly named by cell barcode. Fix the name-alignment bug ",
           "before proceeding (do not silently propagate NAs downstream).")
    }
    seurat_sample$nCount_RNA <- new_ncount
    seurat_sample$cell_type  <- ct
    seurat_sample$individual <- ind

    out_rds <- file.path(OUT_DIR, paste0(sample_label, ".rds"))
    saveRDS(seurat_sample, out_rds)
    log_msg("  Saved: ", out_rds)
  }
}

writeLines(log_lines, LOG_FILE)
message("\nDone. Per-sample CS-CORE input objects written to ", OUT_DIR,
        " | log: ", LOG_FILE)
