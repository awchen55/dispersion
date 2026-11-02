#!/usr/bin/env Rscript

stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(Seurat)
})

# ---------------------------------------------------------------------------
# NAMED CONSTANTS
# ---------------------------------------------------------------------------
RAW_COUNTS_SEURAT_RDS <- "cellranger_cluster-mode_trial/analysis/datasets/100k_full_filtered_norm_seurat.rds"
OUT_DIR <- "figure_reproduction/outputs/coexpression_plasticity_lowfidelity"
BARCODES_TXT <- file.path(OUT_DIR, "epithelial_excl_cluster24_barcodes.txt")
OUT_RDS <- file.path(OUT_DIR, "epithelial_excl_cluster24_seurat_subset.rds")

CELLTYPE_COL <- "new.predicted.id"
TARGET_CELL_TYPE <- "Epithelial_cells"
# Same expectation as script 00a -- both must agree since they describe the
# same underlying cell set from two different files.
EXPECTED_N_CELLS <- 26914

log_msg <- function(...) message("[00_subset_epithelial_excl_cluster24_seurat] ", ...)

stopifnot(
  "raw-counts Seurat RDS not found" = file.exists(RAW_COUNTS_SEURAT_RDS),
  "barcode list not found -- run script 00a first" = file.exists(BARCODES_TXT)
)

barcodes <- readLines(BARCODES_TXT)
log_msg("Read ", length(barcodes), " barcodes from ", BARCODES_TXT)
stopifnot(
  "barcode count mismatch vs. epithelial_subcluster24_excl_sensitivity's known cell count" =
    length(barcodes) == EXPECTED_N_CELLS,
  "duplicate barcodes in barcode list" = !anyDuplicated(barcodes)
)

log_msg("Loading full Seurat object (heavy): ", RAW_COUNTS_SEURAT_RDS)
obj <- readRDS(RAW_COUNTS_SEURAT_RDS)
stopifnot("expected a Seurat object" = inherits(obj, "Seurat"))
log_msg("Total cells: ", ncol(obj), " | total genes: ", nrow(obj))

md <- obj[[]]
stopifnot("cell type column missing" = CELLTYPE_COL %in% colnames(md))

# --- barcode-intersection subset --------------------------------------------
present <- barcodes %in% colnames(obj)
n_present <- sum(present)
n_missing <- sum(!present)
log_msg(n_present, " / ", length(barcodes), " target barcodes found in the Seurat object (",
        n_missing, " missing).")
if (n_missing > 0) {
  stop("Barcode format mismatch or cell-set drift between the h5ad and this Seurat object -- ",
       n_missing, " barcodes from ", BARCODES_TXT, " are not present in ", RAW_COUNTS_SEURAT_RDS,
       ". Do not silently proceed with a partial cell set -- investigate the barcode formats ",
       "in both files before rerunning.")
}

obj_sub <- subset(obj, cells = barcodes)
log_msg("Subset to ", ncol(obj_sub), " cells x ", nrow(obj_sub), " genes.")
stopifnot(
  "subset cell count does not match the target barcode list" = ncol(obj_sub) == length(barcodes)
)

# Confirm every retained cell really is Epithelial_cells in this object's own
# metadata (defensive cross-check, not just trusting the h5ad's own labeling
# carried over correctly).
ct_vals <- unique(as.character(md[barcodes, CELLTYPE_COL]))
if (!identical(ct_vals, TARGET_CELL_TYPE)) {
  stop("Subset cells are not uniformly '", TARGET_CELL_TYPE, "' in this Seurat object's own '",
       CELLTYPE_COL, "' column -- found: ", paste(ct_vals, collapse = ", "),
       ". Investigate before proceeding (barcode misalignment risk).")
}

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(obj_sub, OUT_RDS)
log_msg("Wrote: ", OUT_RDS)
log_msg("DONE.")
