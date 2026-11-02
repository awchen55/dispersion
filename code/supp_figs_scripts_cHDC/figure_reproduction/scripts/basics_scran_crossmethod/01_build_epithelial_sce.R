stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
})

REPO_DIR <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
RAW_COUNTS_SEURAT_RDS <- file.path(
  REPO_DIR, "cellranger_cluster-mode_trial/analysis/datasets/100k_full_filtered_norm_seurat.rds")
BARCODES_CSV <- file.path(
  REPO_DIR, "figure_reproduction/outputs/basics_scran_crossmethod",
  "epithelial_excl_cluster24_barcodes.csv")
MEMENTO_MEAN_CSV <- file.path(REPO_DIR, "figure_reproduction/data/mean_expression_excl_cluster24.csv")
OUT_DIR <- file.path(REPO_DIR, "figure_reproduction/outputs/basics_scran_crossmethod/sce")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CELLTYPE_COL <- "new.predicted.id"
BATCH_COL    <- "vireo.individual"
CELL_TYPE    <- "Epithelial_cells"
INDIVIDUALS  <- c("NA19099", "NA19093", "NA18522")
MIN_BATCHES  <- 2L

RIBO_PATTERN <- "^RP"
MITO_PATTERN <- "^MT-"
MEAN_EXCL_THRESHOLD <- -6

# Minimum fraction of the excl-cluster24 h5ad's barcodes that must be found in
# this Seurat object -- guards against a barcode-format mismatch between the
# two source files (h5ad derived from subcluster_dispersion_robustness's
# export, Seurat RDS from a separate pipeline stage) going undetected. The
# h5ad's own doc records 26,914 retained Epithelial cells; a near-complete
# match is expected, not a partial one.
MIN_BARCODE_MATCH_FRAC <- 0.99

stopifnot(
  "raw-counts Seurat RDS not found" = file.exists(RAW_COUNTS_SEURAT_RDS),
  "barcode CSV not found -- run 00_export_epithelial_excl_barcodes.py first" = file.exists(BARCODES_CSV),
  "Epithelial excl-cluster24 Memento mean CSV not found" = file.exists(MEMENTO_MEAN_CSV)
)

message("Loading Seurat object (heavy): ", RAW_COUNTS_SEURAT_RDS)
obj <- readRDS(RAW_COUNTS_SEURAT_RDS)
stopifnot("expected a Seurat object" = inherits(obj, "Seurat"))
md <- obj[[]]
stopifnot(
  CELLTYPE_COL %in% colnames(md),
  BATCH_COL %in% colnames(md)
)

counts_all <- GetAssayData(obj, assay = "RNA", layer = "counts")
set.seed(1L)
samp_idx <- sample(length(counts_all@x), size = min(2e6L, length(counts_all@x)))
vals <- counts_all@x[samp_idx]
stopifnot(
  "RNA counts slot contains non-integer values" = all(vals == round(vals)),
  "RNA counts slot contains negative values" = all(vals >= 0)
)

excl_barcodes <- read.csv(BARCODES_CSV, stringsAsFactors = FALSE)$barcode
stopifnot("barcode CSV is empty" = length(excl_barcodes) > 0)

cells_ct_all <- rownames(md)[as.character(md[[CELLTYPE_COL]]) == CELL_TYPE]
stopifnot("no Epithelial_cells cells found in Seurat object" = length(cells_ct_all) > 0)

n_match <- length(intersect(excl_barcodes, cells_ct_all))
match_frac <- n_match / length(excl_barcodes)
message("Barcode overlap: ", n_match, " / ", length(excl_barcodes),
        " excl-cluster24 barcodes found among ", length(cells_ct_all),
        " Epithelial_cells cells in the Seurat object (", round(match_frac, 4), ")")
stopifnot(
  "Barcode overlap between the excl-cluster24 h5ad and this Seurat object is too low -- likely a barcode-format mismatch; do NOT proceed" =
    match_frac >= MIN_BARCODE_MATCH_FRAC
)

cells_ct <- intersect(excl_barcodes, cells_ct_all)
batch <- as.character(md[cells_ct, BATCH_COL])
keep_cell <- batch %in% INDIVIDUALS
if (any(!keep_cell)) {
  message("Dropping ", sum(!keep_cell), " cell(s) not in target individuals")
}
cells_ct <- cells_ct[keep_cell]
batch <- batch[keep_cell]
batch_tab <- table(batch)
for (ind in INDIVIDUALS) {
  message("  batch ", ind, ": ", ifelse(ind %in% names(batch_tab), batch_tab[[ind]], 0L), " cells")
}
n_batches <- sum(batch_tab > 0)
stopifnot("BASiCS no-spike-in mode requires >= MIN_BATCHES batches" = n_batches >= MIN_BATCHES)

cts_mat <- counts_all[, cells_ct, drop = FALSE]

ribo <- grepl(RIBO_PATTERN, rownames(cts_mat))
message("Dropping ", sum(ribo), " ribosomal genes")
cts_mat <- cts_mat[!ribo, , drop = FALSE]

mito <- grepl(MITO_PATTERN, rownames(cts_mat))
message("Dropping ", sum(mito), " mitochondrial genes")
cts_mat <- cts_mat[!mito, , drop = FALSE]

memento_mean <- read.csv(MEMENTO_MEAN_CSV, row.names = 1, check.names = FALSE)
ct_mean_cols <- paste0(INDIVIDUALS, "_", CELL_TYPE)
present_cols <- ct_mean_cols %in% colnames(memento_mean)
stopifnot("no Memento mean columns found for Epithelial_cells" = any(present_cols))
ct_mean <- memento_mean[, ct_mean_cols[present_cols], drop = FALSE]
mean_excluded <- rownames(ct_mean)[apply(ct_mean, 1, function(x) any(x > MEAN_EXCL_THRESHOLD))]

in_memento <- rownames(cts_mat) %in% rownames(memento_mean)
above_thresh <- rownames(cts_mat) %in% mean_excluded
message("Genes absent from Memento mean matrix: ", sum(!in_memento))
message("Genes with Memento mean > ", MEAN_EXCL_THRESHOLD, " in >= 1 individual: ", sum(in_memento & above_thresh))

keep_gene <- in_memento & !above_thresh
cts_mat <- cts_mat[keep_gene, , drop = FALSE]
stopifnot("no genes survive filtering" = nrow(cts_mat) > 0)
message("Genes retained: ", nrow(cts_mat))

cts_dense <- as.matrix(cts_mat)
storage.mode(cts_dense) <- "integer"
stopifnot(!anyNA(cts_dense))

sce <- SingleCellExperiment(
  assays = list(counts = cts_dense),
  colData = DataFrame(BatchInfo = batch, row.names = colnames(cts_dense))
)
message("SCE built: ", nrow(sce), " genes x ", ncol(sce), " cells")

out_rds <- file.path(OUT_DIR, paste0(CELL_TYPE, "_sce_raw.rds"))
saveRDS(sce, out_rds)
message("Saved: ", out_rds)
