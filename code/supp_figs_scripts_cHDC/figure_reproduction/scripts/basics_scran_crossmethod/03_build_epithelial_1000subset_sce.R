stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))
suppressPackageStartupMessages(library(SingleCellExperiment))

REPO_DIR <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
SCE_DIR <- file.path(REPO_DIR, "figure_reproduction/outputs/basics_scran_crossmethod/sce")
SCE_IN <- file.path(SCE_DIR, "Epithelial_cells_sce_raw.rds")
BARCODES_CSV <- file.path(REPO_DIR, "figure_reproduction/outputs/basics_scran_crossmethod",
                          "Epithelial_cells_1000subset_barcodes.csv")
SCE_OUT <- file.path(SCE_DIR, "Epithelial_cells_1000subset_sce_raw.rds")

MIN_BATCHES <- 2L
CELLS_PER_INDIVIDUAL <- 1000L
EXPECTED_CELLS <- 3L * CELLS_PER_INDIVIDUAL

stopifnot(
  "input SCE not found -- run 01_build_epithelial_sce.R first" = file.exists(SCE_IN),
  "barcode CSV not found -- run 02_new_epithelial_1000subset_memento.py first" = file.exists(BARCODES_CSV)
)

sce_full <- readRDS(SCE_IN)
message("Full Epithelial_cells SCE: ", nrow(sce_full), " genes x ", ncol(sce_full), " cells")

barcodes <- read.csv(BARCODES_CSV, stringsAsFactors = FALSE)$barcode
stopifnot("barcode CSV is empty" = length(barcodes) > 0)

missing_barcodes <- setdiff(barcodes, colnames(sce_full))
stopifnot("some target barcodes not found in the full SCE -- barcode format mismatch?" =
            length(missing_barcodes) == 0)

sce_sub <- sce_full[, barcodes]
stopifnot(
  "unexpected cell count after subsetting" = ncol(sce_sub) == EXPECTED_CELLS,
  "gene count changed during subsetting" = nrow(sce_sub) == nrow(sce_full)
)

batch_tab <- table(colData(sce_sub)$BatchInfo)
message("BatchInfo composition: ", paste(names(batch_tab), "=", batch_tab, collapse = ", "))
stopifnot(
  "no-spike-in mode requires >= MIN_BATCHES batches" = length(batch_tab) >= MIN_BATCHES,
  "expected exactly CELLS_PER_INDIVIDUAL cells per individual" = all(batch_tab == CELLS_PER_INDIVIDUAL)
)

saveRDS(sce_sub, SCE_OUT)
message("Saved: ", SCE_OUT, " (", nrow(sce_sub), " genes x ", ncol(sce_sub), " cells)")
