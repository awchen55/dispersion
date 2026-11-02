THRESHOLDS <- c(1000, 800, 600, 400, 200)

# Source repo (cHDC_data) — this repo intentionally reuses its already-computed
# downsampled Memento matrices and DD fit objects rather than recomputing the
# Memento step (per 2026-09-01 review discussion: the downsampling/Memento step
# itself was not flagged as buggy; only the model formula used for the
# cross-cohort comparison reference point, and the reproducibility of the
# driver scripts, were). This repo has its own independent git history.
SOURCE_REPO <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"

BRENDAN_DD_FIT_DIR   <- file.path(SOURCE_REPO, "analysis/qtl_analyses/outputs/dd_sweep/brendan")
ERIK_B2_DD_FIT_DIR   <- file.path(SOURCE_REPO, "analysis/qtl_analyses/outputs/dd_sweep_downsampled/b2")
ERIK_B6_DD_FIT_DIR   <- file.path(SOURCE_REPO, "analysis/qtl_analyses/outputs/dd_sweep_downsampled/b6")
BRENDAN_FULL_MEMENTO_DIR <- file.path(SOURCE_REPO, "analysis/qtl_analyses/outputs/brendan_cell_count_sweep/cell_count_full")

BRENDAN_CELL_TYPES <- c(
  "Activated_fibroblasts", "Cardiac_fibroblasts", "Cardiomyocytes",
  "Epithelial_cells", "Foregut_endodermal_cells", "Mesenchymal_cells"
)

# The 4 contrasts shared between Brendan's data and Erik's B2/B6 batches,
# and which batch each is computed in for Erik (matches
# analysis/qtl_analyses/scripts/plot_dd_sweep_downsampled.R's
# shared_contrasts list in the source repo).
SHARED_CONTRASTS <- data.frame(
  contrast = c(
    "compare_Cardiac_fibroblasts_Activated_fibroblasts",
    "compare_Mesenchymal_cells_Activated_fibroblasts",
    "compare_Mesenchymal_cells_Cardiac_fibroblasts",
    "compare_Cardiomyocytes_Cardiac_fibroblasts"
  ),
  erik_batch = c("b2", "b2", "b2", "b6"),
  stringsAsFactors = FALSE
)

erik_dd_fit_dir <- function(batch) {
  if (batch == "b2") ERIK_B2_DD_FIT_DIR else ERIK_B6_DD_FIT_DIR
}
