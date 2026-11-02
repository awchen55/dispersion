stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(scran)
  library(scuttle)
  library(SingleCellExperiment)
})

# ---------------------------------------------------------------------------
# NAMED CONSTANTS
# ---------------------------------------------------------------------------
CLUSTER_SEED <- 42L  # named seed for quickCluster's k-means step (reproducibility)

CELL_TYPES <- c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                "Epithelial_cells", "Activated_fibroblasts",
                "Foregut_endodermal_cells")

# PATCHED (figure_reproduction/PORT_CHANGELOG_crossmethod.md): SCE_DIR now
# points at this fork's mirror sce/ directory (5 cell types symlinked from
# the original analysis, Epithelial_cells rebuilt excl-cluster24 -- see
# setup_mirror_dirs.sh); OUT_DIR is this fork's own output tree. Model code
# below is otherwise byte-identical to
# analysis/scran_dispersion_reanalysis/scripts/01_run_scran_modelgenevar.R.
SCE_DIR <- "figure_reproduction/outputs/basics_scran_crossmethod/sce"
OUT_DIR <- "figure_reproduction/outputs/basics_scran_crossmethod/scran/bio_variance"
LOG_FILE <- file.path(OUT_DIR, "run_scran_log.txt")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

log_lines <- character(0)
# ANALYSIS_OK[within-file-helper]: log_msg is called throughout this script
# (cross-file "dead code" is a false positive -- it has no callers in OTHER
# files by design; each prep script keeps its own logger, matching
# basics_dispersion_reanalysis/scripts/01_prepare_basics_inputs.R's convention).
log_msg <- function(...) {
  line <- paste0(...)
  message(line)
  log_lines[[length(log_lines) + 1L]] <<- line
}

for (ct in CELL_TYPES) {
  log_msg("\n=== ", ct, " ===")
  sce_rds <- file.path(SCE_DIR, paste0(ct, "_sce_raw.rds"))
  stopifnot("prerequisite SCE not found -- run basics_dispersion_reanalysis/scripts/01_prepare_basics_inputs.R first" =
              file.exists(sce_rds))

  sce <- readRDS(sce_rds)
  stopifnot(
    "expected a SingleCellExperiment"    = inherits(sce, "SingleCellExperiment"),
    "counts assay missing"               = "counts" %in% assayNames(sce),
    "BatchInfo missing from colData"     = "BatchInfo" %in% colnames(colData(sce))
  )
  batch <- factor(colData(sce)$BatchInfo)
  n_batches <- nlevels(batch)
  log_msg("  Loaded: ", nrow(sce), " genes x ", ncol(sce), " cells, ",
          n_batches, " individuals (", paste(levels(batch), collapse = ", "), ")")
  stopifnot("need >= 2 individuals to block modelGeneVar/quickCluster" = n_batches >= 2L)

  # ---------------------------------------------------------------------
  # NORMALIZE: pooling-based size factors (scran), blocked on individual.
  # ---------------------------------------------------------------------
  set.seed(CLUSTER_SEED)
  clusters <- quickCluster(sce, block = batch)
  log_msg("  quickCluster: ", length(unique(clusters)), " clusters (blocked on individual)")
  sce <- computeSumFactors(sce, clusters = clusters)
  stopifnot(
    "computeSumFactors produced non-finite size factors" = all(is.finite(sizeFactors(sce))),
    "computeSumFactors produced non-positive size factors" = all(sizeFactors(sce) > 0)
  )
  log_msg("  Size factors: min=", round(min(sizeFactors(sce)), 3),
          " median=", round(median(sizeFactors(sce)), 3),
          " max=", round(max(sizeFactors(sce)), 3))
  sce <- logNormCounts(sce)
  stopifnot("logcounts assay missing after logNormCounts" = "logcounts" %in% assayNames(sce))

  # ---------------------------------------------------------------------
  # MODEL GENE VARIANCE: mean-variance trend fit per individual, combined.
  # ---------------------------------------------------------------------
  var_out <- modelGeneVar(sce, block = batch)
  var_df <- as.data.frame(var_out)
  var_df$gene <- rownames(var_out)
  var_df <- var_df[, c("gene", setdiff(colnames(var_df), "gene"))]
  stopifnot(
    "modelGeneVar output has no rows"      = nrow(var_df) > 0,
    "modelGeneVar output missing 'bio' column" = "bio" %in% colnames(var_df),
    "gene count changed unexpectedly"      = nrow(var_df) == nrow(sce)
  )
  n_finite_bio <- sum(is.finite(var_df$bio))
  log_msg("  modelGeneVar: ", nrow(var_df), " genes | finite bio: ", n_finite_bio,
          " | non-finite bio (excluded downstream): ", nrow(var_df) - n_finite_bio)

  out_csv <- file.path(OUT_DIR, paste0(ct, "_bio_variance.csv"))
  write.csv(var_df, out_csv, row.names = FALSE)
  log_msg("  Saved: ", out_csv)
}

writeLines(log_lines, LOG_FILE)
message("\nDone. Per-cell-type scran bio-variance tables written to ", OUT_DIR, " | log: ", LOG_FILE)
