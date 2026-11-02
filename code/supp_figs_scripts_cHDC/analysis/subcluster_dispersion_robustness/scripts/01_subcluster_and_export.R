#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
})

# ---------------------------------------------------------------------------
# COMPATIBILITY PATCH: SeuratDisk's WriteH5Group(Assay-method) (used inside
# SaveH5Seurat) calls GetAssayData(object = x, slot = i), which is defunct on
# the installed SeuratObject (>= 5.0). This is an unfixed upstream bug (open
# since 2026-02, mojaveazure/seurat-disk#199, no official patch as of
# 2026-07-30) -- SeuratDisk is not actively maintained against SeuratObject's
# slot->layer migration. Rather than downgrade SeuratObject repo-wide (would
# break every other script's layer= calls), override just this one S4 method
# for the duration of THIS script's session: an exact copy of the original
# method body with slot= -> layer= in the one affected GetAssayData call.
# environment() is set to the SeuratDisk namespace so the copy can still see
# SeuratDisk's unexported helpers (GuessDType, IsMatrixEmpty, etc.) by lexical
# scoping. Verified end-to-end (SaveH5Seurat + Convert to h5ad) against a
# small dummy Seurat object before use here -- see .living/learnings.md,
# 2026-07-30.
patched_write_h5group_assay <- function(x, name, hgroup, verbose = TRUE) {
  xgroup <- hgroup$create_group(name = name)
  for (i in c("counts", "data", "scale.data")) {
    dat <- GetAssayData(object = x, layer = i)
    if (!IsMatrixEmpty(x = dat)) {
      if (verbose) message("Adding ", i, " for ", name)
      WriteH5Group(x = dat, name = i, hgroup = xgroup, verbose = verbose)
    }
    if (i == "scale.data") {
      WriteH5Group(x = rownames(x = dat), name = "scaled.features", hgroup = xgroup, verbose = verbose)
    }
  }
  WriteH5Group(x = rownames(x = x), name = "features", hgroup = xgroup, verbose = verbose)
  xgroup$create_attr(attr_name = "key", robj = Key(object = x), dtype = GuessDType(x = Key(object = x)))
  if (length(x = VariableFeatures(object = x))) {
    if (verbose) message("Adding variable features for ", name)
    WriteH5Group(x = VariableFeatures(object = x), name = "variable.features", hgroup = xgroup, verbose = verbose)
  } else if (verbose) {
    message("No variable features found for ", name)
  }
  if (ncol(x = x[[]])) {
    if (verbose) message("Adding feature-level metadata for ", name)
    WriteH5Group(x = x[[]], name = "meta.features", hgroup = xgroup, verbose = verbose)
  } else if (verbose) {
    message("No feature-level metadata found for ", name)
  }
  WriteH5Group(x = Misc(object = x), name = "misc", hgroup = xgroup, verbose = verbose)
  if (class(x = x)[1] != "Assay") {
    extclass <- attr(x = SeuratObject::S4ToList(object = x), which = "classDef")
    xgroup$create_attr(attr_name = "s4class", robj = extclass, dtype = GuessDType(x = extclass))
    slots.extended <- setdiff(x = slotNames(x = x), y = slotNames(x = methods::getClassDef(Class = "Assay")))
    for (slot in slots.extended) {
      if (verbose) message("Writing out ", slot, " for ", name)
      WriteH5Group(x = slot(object = x, name = slot), name = slot, hgroup = xgroup, verbose = verbose)
    }
  }
  return(invisible(x = NULL))
}
environment(patched_write_h5group_assay) <- asNamespace("SeuratDisk")
setMethod("WriteH5Group", "Assay", patched_write_h5group_assay)

# ---- Named constants (lowercase per R convention) ---------------------------
repo_dir     <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
trial_dir    <- file.path(repo_dir, "cellranger_cluster-mode_trial")
analysis_dir <- file.path(repo_dir, "analysis/subcluster_dispersion_robustness")
out_dir      <- file.path(analysis_dir, "outputs")

seurat_rds <- file.path(trial_dir, "analysis/datasets/100k_full_filtered_norm_seurat.rds")

cell_type_col  <- "new.predicted.id"
individual_col <- "vireo.individual"

cell_types <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)
# ANALYSIS_OK[fixed-design]: the 3 pooled cHDC iPSC lines (the entire study
# cohort), not an exclusion filter; presence is asserted at the setequal() check
# below.
individuals <- c("NA19099", "NA19093", "NA18522")

# TODO(assumption): resolutions per REVISION_PLAN sec 3.1 -- confirm
resolutions <- c(0.5, 1.0, 2.0)

# Deterministic clustering.
seed <- 1L

# Raw counts must be integer within floating-point tolerance before export to
# Memento; a larger deviation means the counts slot holds normalized data.
max_int_deviation_tol <- 1e-8

# SeuratDisk maps the RNA v3 "data" slot to adata.X. We place raw counts there.
h5seurat_path <- file.path(out_dir, "subcluster_raw_counts.h5Seurat")
h5ad_path     <- file.path(out_dir, "subcluster_raw_counts.h5ad")
counts_log    <- file.path(out_dir, "cell_counts_per_individual_subcluster.csv")

# Resolutions map to Seurat metadata columns "SCT_snn_res.<res>". Sanitize the
# resolution into a stable suffix used for the exported sub-cluster column.
# ANALYSIS_OK[within-file-helper]: res_suffix/subcluster_col are used within this
# script (subcluster_col in the FindClusters loop and the export block below;
# res_suffix inside subcluster_col); R026's cross-file rule cannot see intra-file
# callers.
res_suffix <- function(res) gsub("\\.", "p", format(res, trim = TRUE))
subcluster_col <- function(res) paste0("subcluster_res_", res_suffix(res))

# ---- Helpers ----------------------------------------------------------------
# Logging logic lives in scripts/utils.R (shared with 03) to avoid duplicate
# definitions drifting apart (scilintr R025).
source(file.path(repo_dir, "analysis/subcluster_dispersion_robustness/scripts/utils.R"))
log_msg <- make_logger("01_subcluster")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Load + subset ----------------------------------------------------------
log_msg("Loading Seurat object: %s", seurat_rds)
stopifnot(file.exists(seurat_rds))
obj <- readRDS(seurat_rds)
log_msg("Loaded object with %d cells x %d features.", ncol(obj), nrow(obj))

# Fail loudly if expected metadata columns are missing.
required_cols <- c(cell_type_col, individual_col)
missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required metadata columns: %s",
               paste(missing_cols, collapse = ", ")))
}

# Confirm the annotated cell types / individuals actually exist.
present_cts  <- intersect(cell_types, unique(obj@meta.data[[cell_type_col]]))
present_inds <- intersect(individuals, unique(obj@meta.data[[individual_col]]))
if (!setequal(present_cts, cell_types)) {
  stop(sprintf("Expected cell types not all present. Missing: %s",
               paste(setdiff(cell_types, present_cts), collapse = ", ")))
}
if (!setequal(present_inds, individuals)) {
  stop(sprintf("Expected individuals not all present. Missing: %s",
               paste(setdiff(individuals, present_inds), collapse = ", ")))
}

keep <- obj@meta.data[[cell_type_col]] %in% cell_types &
        obj@meta.data[[individual_col]] %in% individuals
log_msg("Subsetting to 6 cell types x 3 individuals: %d of %d cells retained.",
        sum(keep), ncol(obj))
# ANALYSIS_OK[subset-to-design]: subset to the 6 manuscript cell types x 3
# individuals; retained-cell count logged just above and asserted just below.
obj <- obj[, keep]
stopifnot(ncol(obj) == sum(keep))

# ---- Sub-cluster at each resolution -----------------------------------------
# FindClusters on the SCT graph produces a global fine clustering; intersecting
# it with the annotated cell type gives sub-clusters within each cell type via
# the combined "{cell_type}::{cluster}" label below.
#
# TODO(assumption): We cluster the full 6-cell-type subset on the existing
# SCT_snn graph rather than recomputing neighbors within each cell type. This
# matches variability_and_cluster_resolution.Rmd. Confirm this is the intended
# sub-clustering strategy (vs. per-cell-type FindNeighbors+FindClusters).
DefaultAssay(obj) <- "SCT"
graph_ok <- "SCT_snn" %in% names(obj@graphs)
if (!graph_ok) {
  stop(paste0("SCT_snn graph not found on the subset object. FindClusters ",
              "needs an existing SNN graph. Recompute with FindNeighbors on ",
              "the SCT assay before sub-clustering (not done here to avoid ",
              "silently changing the manuscript neighbor graph)."))
}

for (res in resolutions) {
  log_msg("FindClusters at resolution %.3f (seed=%d) ...", res, seed)
  obj <- FindClusters(obj, resolution = res, graph.name = "SCT_snn",
                      random.seed = seed, verbose = TRUE)
  # FindClusters writes to obj$seurat_clusters and obj$SCT_snn_res.<res>.
  res_col <- paste0("SCT_snn_res.", format(res, trim = TRUE))
  if (!res_col %in% colnames(obj@meta.data)) {
    # Fall back to seurat_clusters if the naming differs across Seurat versions.
    res_col_alt <- "seurat_clusters"
    if (!res_col_alt %in% colnames(obj@meta.data)) {
      stop(sprintf("Could not locate cluster assignments for resolution %s.", res))
    }
    res_col <- res_col_alt
  }
  cluster_ids <- as.character(obj@meta.data[[res_col]])
  # Combined sub-cluster label: {cell_type}::{cluster}
  combined <- paste(obj@meta.data[[cell_type_col]], cluster_ids, sep = "::")
  obj@meta.data[[subcluster_col(res)]] <- combined
  log_msg("  resolution %.3f -> %d global clusters, %d combined sub-clusters.",
          res, length(unique(cluster_ids)), length(unique(combined)))
}

# ---- Log cell counts per (individual x sub-cluster) per resolution ----------
count_rows <- list()
for (res in resolutions) {
  sc_col <- subcluster_col(res)
  tab <- as.data.frame(table(
    individual = obj@meta.data[[individual_col]],
    subcluster = obj@meta.data[[sc_col]]
  ))
  # ANALYSIS_OK[drop-empty-crosstab-cells]: removes zero-count (individual x
  # sub-cluster) cells of the contingency table, not analysis-population rows;
  # the retained rows are written to counts_log and their count logged below.
  tab <- tab[tab$Freq > 0, ]
  tab$resolution <- format(res, trim = TRUE)
  count_rows[[as.character(res)]] <- tab
}
counts_df <- do.call(rbind, count_rows)
# ANALYSIS_OK[column-reorder]: selects/reorders columns only (no row drop); row
# count is unchanged and logged at the write.csv() below.
counts_df <- counts_df[, c("resolution", "individual", "subcluster", "Freq")]
colnames(counts_df)[colnames(counts_df) == "Freq"] <- "n_cells"
write.csv(counts_df, counts_log, row.names = FALSE)
log_msg("Wrote per-(individual x sub-cluster) counts: %s (%d rows).",
        counts_log, nrow(counts_df))
log_msg("Sub-cluster group size summary (cells per individual x sub-cluster): min=%d median=%.0f max=%d.",
        min(counts_df$n_cells), median(counts_df$n_cells), max(counts_df$n_cells))

# ---- Export RAW COUNTS to h5ad ----------------------------------------------
# Raw counts live in the RNA assay counts slot.
DefaultAssay(obj) <- "RNA"
raw_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
log_msg("RNA counts matrix: %d genes x %d cells.", nrow(raw_counts), ncol(raw_counts))

# Assert integrality of raw counts (fail loudly, do NOT coerce silently).
nz <- raw_counts@x
if (length(nz) == 0) stop("RNA counts slot is empty -- nothing to export.")
if (any(nz < 0)) stop("RNA counts contain negative values -- not raw counts.")
max_int_err <- max(abs(nz - round(nz)))
if (max_int_err > max_int_deviation_tol) {
  stop(sprintf(paste0("RNA counts are NOT integer-valued (max deviation %.3g). ",
                      "The counts slot appears to hold normalized data. Aborting ",
                      "to avoid feeding non-raw counts to Memento."), max_int_err))
}
log_msg("Raw-count integrality assertion passed (max deviation from integer = %.3g).",
        max_int_err)

# Build a minimal export object whose data slot == raw counts so that the
# SeuratDisk-produced adata.X is the raw counts (Memento reads adata.X).
export_meta <- obj@meta.data[, c(individual_col, cell_type_col,
                                 vapply(resolutions, subcluster_col, character(1))),
                             drop = FALSE]
# CreateAssayObject() builds a v3 Assay directly (not the default Assay5),
# with @data == @counts by default -- this both satisfies "data slot == raw
# counts" and sidesteps as(Assay5, "Assay") coercion, which errors on the
# installed SeuratObject 5.4.0 ("subscript out of bounds" from an empty
# scale.data layer during the coercion) -- see .living/learnings.md, 2026-07-30.
export_obj <- CreateSeuratObject(counts = CreateAssayObject(counts = raw_counts),
                                 meta.data = export_meta)
DefaultAssay(export_obj) <- "RNA"

# Sanity: confirm data slot integrality survived the round-trip into export_obj.
data_nz <- GetAssayData(export_obj, assay = "RNA", layer = "data")@x
if (max(abs(data_nz - round(data_nz))) > max_int_deviation_tol) {
  stop("Export object data slot is not integer-valued after assignment. Aborting.")
}

if (file.exists(h5seurat_path)) file.remove(h5seurat_path)
if (file.exists(h5ad_path))     file.remove(h5ad_path)
log_msg("Writing h5Seurat -> h5ad ...")
SaveH5Seurat(export_obj, filename = h5seurat_path, overwrite = TRUE)
Convert(h5seurat_path, dest = "h5ad", overwrite = TRUE)
log_msg("Wrote h5ad: %s", h5ad_path)
log_msg("DONE. Next: scripts/02_memento_per_subcluster.py")
