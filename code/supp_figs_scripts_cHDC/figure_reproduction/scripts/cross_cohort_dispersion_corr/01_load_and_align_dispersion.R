stopifnot(file.exists("CLAUDE.md"))  # must run from project root

# Shared filter_cell_types_by_erik_coverage()/read_active_cell_types() live
# in utils.R (single canonical definition; see R025).
# PATCHED for figure_reproduction (see PORT_CHANGELOG_crosscohort.md): points
# at this folder's own copy of utils.R, not the source analysis's.
source(file.path("figure_reproduction", "scripts", "cross_cohort_dispersion_corr", "utils.R"))

# ---------------------------------------------------------------------------
# CONSTANTS
# ---------------------------------------------------------------------------
# Fixed candidate list (the 6 cell types this analysis covers). The set
# actually used downstream (CELL_TYPES, below) is a per-threshold subset of
# this — see filter_cell_types_by_erik_coverage() in utils.R.
ALL_CELL_TYPES <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)

# Brendan discovery cohort: Memento residual variance matrix.
# PATCHED for figure_reproduction (see PORT_CHANGELOG_crosscohort.md): reads
# the Epithelial_cells sub-cluster-24-excluded matrix (only Epithelial's 3
# columns differ from the original; identical shape/column order otherwise)
# instead of the published canonical matrix.
BRENDAN_DISP_CSV <- file.path("figure_reproduction", "data", "res_var_excl_cluster24.csv")
# Brendan discovery cohort: Memento mean-expression matrix, same shape/column
# order as BRENDAN_DISP_CSV — feeds the -6 mean-threshold exclusion
# (exclude_top_expressed() in utils.R; see REVISION_ANALYSIS_PLAN.md item 1).
BRENDAN_MEAN_CSV <- file.path("figure_reproduction", "data", "mean_expression_excl_cluster24.csv")
# ANALYSIS_OK[fixed-design]: the 3 pooled cHDC iPSC lines (NA19099, NA19093,
# NA18522) are the discovery cohort's fixed experimental design, not an
# exclusion filter; used to build the expected {IND}_{CT} column set that is
# asserted present (setdiff/stopifnot) at load below.
BRENDAN_INDS <- c("NA19099", "NA19093", "NA18522")  # the 3 discovery individuals

# Erik cohort cell-count-threshold sweep. Threshold is env-overridable so the
# sensitivity re-run (cell_count_50 vs cell_count_150; see folder Open
# questions) is a one-line change plus a re-run of 01 then 02 — outputs are
# namespaced by threshold below.
ERIK_CELL_COUNT_THRESHOLD <- as.integer(
  Sys.getenv("CCDC_CELL_COUNT_THRESHOLD", unset = "150")
)
# ANALYSIS_OK[pinned-input]: analysis/qtl_analyses/outputs/cell_count_threshold_sweep
# is Erik's cell-count-threshold sweep output; the specific release is pinned by
# ERIK_CELL_COUNT_THRESHOLD (cell_count_150 by default) in the path built below,
# and each per-CT file's existence is asserted at read time (stopifnot file.exists).
ERIK_SWEEP_BASE <- file.path(
  "analysis", "qtl_analyses", "outputs", "cell_count_threshold_sweep"
)
# TODO(assumption): The task brief stated that new_dispersion_by_ind_*.csv
# exists for only 4 cell types. That is INCORRECT as of 2026-07-20 — verified
# that BOTH dispersion_by_ind_*.csv and new_dispersion_by_ind_*.csv exist for
# all 6 cell types at cell_count_150 AND cell_count_50, and hold identical
# dispersion values (new_* is the same data re-saved with gene as the row name
# and no leading integer index). new_dispersion_by_ind_* is chosen here because
# (a) it covers all 6 cell types, (b) its rownames=gene layout matches Brendan's
# matrix, and (c) existing repo scripts (analysis/qtl_analyses/scripts/
# compute_and_plot_correlations.R, analysis/dqtl_mapping/scripts/
# 01_prepare_dispersion_inputs.R) already consume this exact file set. To switch
# to the dispersion_by_ind_* layout the loader below must change (that format
# has a leading integer index column and a separate 'gene' column).
ERIK_DISP_PREFIX <- "new_dispersion_by_ind"
# Erik's per-CT mean-expression matrix, same shape/sample-columns as
# ERIK_DISP_PREFIX's file for that cell type — feeds the -6 mean-threshold
# exclusion (exclude_top_expressed() in utils.R).
ERIK_MEAN_PREFIX <- "new_mean_by_ind"

# Below this many shared genes for a cell type, warn (gene-ID overlap likely
# too low to trust the correlation). Soft check — analysis proceeds regardless.
MIN_EXPECTED_SHARED_GENES <- 1000L

# Drop ribosomal genes first (repo convention). NOTE: "^RP" is a broad prefix —
# it removes RPL*/RPS* ribosomal genes but also incidentally a handful of
# non-ribosomal symbols (e.g. RPTOR, RPGR). This is the documented repo
# convention and is applied identically to both cohorts.
RIBO_PATTERN <- "^RP"

# PATCHED for figure_reproduction: outputs land under this folder, not the
# source analysis's own outputs/.
OUT_ROOT <- file.path("figure_reproduction", "outputs", "cross_cohort_dispersion_corr")
OUT_DIR  <- file.path(OUT_ROOT, paste0("cc", ERIK_CELL_COUNT_THRESHOLD))
ALIGN_DIR <- file.path(OUT_DIR, "aligned")

dir.create(ALIGN_DIR, recursive = TRUE, showWarnings = FALSE)

message("=== 01_load_and_align_dispersion.R ===")
message("Erik cell-count threshold: ", ERIK_CELL_COUNT_THRESHOLD)
message("Output dir: ", OUT_DIR)

# Cell types actually used below: ALL_CELL_TYPES filtered down to those with
# >= 3 individuals in Erik's data at this threshold (default in
# filter_cell_types_by_erik_coverage). Persisted to
# OUT_DIR/active_cell_types.txt so 02/03/04 read back the same decision —
# re-running at a different CCDC_CELL_COUNT_THRESHOLD needs no script edits.
message("\n-- Cell-type coverage filter (Erik, this threshold) --")
CELL_TYPES <- filter_cell_types_by_erik_coverage(
  ALL_CELL_TYPES, ERIK_SWEEP_BASE, ERIK_CELL_COUNT_THRESHOLD, ERIK_DISP_PREFIX,
  OUT_DIR
)

# ---------------------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------------------

# Row-wise median over finite entries only. -Inf/NaN appear in dispersion
# matrices for genes too sparsely expressed in some individuals; we summarise
# over the finite entries explicitly (NOT a silent fillna) and return NA for a
# gene with zero finite entries so it can be dropped and logged downstream.
# ANALYSIS_OK[within-file-helper]: called by collapse_to_gene_median below; R026
# cross-file rule cannot see intra-file callers.
row_median_finite <- function(mat) {
  stopifnot(is.matrix(mat))
  apply(mat, 1, function(x) {
    xf <- x[is.finite(x)]
    if (length(xf) == 0L) return(NA_real_)
    median(xf)
  })
}

# Collapse a per-sample dispersion matrix (genes x samples) to one value per
# gene, logging how many genes carried any non-finite entry and how many are
# dropped entirely (no finite entry). Returns a named numeric vector (finite).
# ANALYSIS_OK[within-file-helper]: called in both cohort-load loops below; R026
# cross-file rule cannot see intra-file callers.
collapse_to_gene_median <- function(mat, label) {
  n_gene0 <- nrow(mat)
  n_nonfinite_gene <- sum(apply(mat, 1, function(x) any(!is.finite(x))))
  med <- row_median_finite(mat)
  n_dropped <- sum(is.na(med))
  # ANALYSIS_OK[finite-filter]: drop genes with no finite dispersion entry in any
  # sample (nothing to summarise); n_dropped is counted above and logged in the
  # message() immediately below, so the drop is recorded, not silent.
  med <- med[!is.na(med)]
  message(sprintf(
    "  [%s] %d genes x %d samples; %d gene(s) had >=1 non-finite entry; "
    , label, n_gene0, ncol(mat), n_nonfinite_gene),
    sprintf("%d gene(s) dropped (no finite entry) -> %d genes retained",
            n_dropped, length(med)))
  med
}

# Drop ^RP genes from a named numeric vector; log the count removed.
# ANALYSIS_OK[within-file-helper]: called in both cohort-load loops below; R026
# cross-file rule cannot see intra-file callers.
drop_ribo <- function(vec, label) {
  is_ribo <- grepl(RIBO_PATTERN, names(vec))
  message(sprintf("  [%s] dropping %d ribosomal (^RP) gene(s): %d -> %d",
                  label, sum(is_ribo), length(vec), sum(!is_ribo)))
  vec[!is_ribo]
}

# ---------------------------------------------------------------------------
# LOAD BRENDAN COHORT
# ---------------------------------------------------------------------------
message("\n-- Brendan cohort --")
stopifnot(file.exists(BRENDAN_DISP_CSV))
brendan_mat <- as.matrix(read.csv(BRENDAN_DISP_CSV, row.names = 1,
                                  check.names = FALSE))
stopifnot(is.numeric(brendan_mat))  # fail loud if a col parsed as character
message("Loaded Brendan matrix: ", nrow(brendan_mat), " genes x ",
        ncol(brendan_mat), " columns")

stopifnot(file.exists(BRENDAN_MEAN_CSV))
brendan_mean_mat <- as.matrix(read.csv(BRENDAN_MEAN_CSV, row.names = 1,
                                       check.names = FALSE))
stopifnot(is.numeric(brendan_mean_mat))
stopifnot(identical(rownames(brendan_mat), rownames(brendan_mean_mat)),
          identical(colnames(brendan_mat), colnames(brendan_mean_mat)))
message("Loaded Brendan mean matrix: ", nrow(brendan_mean_mat), " genes x ",
        ncol(brendan_mean_mat), " columns")

# Brendan's raw file always has the full fixed 18-column design (3
# individuals x ALL_CELL_TYPES) regardless of which cell types survive
# Erik's per-threshold coverage filter (CELL_TYPES, possibly a subset) — that
# full shape is asserted here, independent of the active-subset logic below.
full_expected_cols <- as.vector(outer(BRENDAN_INDS, ALL_CELL_TYPES, paste, sep = "_"))
stopifnot(ncol(brendan_mat) == length(full_expected_cols))
stopifnot(!anyDuplicated(rownames(brendan_mat)))

# Expect (at least) the {IND}_{CellType} columns for the active CELL_TYPES.
expected_cols <- as.vector(outer(BRENDAN_INDS, CELL_TYPES, paste, sep = "_"))
missing_cols <- setdiff(expected_cols, colnames(brendan_mat))
if (length(missing_cols) > 0) {
  stop("Brendan matrix missing expected columns: ",
       paste(missing_cols, collapse = ", "))
}

brendan_by_ct <- list()
for (ct in CELL_TYPES) {
  cols <- paste(BRENDAN_INDS, ct, sep = "_")
  stopifnot(all(cols %in% colnames(brendan_mat)))
  message("Cell type ", ct, ": median across ", length(cols),
          " individual columns")
  ct_mat <- exclude_top_expressed(brendan_mat[, cols, drop = FALSE],
                                  brendan_mean_mat[, cols, drop = FALSE],
                                  paste0("brendan:", ct))
  vec <- collapse_to_gene_median(ct_mat, paste0("brendan:", ct))
  vec <- drop_ribo(vec, paste0("brendan:", ct))
  brendan_by_ct[[ct]] <- vec
}

# ---------------------------------------------------------------------------
# LOAD ERIK COHORT
# ---------------------------------------------------------------------------
message("\n-- Erik cohort --")
erik_by_ct <- list()
for (ct in CELL_TYPES) {
  f <- file.path(ERIK_SWEEP_BASE, paste0("cell_count_", ERIK_CELL_COUNT_THRESHOLD),
                 paste0(ERIK_DISP_PREFIX, "_", ct, ".csv"))
  stopifnot(file.exists(f))
  mat <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE))
  stopifnot(is.numeric(mat))  # fail loud if a col parsed as character
  # Every column name must end with the cell type (columns are {IND}_{B#}_{CT}).
  bad <- colnames(mat)[!grepl(paste0("_", ct, "$"), colnames(mat))]
  if (length(bad) > 0) {
    stop("Erik matrix for ", ct, " has columns not ending in _", ct, ": ",
         paste(head(bad, 5), collapse = ", "))
  }
  stopifnot(!anyDuplicated(rownames(mat)))
  message("Cell type ", ct, ": loaded ", nrow(mat), " genes x ", ncol(mat),
          " individual-batch samples from ", basename(f))

  mean_f <- file.path(ERIK_SWEEP_BASE, paste0("cell_count_", ERIK_CELL_COUNT_THRESHOLD),
                      paste0(ERIK_MEAN_PREFIX, "_", ct, ".csv"))
  stopifnot(file.exists(mean_f))
  mean_mat <- as.matrix(read.csv(mean_f, row.names = 1, check.names = FALSE))
  stopifnot(is.numeric(mean_mat))
  stopifnot(identical(rownames(mat), rownames(mean_mat)),
            identical(colnames(mat), colnames(mean_mat)))
  mat <- exclude_top_expressed(mat, mean_mat, paste0("erik:", ct))

  # TODO(assumption): each column is one individual-batch sample; an individual
  # present in multiple batches is weighted by its batch count. We take the
  # median across ALL sample columns (matches existing repo practice in
  # compute_and_plot_correlations.R, which treats each ind-batch as a sample).
  # Sensitivity alternative: collapse batches within individual first.
  vec <- collapse_to_gene_median(mat, paste0("erik:", ct))
  vec <- drop_ribo(vec, paste0("erik:", ct))
  erik_by_ct[[ct]] <- vec
}

# ---------------------------------------------------------------------------
# WIDE PER-COHORT MEDIAN TABLES (gene x cell type)
# ---------------------------------------------------------------------------
# ANALYSIS_OK[within-file-helper]: called just below for both cohorts; R026
# cross-file rule cannot see intra-file callers.
build_wide <- function(vec_list) {
  all_genes <- sort(unique(unlist(lapply(vec_list, names))))
  m <- matrix(NA_real_, nrow = length(all_genes), ncol = length(CELL_TYPES),
              dimnames = list(all_genes, CELL_TYPES))
  for (ct in CELL_TYPES) {
    v <- vec_list[[ct]]
    m[names(v), ct] <- v
  }
  data.frame(gene = all_genes, m, check.names = FALSE, row.names = NULL)
}

brendan_wide <- build_wide(brendan_by_ct)
erik_wide    <- build_wide(erik_by_ct)
write.csv(brendan_wide, file.path(OUT_DIR, "brendan_median_dispersion.csv"),
          row.names = FALSE)
write.csv(erik_wide, file.path(OUT_DIR, "erik_median_dispersion.csv"),
          row.names = FALSE)
message("\nWrote wide median tables: brendan ", nrow(brendan_wide),
        " genes, erik ", nrow(erik_wide), " genes")

# ---------------------------------------------------------------------------
# ALIGN MATCHED CELL TYPES (intersect gene sets)
# ---------------------------------------------------------------------------
message("\n-- Aligning matched cell types --")
counts <- data.frame(cell_type = character(), n_brendan = integer(),
                     n_erik = integer(), n_shared = integer(),
                     stringsAsFactors = FALSE)
for (ct in CELL_TYPES) {
  b <- brendan_by_ct[[ct]]
  e <- erik_by_ct[[ct]]
  shared <- intersect(names(b), names(e))
  n_shared <- length(shared)
  message(sprintf("  %s: brendan=%d, erik=%d, shared=%d",
                  ct, length(b), length(e), n_shared))
  if (n_shared == 0) {
    stop("No shared genes for ", ct, " — check gene-ID namespaces ",
         "(Brendan mixes ~376 ENSG IDs with symbols; Erik is symbol-only).")
  }
  if (n_shared < MIN_EXPECTED_SHARED_GENES) {
    warning("Only ", n_shared, " shared genes for ", ct,
            " — fewer than expected; inspect gene-ID overlap.")
  }
  df <- data.frame(gene = shared,
                   brendan_disp = b[shared],
                   erik_disp = e[shared],
                   row.names = NULL)
  stopifnot(!anyNA(df$brendan_disp), !anyNA(df$erik_disp))
  stopifnot(!anyDuplicated(df$gene))
  write.csv(df, file.path(ALIGN_DIR, paste0(ct, ".csv")), row.names = FALSE)
  counts <- rbind(counts, data.frame(cell_type = ct, n_brendan = length(b),
                                      n_erik = length(e), n_shared = n_shared))
}

write.csv(counts, file.path(OUT_DIR, "shared_gene_counts.csv"),
          row.names = FALSE)
message("\nShared-gene counts written. Total shared (summed): ",
        sum(counts$n_shared))
message("=== 01 complete ===")
