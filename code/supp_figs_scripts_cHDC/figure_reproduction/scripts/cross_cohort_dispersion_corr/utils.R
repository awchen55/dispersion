read_wide <- function(name) {
  f <- file.path(OUT_DIR, name)
  stopifnot(file.exists(f))
  m <- as.matrix(read.csv(f, check.names = FALSE, row.names = 1))
  stopifnot(all(CELL_TYPES %in% colnames(m)))
  m[, CELL_TYPES, drop = FALSE]
}

# Named constants (not smuggled function-signature defaults) for the
# cell-type coverage filter below.
MIN_INDIVIDUALS_PER_CT <- 3L        # below this, a cell type is dropped
MIN_CELL_TYPES_FOR_COMPARISON <- 2L # need >=2 surviving CTs for any
                                    # matched-vs-mismatched-CT comparison

# Determine which of all_cell_types have >= min_individuals unique
# individuals in Erik's per-cell-type dispersion file at the given
# threshold, so the pipeline can be re-run at any CCDC_CELL_COUNT_THRESHOLD
# purely by changing that env var -- no script edits needed when a stricter
# threshold leaves too few individuals in some cell type. Called once, by
# 01_load_and_align_dispersion.R (the only script that touches Erik's raw
# per-CT files); the result is persisted to out_dir/active_cell_types.txt so
# 02/03/04 can read the same decision back via read_active_cell_types()
# below instead of re-deriving it.
#
# Individual IDs are recovered from sample columns named
# "{individual}_B{batch}_{cell_type}" by stripping from the first
# "_B<digits>_" marker onward. Reads only the CSV header (nrows = 0), not
# the full file.
filter_cell_types_by_erik_coverage <- function(all_cell_types, erik_sweep_base,
                                               threshold, disp_prefix, out_dir,
                                               min_individuals = MIN_INDIVIDUALS_PER_CT) {
  keep <- logical(length(all_cell_types))
  names(keep) <- all_cell_types
  for (ct in all_cell_types) {
    f <- file.path(erik_sweep_base, paste0("cell_count_", threshold),
                  paste0(disp_prefix, "_", ct, ".csv"))
    stopifnot(file.exists(f))
    sample_cols <- colnames(read.csv(f, row.names = 1, check.names = FALSE,
                                     nrows = 0))
    individuals <- unique(sub("_B[0-9]+_.*$", "", sample_cols))
    keep[[ct]] <- length(individuals) >= min_individuals
    message(sprintf(
      "  [coverage] %s @ cell_count_%s: %d unique individual(s) in Erik's data (%s)",
      ct, threshold, length(individuals),
      if (keep[[ct]]) "kept" else sprintf("DROPPED, < %d", min_individuals)))
  }
  dropped <- all_cell_types[!keep]
  if (length(dropped) > 0) {
    message("Cell type(s) dropped for insufficient Erik individual coverage ",
            "at threshold ", threshold, ": ", paste(dropped, collapse = ", "))
  }
  active <- all_cell_types[keep]
  stopifnot(length(active) >= MIN_CELL_TYPES_FOR_COMPARISON)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(active, file.path(out_dir, "active_cell_types.txt"))
  active
}

# Read back the cell-type list 01_load_and_align_dispersion.R determined
# survives Erik's per-CT individual-coverage filter at this threshold (see
# filter_cell_types_by_erik_coverage above).
read_active_cell_types <- function(out_dir) {
  f <- file.path(out_dir, "active_cell_types.txt")
  if (!file.exists(f)) {
    stop("Missing ", f, " — run 01_load_and_align_dispersion.R first ",
         "(same CCDC_CELL_COUNT_THRESHOLD).")
  }
  cts <- readLines(f)
  stopifnot(length(cts) >= MIN_CELL_TYPES_FOR_COMPARISON)
  cts
}

# Alexander's LOWESS-derived Memento mean-correction threshold (natural-log
# relative-abundance scale, same scale/value as
# analysis/dd_excl_topexpressed/scripts/01_exclude_topexpressed_genes.R):
# dispersion for genes with mean expression above this value is not reliably
# mean-corrected by Memento. Applies to any cohort quantified by Memento, not
# just Brendan's discovery cohort — used identically for both cohorts here.
MEAN_EXCL_THRESHOLD <- -6

# Exclude genes whose mean expression exceeds MEAN_EXCL_THRESHOLD in ANY
# sample column of mean_mat (generalizes the dd_excl_topexpressed
# "any-individual" rule to "any-sample", since Erik's per-CT matrices carry
# more than one sample per individual). disp_mat and mean_mat must have
# identical dimnames (same genes, same sample columns), asserted here so a
# misaligned mean matrix fails loudly rather than silently mis-excluding.
# Returns disp_mat subset to retained genes; logs excluded/retained counts.
exclude_top_expressed <- function(disp_mat, mean_mat, label,
                                  threshold = MEAN_EXCL_THRESHOLD) {
  stopifnot(identical(rownames(disp_mat), rownames(mean_mat)),
            identical(colnames(disp_mat), colnames(mean_mat)))
  stopifnot(!anyNA(mean_mat))  # NaN/NA would silently propagate through any()
  excluded <- apply(mean_mat, 1, function(x) any(x > threshold))
  message(sprintf("  [%s] mean-threshold excl (Memento mean > %s in any sample): %d excluded / %d retained",
                  label, threshold, sum(excluded), sum(!excluded)))
  disp_mat[!excluded, , drop = FALSE]
}
