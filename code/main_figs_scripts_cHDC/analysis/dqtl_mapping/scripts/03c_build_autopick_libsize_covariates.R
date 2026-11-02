stopifnot(file.exists("CLAUDE.md"))

DQTL_BASE    <- "analysis/dqtl_mapping"
META_BASE    <- file.path(DQTL_BASE, "inputs/metadata")
DISP_BASE    <- file.path(DQTL_BASE, "inputs/dispersion_matrices")
GENESET_BASE <- file.path(DQTL_BASE, "inputs/gene_sets")
RESULTS_BASE <- file.path(DQTL_BASE, "outputs/covariate_optimization")

RUN_ID <- Sys.getenv("DQTL_RUN_ID", unset = format(Sys.time(), "%Y%m%d_%H%M%S"))
message("Run ID: ", RUN_ID)

SOURCE_RUN_ID <- Sys.getenv("DQTL_SOURCE_RUN_ID", unset = "")
stopifnot("DQTL_SOURCE_RUN_ID must be set to the RUN_ID whose inputs/covariate snapshot to reuse" =
            nzchar(SOURCE_RUN_ID))

SOURCE_META_DIR <- file.path(META_BASE, SOURCE_RUN_ID)
stopifnot("inputs/metadata/{DQTL_SOURCE_RUN_ID}/ not found" = dir.exists(SOURCE_META_DIR))
SOURCE_DISP_DIR <- file.path(DISP_BASE, SOURCE_RUN_ID)
stopifnot("inputs/dispersion_matrices/{DQTL_SOURCE_RUN_ID}/ not found" = dir.exists(SOURCE_DISP_DIR))
SOURCE_GENESET_DIR <- file.path(GENESET_BASE, SOURCE_RUN_ID)
stopifnot("inputs/gene_sets/{DQTL_SOURCE_RUN_ID}/ not found" = dir.exists(SOURCE_GENESET_DIR))

# ---- Symlink this RUN_ID's dispersion_matrices/gene_sets to the source RUN_ID
# so 04_map_dqtls.R's unmodified RUN_ID-versioned input resolution works
# unchanged. Relative symlink (target is a sibling directory), created only if
# absent so re-running this script for an already-set-up RUN_ID is a no-op here.
link_to_source <- function(base_dir) {
  link_path <- file.path(base_dir, RUN_ID)
  if (dir.exists(link_path) || file.exists(link_path)) {
    stopifnot("RUN_ID input symlink already exists but points elsewhere" =
                identical(normalizePath(link_path), normalizePath(file.path(base_dir, SOURCE_RUN_ID))))
    message("  ", link_path, " already -> ", SOURCE_RUN_ID)
  } else {
    stopifnot(file.symlink(SOURCE_RUN_ID, link_path))
    message("  Linked ", link_path, " -> ", SOURCE_RUN_ID)
  }
}
link_to_source(DISP_BASE)
link_to_source(GENESET_BASE)

OUT_DIR <- file.path(RESULTS_BASE, RUN_ID)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CELL_TYPES <- strsplit(Sys.getenv("DQTL_CELL_TYPES", unset = ""), ",")[[1]]
stopifnot("DQTL_CELL_TYPES must be set explicitly (no default)" = length(CELL_TYPES) > 0)

no_libsize_str <- Sys.getenv("DQTL_NO_LIBSIZE_CELL_TYPES", unset = "")
NO_LIBSIZE_CTS <- if (nzchar(no_libsize_str)) strsplit(no_libsize_str, ",")[[1]] else character(0)
stopifnot("DQTL_NO_LIBSIZE_CELL_TYPES entries must all be in DQTL_CELL_TYPES" =
            all(NO_LIBSIZE_CTS %in% CELL_TYPES))

# ---- Resolve per-cell-type PC count (see header for the two paths) ----
force_pc_str <- Sys.getenv("DQTL_FORCE_PC_COUNT", unset = "")
if (nzchar(force_pc_str)) {
  # ANALYSIS_OK[warning-suppression]: as.integer() on a non-numeric env var
  # emits "NAs introduced by coercion" -- narrowly suppressed here because the
  # very next line turns that same NA into a loud stopifnot() with an
  # actionable message, so nothing is silently swallowed.
  force_k <- suppressWarnings(as.integer(force_pc_str))
  stopifnot("DQTL_FORCE_PC_COUNT must be a non-negative integer" = !is.na(force_k), force_k >= 0)
  AUTOPICK_PCS <- setNames(rep(force_k, length(CELL_TYPES)), CELL_TYPES)
  message("Forcing ", force_k, " PC(s) for every cell type (DQTL_FORCE_PC_COUNT set) -- ",
          "not reading source RUN_ID ", SOURCE_RUN_ID, "'s auto-pick suggestion")
} else {
  source_manual_f <- file.path(RESULTS_BASE, SOURCE_RUN_ID, "manual_pc_selection.txt")
  stopifnot(
    "Source RUN_ID has already been through manual PC review (manual_pc_selection.txt exists) -- its metadata_ready.tsv now holds the manual choice, not 03's auto-pick. This script only supports a source RUN_ID whose auto-pick is still on disk unreviewed (or set DQTL_FORCE_PC_COUNT to bypass reading it entirely)." =
      !file.exists(source_manual_f)
  )
  AUTOPICK_PCS <- setNames(vapply(CELL_TYPES, function(ct) {
    ready_f <- file.path(RESULTS_BASE, SOURCE_RUN_ID, paste0(ct, "_metadata_ready.tsv"))
    stopifnot("Source RUN_ID's auto-pick metadata_ready.tsv missing for a cell type -- run 03_optimize_covariates.R for it first" =
                file.exists(ready_f))
    ids <- read.table(ready_f, header = TRUE, sep = "\t", check.names = FALSE)$id
    sum(grepl("^PC[0-9]+$", ids))
  }, integer(1)), CELL_TYPES)
  message("Read auto-pick PC counts from source RUN_ID ", SOURCE_RUN_ID,
          "'s own (unreviewed) metadata_ready.tsv snapshots")
}

model_lines  <- character(0)
manual_lines <- character(0)

for (ct in CELL_TYPES) {
  k <- AUTOPICK_PCS[[ct]]
  use_libsize <- !(ct %in% NO_LIBSIZE_CTS)

  meta_f <- file.path(SOURCE_META_DIR, paste0(ct, "_metadata_sex_pcs.tsv"))
  stopifnot(file.exists(meta_f))
  meta_df <- read.table(meta_f, header = TRUE, sep = "\t", check.names = FALSE)
  stopifnot("library_size" %in% meta_df$id, "Sex" %in% meta_df$id)

  rows_keep <- c("Sex", if (use_libsize) "library_size", if (k > 0) paste0("PC", seq_len(k)))
  meta_ready     <- meta_df[meta_df$id %in% rows_keep, , drop = FALSE]
  meta_ready_idx <- match(rows_keep, meta_ready$id)
  stopifnot("Chosen covariate rows missing from the source snapshot" = !anyNA(meta_ready_idx))
  # ANALYSIS_OK[self-assignment]: reorder only, not a narrowing filter -- the
  # stopifnot immediately above guarantees meta_ready_idx has no NA, so every
  # name in rows_keep has a matching row here. Same pattern as
  # 03_optimize_covariates.R/03b_apply_manual_pcs.R.
  meta_ready <- meta_ready[meta_ready_idx, ]

  ready_out <- file.path(OUT_DIR, paste0(ct, "_metadata_ready.tsv"))
  write.table(meta_ready, ready_out, sep = "\t", quote = FALSE, row.names = FALSE)

  model_desc <- paste(c("Sex", if (use_libsize) "library_size", if (k > 0) paste0("PC1-", k)), collapse = "+")
  model_lines  <- c(model_lines, paste0(ct, "\t", model_desc))
  manual_lines <- c(manual_lines, paste0(ct, "\t", k))
  message(ct, ": Sex + ", if (use_libsize) "library_size + " else "", k,
          " PC(s) (source RUN_ID ", SOURCE_RUN_ID, ") -> ", ready_out)
}

# Written so 04_map_dqtls.R's PC-count provenance gate passes (it only checks
# PC count, so the added library_size covariate doesn't affect that check).
writeLines(manual_lines, file.path(OUT_DIR, "manual_pc_selection.txt"))

# Durable record of the actual covariate model used -- manual_pc_selection.txt
# alone doesn't capture the added library_size covariate.
pc_source_note <- if (nzchar(force_pc_str)) {
  paste0("# PC count forced to ", force_k, " for every cell type (DQTL_FORCE_PC_COUNT)")
} else {
  paste0("# Auto-pick PC counts reused from source RUN_ID ", SOURCE_RUN_ID,
         "'s own 03_optimize_covariates.R sweep")
}
libsize_note <- if (length(NO_LIBSIZE_CTS) == 0) {
  "# library_size added for every cell type, per dqtl_batch_structure F-019 (.living/findings/dd-qtl-mapping.md)"
} else {
  paste0("# library_size added per dqtl_batch_structure F-019 (.living/findings/dd-qtl-mapping.md), EXCEPT for: ",
         paste(NO_LIBSIZE_CTS, collapse = ", "), " (DQTL_NO_LIBSIZE_CELL_TYPES)")
}
writeLines(
  c(paste0("# Covariate model for RUN_ID ", RUN_ID),
    pc_source_note,
    libsize_note,
    model_lines),
  file.path(OUT_DIR, "covariate_model.txt")
)

message("\nDone. Covariate model written for RUN_ID ", RUN_ID, " (source: ", SOURCE_RUN_ID, ")")
