stopifnot(file.exists("CLAUDE.md"))  # must run from project root

library(tidyverse)

# ---------------------------------------------------------------------------
# CONFIG (env-var overridable -- see DQTL_MAPPING.md "Run configuration";
# defaults match the original script's)
# ---------------------------------------------------------------------------
CELL_COUNT_THRESHOLD <- as.integer(Sys.getenv("DQTL_CELL_COUNT_THRESHOLD", unset = "150"))
stopifnot("This q=0.2 rerun only has Memento output at cell_count_threshold=150 -- see analysis/dqtl_mapping_q02/scripts/01_memento_q02_threshold150.py" =
            CELL_COUNT_THRESHOLD == 150L)
SD_CUTOFF_FRAC       <- as.numeric(Sys.getenv("DQTL_SD_CUTOFF_FRAC", unset = "0.30"))
CELL_TYPES <- strsplit(
  Sys.getenv("DQTL_CELL_TYPES",
             unset = "Activated_fibroblasts,Cardiac_fibroblasts,Cardiomyocytes"),
  ","
)[[1]]
GENE_SELECTION_STRATEGY <- Sys.getenv("DQTL_GENE_SELECTION_STRATEGY", unset = "union")
stopifnot("DQTL_GENE_SELECTION_STRATEGY must be 'union' or 'intersection'" =
            GENE_SELECTION_STRATEGY %in% c("union", "intersection"))
RUN_ID <- Sys.getenv("DQTL_RUN_ID", unset = format(Sys.time(), "%Y%m%d_%H%M%S"))
message("Run ID: ", RUN_ID)

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------
PROJECT_ROOT <- getwd()
# q=0.2 rerun's own Memento output (reformatted by 00_reformat_memento_q02.R),
# in place of the original script's
# analysis/qtl_analyses/outputs/cell_count_threshold_sweep/cell_count_{N}/ --
# see .living/decisions.md, 2026-09-24.
ERIK_SWEEP <- file.path(PROJECT_ROOT,
                        "analysis/dqtl_mapping_q02/outputs/cell_count_threshold_sweep_q02",
                        paste0("cell_count_", CELL_COUNT_THRESHOLD))
ERIK_META_DIR  <- "/project/gilad/emcintire/chemo_project_2023/pipeline/pseudobulk_and_normalize/metadata"
# Filtered raw pseudobulk counts (pre-TMM/INV) -- library size isn't saved
# anywhere by Erik's pipeline, so it's derived here the same way as the
# original script: colSums() of these raw counts. Unaffected by Memento's q,
# so reused as-is from Erik's directory (read-only).
ERIK_RAW_COUNTS_DIR <- "/project/gilad/emcintire/chemo_project_2023/pipeline/pseudobulk_and_normalize/filtered_raw_counts"
PROTEIN_CODING <- "/project/gilad/emcintire/chemo_project_2023/pipeline/pseudobulk_and_normalize/protein_coding_genes_19268_symbols_27MAR25_EM.tsv"
OUT_DISP       <- file.path(PROJECT_ROOT, "analysis/dqtl_mapping/inputs/dispersion_matrices", RUN_ID)
OUT_META       <- file.path(PROJECT_ROOT, "analysis/dqtl_mapping/inputs/metadata", RUN_ID)
OUT_GENES      <- file.path(PROJECT_ROOT, "analysis/dqtl_mapping/inputs/gene_sets", RUN_ID)
dir.create(OUT_DISP,  recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_META,  recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_GENES, recursive = TRUE, showWarnings = FALSE)

stopifnot("ERIK_SWEEP directory not found -- run 00_reformat_memento_q02.R (and 01_memento_q02_threshold150.py before it) first" =
            dir.exists(ERIK_SWEEP))
stopifnot(dir.exists(ERIK_META_DIR), dir.exists(ERIK_RAW_COUNTS_DIR),
          file.exists(PROTEIN_CODING))

# ---------------------------------------------------------------------------
# LOAD INPUTS
# ---------------------------------------------------------------------------
protein_coding <- read.table(PROTEIN_CODING, header = TRUE,
                              stringsAsFactors = FALSE)$symbol
stopifnot(length(protein_coding) > 0)
message("Protein-coding genes loaded: ", length(protein_coding))

# ---------------------------------------------------------------------------
# STEP 1: Load all cell types' q=0.2 matrices; find genes finite in all
#         samples per cell type; take the cross-cell-type intersection
# ---------------------------------------------------------------------------
erik_mats <- lapply(setNames(CELL_TYPES, CELL_TYPES), function(ct) {
  f <- file.path(ERIK_SWEEP, paste0("new_dispersion_by_ind_", ct, ".csv"))
  stopifnot(file.exists(f))
  mat <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE))
  colnames(mat) <- gsub('^"|"$', "", colnames(mat))
  message(ct, ": ", nrow(mat), " genes, ", ncol(mat), " samples")
  mat
})

finite_genes <- lapply(CELL_TYPES, function(ct) {
  mat <- erik_mats[[ct]]
  keep <- apply(mat, 1, function(x) all(is.finite(x)))
  message(ct, ": ", sum(keep), " genes with all-finite values (",
          sum(!keep), " dropped for -Inf)")
  rownames(mat)[keep]
})
names(finite_genes) <- CELL_TYPES

genes_all_ct <- Reduce(intersect, finite_genes)
message("\nGenes finite in all ", length(CELL_TYPES), " cell types: ", length(genes_all_ct))

# ---------------------------------------------------------------------------
# STEP 2: Filter to protein-coding
# ---------------------------------------------------------------------------
genes_pc <- intersect(genes_all_ct, protein_coding)
message("After protein-coding filter: ", length(genes_pc),
        " (dropped ", length(genes_all_ct) - length(genes_pc), ")")

# ---------------------------------------------------------------------------
# STEP 3: per cell type -- rank by inter-individual SD of dispersion, take the
# top SD_CUTOFF_FRAC, computed directly from erik_mats (this run's own
# dispersion-by-individual matrices, the same ones used as the QTL-mapping
# phenotype) -- matches the original script's method (see its header,
# .living/decisions.md 2026-07-14), just against q=0.2 data.
# ---------------------------------------------------------------------------
top30_per_ct <- lapply(setNames(CELL_TYPES, CELL_TYPES), function(ct) {
  mat_pc  <- erik_mats[[ct]][genes_pc, , drop = FALSE]
  sd_vals <- apply(mat_pc, 1, sd)
  n_keep  <- ceiling(length(sd_vals) * SD_CUTOFF_FRAC)
  # ANALYSIS_OK[top-n-selection]: seq_len(n_keep) takes the top n_keep entries
  # of this same already-sorted vector (rank-based selection), not a
  # cross-structure positional alignment -- n_keep is computed from and
  # applied to the identical sd_vals vector.
  genes_use <- names(sort(sd_vals, decreasing = TRUE))[seq_len(n_keep)]
  message(ct, ": top ", SD_CUTOFF_FRAC * 100, "% by inter-individual SD (q=0.2 data) = ",
          length(genes_use), " genes")
  genes_use
})

message("\nGenes shared by all ", length(CELL_TYPES), " per-cell-type top-",
        SD_CUTOFF_FRAC * 100, "% sets: ", length(Reduce(intersect, top30_per_ct)),
        " (of ", length(Reduce(union, top30_per_ct)), " total across the union)")

# ---------------------------------------------------------------------------
# Resolve which gene set each cell type actually gets tested on (identical
# logic to the original script -- see its header)
# ---------------------------------------------------------------------------
gene_set_per_ct <- if (GENE_SELECTION_STRATEGY == "union") {
  union_genes <- Reduce(union, top30_per_ct)
  message("\nGene selection strategy: union -- testing the same ",
          length(union_genes), "-gene union of all per-cell-type top-",
          SD_CUTOFF_FRAC * 100, "% sets in every cell type.")
  setNames(replicate(length(CELL_TYPES), union_genes, simplify = FALSE), CELL_TYPES)
} else {
  intersection_genes <- Reduce(intersect, top30_per_ct)
  stopifnot("intersection of per-cell-type top-SD_CUTOFF_FRAC sets is empty -- no genes shared a high-SD signature across all cell types" =
              length(intersection_genes) > 0)
  message("\nGene selection strategy: intersection -- testing the same ",
          length(intersection_genes), "-gene intersection of all per-cell-type top-",
          SD_CUTOFF_FRAC * 100, "% sets in every cell type.")
  setNames(replicate(length(CELL_TYPES), intersection_genes, simplify = FALSE), CELL_TYPES)
}

# ---------------------------------------------------------------------------
# Output: write each cell type's resolved gene set (identical to original)
# ---------------------------------------------------------------------------
for (ct in CELL_TYPES) {
  message("\n=== ", ct, " ===")

  sub_mat <- erik_mats[[ct]][gene_set_per_ct[[ct]], , drop = FALSE]
  # ANALYSIS_OK[magic-threshold]: 100 is a sanity floor guarding against a
  # degenerate/near-empty gene set, not a scientifically meaningful cutoff --
  # the actual selection threshold is DQTL_SD_CUTOFF_FRAC above, recorded in
  # prep_config.txt. Matches the original script's identical floor.
  stopifnot(nrow(sub_mat) > 100, all(is.finite(sub_mat)))

  # Write MatrixEQTL dispersion TSV
  out_disp <- file.path(OUT_DISP, paste0(ct, "_dispersion.tsv"))
  con <- file(out_disp, "w")
  writeLines(paste(colnames(sub_mat), collapse = "\t"), con)
  write.table(sub_mat, con, sep = "\t", quote = FALSE, col.names = FALSE)
  close(con)
  message("  Written: ", out_disp)

  writeLines(rownames(sub_mat), file.path(OUT_GENES, paste0(ct, "_genes.txt")))

  # -------------------------------------------------------------------------
  # Extract sex metadata aligned to dispersion sample IDs (identical logic;
  # Erik's Sex metadata doesn't depend on Memento's q)
  # -------------------------------------------------------------------------
  meta_f <- file.path(ERIK_META_DIR,
                       paste0(ct, "_Control_PB_FILT_TMM_INV_metadata.tsv"))
  stopifnot(file.exists(meta_f))
  meta <- read.table(meta_f, header = TRUE, sep = "\t",
                      row.names = 1, check.names = FALSE)

  disp_ind_batch <- sub(paste0("_", ct, "$"), "", colnames(sub_mat))
  meta_ind_batch <- sub("_Control$", "", colnames(meta))

  idx <- match(disp_ind_batch, meta_ind_batch)

  fallback <- is.na(idx)
  if (any(fallback)) {
    disp_ind_only <- sub("_B[0-9]+$", "", disp_ind_batch)
    meta_ind_only <- sub("_B[0-9]+$", "", meta_ind_batch)
    idx[fallback] <- match(disp_ind_only[fallback], meta_ind_only)
    if (anyNA(idx)) {
      miss <- disp_ind_batch[is.na(idx)]
      stop(ct, ": metadata missing for: ", paste(miss, collapse = ", "))
    }
    message("  Resolved ", sum(fallback), " sample(s) via individual-only match: ",
            paste(disp_ind_batch[fallback], collapse = ", "))
  }

  sex_vals <- meta["Sex", idx, drop = FALSE]
  colnames(sex_vals) <- colnames(sub_mat)

  # -------------------------------------------------------------------------
  # Extract library size aligned to dispersion sample IDs (identical logic;
  # raw pseudobulk counts don't depend on Memento's q)
  # -------------------------------------------------------------------------
  raw_f <- file.path(ERIK_RAW_COUNTS_DIR, paste0("AllConditions_", ct, "_PB_FILT.txt"))
  stopifnot(file.exists(raw_f))
  raw_counts  <- as.matrix(read.table(raw_f, header = TRUE, row.names = 1,
                                       check.names = FALSE, sep = "\t"))
  raw_control <- raw_counts[, grepl("_Control$", colnames(raw_counts)), drop = FALSE]
  lib_size_by_ind_batch <- setNames(colSums(raw_control),
                                     sub("_Control$", "", colnames(raw_control)))
  stopifnot(
    "raw-count Control samples don't match Sex-metadata Control samples" =
      setequal(names(lib_size_by_ind_batch), meta_ind_batch)
  )
  lib_size_vals <- lib_size_by_ind_batch[meta_ind_batch[idx]]
  names(lib_size_vals) <- colnames(sub_mat)

  out_meta <- file.path(OUT_META, paste0(ct, "_metadata_sex.tsv"))
  covariate_df <- data.frame(
    id = c("Sex", "library_size"),
    rbind(as.numeric(sex_vals), as.numeric(lib_size_vals)),
    check.names = FALSE
  )
  colnames(covariate_df)[-1] <- colnames(sub_mat)
  write.table(covariate_df, out_meta, sep = "\t", quote = FALSE, row.names = FALSE)
  message("  Written metadata: ", out_meta)
}

# ---------------------------------------------------------------------------
# Write config provenance
# ---------------------------------------------------------------------------
writeLines(
  c(paste0("cell_count_threshold: ", CELL_COUNT_THRESHOLD),
    paste0("sd_cutoff_frac: ", SD_CUTOFF_FRAC),
    paste0("gene_selection_strategy: ", GENE_SELECTION_STRATEGY),
    paste0("n_cell_types: ", length(CELL_TYPES)),
    paste0("cell_types: ", paste(CELL_TYPES, collapse = ",")),
    paste0("n_genes_selected: ", paste(lengths(gene_set_per_ct[CELL_TYPES]), collapse = ",")),
    "memento_capture_rate: 0.2",
    paste0("created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))),
  file.path(OUT_GENES, "prep_config.txt")
)
message("Written config: ", file.path(OUT_GENES, "prep_config.txt"))

message("\nDone. Inputs written to analysis/dqtl_mapping/inputs/*/", RUN_ID, "/")
