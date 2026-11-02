stopifnot(file.exists("CLAUDE.md"))

DQTL_BASE <- "analysis/dqtl_mapping"
RUN_ID    <- Sys.getenv("DQTL_RUN_ID", unset = format(Sys.time(), "%Y%m%d_%H%M%S"))
message("Run ID: ", RUN_ID)
DISP_DIR  <- file.path(DQTL_BASE, "inputs/dispersion_matrices", RUN_ID)
META_DIR  <- file.path(DQTL_BASE, "inputs/metadata", RUN_ID)
OUT_DIR   <- file.path(DQTL_BASE, "outputs/dispersion_pcs", RUN_ID)
MAX_N_PCS <- 20

stopifnot("inputs/dispersion_matrices/{RUN_ID}/ not found -- run 01_prepare_dispersion_inputs.R with the same DQTL_RUN_ID first" =
            dir.exists(DISP_DIR))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Env-var overridable to match 01's config (see DQTL_MAPPING.md "Run
# configuration"); defaults match 01's defaults.
CELL_TYPES <- strsplit(
  Sys.getenv("DQTL_CELL_TYPES",
             unset = "Activated_fibroblasts,Cardiac_fibroblasts,Cardiomyocytes"),
  ","
)[[1]]

for (ct in CELL_TYPES) {
  message("=== ", ct, " ===")

  # ---- Read dispersion matrix (MatrixEQTL format) ----
  # Line 1: sample IDs (no leading field); lines 2+: gene\tval1\tval2...
  disp_f  <- file.path(DISP_DIR, paste0(ct, "_dispersion.tsv"))
  stopifnot(file.exists(disp_f))
  samples <- strsplit(readLines(disp_f, n = 1), "\t")[[1]]
  mat     <- read.table(disp_f, header = FALSE, skip = 1, sep = "\t",
                         row.names = 1, check.names = FALSE)
  colnames(mat) <- samples
  mat <- as.matrix(mat)
  stopifnot(ncol(mat) == length(samples), all(is.finite(mat)))
  message("  Matrix: ", nrow(mat), " genes x ", ncol(mat), " samples")

  # ---- PCA: samples as rows, genes as columns ----
  pca_res <- prcomp(t(mat), center = TRUE, scale. = FALSE)

  # Cap PCs at the smaller of MAX_N_PCS and the matrix's actual rank --
  # centered PCA on n samples yields at most n-1 non-trivial components, so
  # cell types with <=20 samples (e.g. Mesenchymal_cells, Epithelial_cells)
  # would otherwise index past the end of pca_res$x ("subscript out of
  # bounds"). See .living/decisions.md, 2026-07-15.
  n_pcs <- min(MAX_N_PCS, ncol(mat) - 1)
  if (n_pcs < MAX_N_PCS) {
    message("  Capping at ", n_pcs, " PCs (", ncol(mat), " samples < ", MAX_N_PCS + 1, ")")
  }

  # Save variance explained
  var_exp <- summary(pca_res)$importance["Proportion of Variance", 1:n_pcs]
  write.csv(
    data.frame(PC = paste0("PC", 1:n_pcs), variance_explained = var_exp),
    file.path(OUT_DIR, paste0(ct, "_variance_explained.csv")),
    row.names = FALSE
  )
  message("  PC1 explains ", round(var_exp[1] * 100, 1), "% of variance")

  # ---- Append PCs to sex metadata (in Brendan's inputs/metadata/) ----
  meta_f <- file.path(META_DIR, paste0(ct, "_metadata_sex.tsv"))
  stopifnot(file.exists(meta_f))
  meta <- read.table(meta_f, header = TRUE, sep = "\t", check.names = FALSE)
  # Columns: id, sample_1, sample_2, ...  |  id column contains: "Sex"

  sample_cols <- colnames(meta)[-1]
  stopifnot(identical(sample_cols, samples))

  # PC scores matrix: n_pcs x N_samples
  pc_scores <- t(pca_res$x[, 1:n_pcs, drop = FALSE])
  pc_df <- data.frame(id = paste0("PC", 1:n_pcs), pc_scores,
                       check.names = FALSE, stringsAsFactors = FALSE)
  colnames(pc_df)[-1] <- samples

  meta_with_pcs <- rbind(meta, pc_df)

  out_f <- file.path(META_DIR, paste0(ct, "_metadata_sex_pcs.tsv"))
  write.table(meta_with_pcs, out_f, sep = "\t", quote = FALSE, row.names = FALSE)
  message("  Written: ", out_f)
}

message("\nDone.")
