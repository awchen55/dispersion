stopifnot("Run from the repo root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(mashr)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

source("analysis/dqtl_genotype_dosage_plots/scripts/utils.R")

RUN_ID        <- "20260925_143116"
SOURCE_RUN_ID <- "20260924_150719"
DQTL_MAPPING_DIR <- "analysis/dqtl_mapping"
MASH_RDS      <- file.path(DQTL_MAPPING_DIR, "outputs/mash/m2_dQTL_250kb_20260925_143116.rds")
GENO_DIR      <- file.path(DQTL_MAPPING_DIR, "inputs/genotypes")
EXPECTED_N_CELL_TYPES <- 4L

# Matches the primary threshold used throughout dqtl_mapping/dqtl_egene_overlap.
LFSR_THRESHOLD <- 0.1

# This run's own q=0.2 Memento output (threshold=150, the only threshold this
# rerun computed) -- NOT Erik's q=0.25 matrices. Raw (non-"new_"-reformatted)
# CSVs, same as dqtl_genotype_dosage_plots' own scripts read: `gene` column +
# raw "sg^{ind}_{batch}^{celltype}" sample columns, matching get_cohort()'s
# mat_cols convention.
CELL_COUNT_DIR <- "analysis/dqtl_mapping_q02/outputs/cell_count_threshold_sweep_q02/cell_count_150"

OUT_DIR  <- "analysis/dqtl_mapping_q02/outputs/genotype_dosage_plots"
OUT_DATA <- file.path(OUT_DIR, "data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DATA, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load mash; build gene -> cell types (dGene at LFSR_THRESHOLD) membership.
# ---------------------------------------------------------------------------

stopifnot(file.exists(MASH_RDS))
m2       <- readRDS(MASH_RDS)
lfsr_mat <- m2$result$lfsr
stopifnot(is.matrix(lfsr_mat))

cell_types <- sub("^dQTL_(.+)_AllChroms\\.out$", "\\1", colnames(lfsr_mat))
stopifnot("Expected 4 cell types for RUN_ID 20260925_143116" = length(cell_types) == EXPECTED_N_CELL_TYPES)

rn_genes <- sub("_[0-9]+:.*$", "", rownames(lfsr_mat))
stopifnot("mash rownames are not 1 row per gene" = !anyDuplicated(rn_genes))

dgene_sets <- lapply(cell_types, function(ct) {
  col <- paste0("dQTL_", ct, "_AllChroms.out")
  rn_genes[lfsr_mat[, col] < LFSR_THRESHOLD]
})
names(dgene_sets) <- cell_types

all_dgenes <- sort(unique(unlist(dgene_sets)))
gene_celltypes <- setNames(
  lapply(all_dgenes, function(g) cell_types[vapply(dgene_sets, function(s) g %in% s, logical(1))]),
  all_dgenes
)
message(length(all_dgenes), " unique dGenes (lfsr < ", LFSR_THRESHOLD, ") across ",
        length(cell_types), " cell types")

# ---------------------------------------------------------------------------
# Per-cell-type caches: cohort, this run's own mean/dispersion matrices.
# ---------------------------------------------------------------------------

cohorts <- setNames(lapply(cell_types, function(ct) get_cohort(DQTL_MAPPING_DIR, SOURCE_RUN_ID, ct)),
                     cell_types)

mean_mats <- setNames(lapply(cell_types, function(ct) {
  read.csv(file.path(CELL_COUNT_DIR, paste0("mean_by_ind_", ct, ".csv")),
           check.names = FALSE, row.names = 1)
}), cell_types)

disp_mats <- setNames(lapply(cell_types, function(ct) {
  read.csv(file.path(CELL_COUNT_DIR, paste0("dispersion_by_ind_", ct, ".csv")),
           check.names = FALSE, row.names = 1)
}), cell_types)

# ---------------------------------------------------------------------------
# Per (gene, cell type) data frame: dosage + Memento mean + dispersion, for
# that cell type's dqtl_mapping cohort.
# ---------------------------------------------------------------------------

build_gene_celltype_df <- function(gene, ct, dosage_by_donor) {
  cohort <- cohorts[[ct]]
  donor_batch <- cohort$donor_batch
  donor       <- cohort$donor
  mat_cols    <- cohort$mat_cols
  stopifnot("A cohort donor is missing a genotype column for this SNP" =
              all(donor %in% names(dosage_by_donor)))

  mean_mat <- mean_mats[[ct]]
  disp_mat <- disp_mats[[ct]]
  gene_row <- which(mean_mat$gene == gene)
  stopifnot("GENE not found (or duplicated) in this run's q=0.2 cell_count matrix for this cell type" =
              length(gene_row) == 1)
  stopifnot("A dqtl_mapping cohort sample is missing from this run's q=0.2 cell_count matrix" =
              all(mat_cols %in% colnames(mean_mat)))
  mean_vals <- as.numeric(mean_mat[gene_row, mat_cols])
  disp_vals <- as.numeric(disp_mat[gene_row, mat_cols])
  stopifnot("Non-finite Memento mean/dispersion value(s) for this gene/cohort" =
              all(is.finite(mean_vals)), all(is.finite(disp_vals)))

  data.frame(
    cell_type   = ct,
    donor_batch = donor_batch,
    dosage      = factor(dosage_by_donor[donor], levels = c("0", "1", "2")),
    mean_memento = mean_vals,
    dispersion   = disp_vals
  )
}

# ---------------------------------------------------------------------------
# Per (gene, cell type): stacked mean (top) / dispersion (bottom) figure.
# ---------------------------------------------------------------------------

n_written <- 0L
for (gene in all_dgenes) {
  snp_info <- resolve_gene_snp(lfsr_mat, gene)
  stopifnot("dGene missing from mash rownames -- should be impossible by construction" =
              !is.null(snp_info))
  chr <- snp_info$chr; ref <- snp_info$ref; alt <- snp_info$alt
  geno_snp_id <- snp_info$geno_snp_id
  dosage_by_donor <- lookup_genotype_dosage(GENO_DIR, chr, geno_snp_id)

  for (ct in gene_celltypes[[gene]]) {
    gene_ct_df <- build_gene_celltype_df(gene, ct, dosage_by_donor)

    out_csv <- file.path(OUT_DATA, paste0(gene, "_", ct, "_dosage_data_", RUN_ID, ".csv"))
    write.csv(gene_ct_df, out_csv, row.names = FALSE)

    p_mean <- dosage_boxplot_panel(gene_ct_df, "mean_memento", "Mean expression", gene, ref, alt)
    p_disp <- dosage_boxplot_panel(gene_ct_df, "dispersion", "Dispersion", NULL, ref, alt)

    combined <- (p_mean / p_disp) +
      plot_annotation(
        title = paste0(gsub("_", " ", ct), ", SNP ", geno_snp_id),
        theme = theme(plot.title = element_text(size = 6, family = "Arial"))
      )

    out_pdf <- file.path(OUT_DIR, paste0(gene, "_", ct, "_dosage_", RUN_ID, ".pdf"))
    ggsave(out_pdf, combined, width = 1.8, height = 4.0, units = "in", device = cairo_pdf)

    n_written <- n_written + 1L
  }
  message(gene, ": ", length(gene_celltypes[[gene]]), " cell type(s) -- ",
          paste(gene_celltypes[[gene]], collapse = ", "))
}

message("\nDone: ", n_written, " (gene, cell type) figures written under ", OUT_DIR)
