stopifnot(file.exists("CLAUDE.md"))

library(MatrixEQTL)
library(Matrix)
library(doParallel)
library(foreach)

registerDoParallel(cores = as.integer(Sys.getenv("SLURM_NTASKS_PER_NODE",
                                                 unset = parallel::detectCores() - 1)))

DQTL_BASE    <- "analysis/dqtl_mapping"
GENO_DIR     <- file.path(DQTL_BASE, "inputs/genotypes")
SNPPOS_DIR   <- file.path(DQTL_BASE, "inputs/snp_positions")
GENEPOS_FILE <- file.path(DQTL_BASE, "inputs/gene_positions/concise_GRCh38-2020_gtf_strand_considered.csv")
RESULTS_BASE <- file.path(DQTL_BASE, "outputs/covariate_optimization")

RUN_ID  <- Sys.getenv("DQTL_RUN_ID", unset = format(Sys.time(), "%Y%m%d_%H%M%S"))
OUT_DIR <- file.path(RESULTS_BASE, RUN_ID)
message("Run ID: ", RUN_ID, " -> ", OUT_DIR)

# As of 2026-07-31, inputs/dispersion_matrices and inputs/metadata are
# versioned per RUN_ID (must match whatever 01/02 were run with in this same
# pipeline invocation) -- see 01_prepare_dispersion_inputs.R's header and
# .living/decisions.md, 2026-07-31.
DISP_DIR <- file.path(DQTL_BASE, "inputs/dispersion_matrices", RUN_ID)
META_DIR <- file.path(DQTL_BASE, "inputs/metadata", RUN_ID)
stopifnot("inputs/dispersion_matrices/{RUN_ID}/ not found -- run 01_prepare_dispersion_inputs.R (and 02) with the same DQTL_RUN_ID first" =
            dir.exists(DISP_DIR))
stopifnot("inputs/metadata/{RUN_ID}/ not found -- run 01_prepare_dispersion_inputs.R (and 02) with the same DQTL_RUN_ID first" =
            dir.exists(META_DIR))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Env-var overridable to match 01's config (see DQTL_MAPPING.md "Run
# configuration"); defaults match 01's defaults.
CELL_TYPES <- strsplit(
  Sys.getenv("DQTL_CELL_TYPES",
             unset = "Activated_fibroblasts,Cardiac_fibroblasts,Cardiomyocytes"),
  ","
)[[1]]
# Env-var overridable (default unchanged) so a threshold sweep doesn't require
# editing this file -- see ENV_VARS.md. Recomputing at a *different* threshold
# after the fact doesn't need this override at all: use the raw p-value cache
# this script writes below, via recompute_covariate_sweep_at_threshold.R.
BH_FDR_THRESHOLD <- as.numeric(Sys.getenv("DQTL_COVARIATE_BH_THRESHOLD", unset = "0.1"))
CIS_DIST         <- 250e3

# Strip {cell_type} and {_B#} from sample IDs to get bare donor ID for genotype matching
# NOTE: duplicated in 04_map_dqtls.R — keep both copies in sync
strip_celltype_batch <- function(x, ct) {
  x <- sub(paste0("_", ct, "$"), "", x)
  sub("_B[0-9]+$", "", x)
}

# ---- Gene positions (shared across all cell types) ----
genepos_all <- read.csv(GENEPOS_FILE, stringsAsFactors = FALSE) |>
  dplyr::select(gene, chr, s1 = start, s2 = start)  # s2 = start: TSS-centred +/-250kb window (intentional)

for (ct in CELL_TYPES) {
  message("\n=== ", ct, " ===")

  disp_f <- file.path(DISP_DIR, paste0(ct, "_dispersion.tsv"))
  meta_f <- file.path(META_DIR, paste0(ct, "_metadata_sex_pcs.tsv"))
  stopifnot(file.exists(disp_f), file.exists(meta_f))

  # Snapshot the full Sex+PC1..N covariate matrix now, before it's overwritten
  # in place by a future run under a different config -- this is what lets
  # recompute_covariate_sweep_at_threshold.R rebuild a metadata_ready file for
  # any PC count later, without needing inputs/ restored again.
  file.copy(meta_f, file.path(OUT_DIR, paste0(ct, "_metadata_sex_pcs.tsv")), overwrite = TRUE)

  # ---- Load dispersion phenotype ----
  gene <- SlicedData$new()
  gene$fileDelimiter      <- "\t"
  gene$fileOmitCharacters <- "NA"
  gene$fileSkipRows       <- 1
  gene$fileSkipColumns    <- 1
  gene$fileSliceSize      <- 2000
  gene$LoadFile(disp_f)
  pheno_id <- gene$columnNames

  # ---- Load full covariate file (Sex + PC1..20) ----
  cvrt <- SlicedData$new()
  cvrt$fileDelimiter      <- "\t"
  cvrt$fileOmitCharacters <- "NA"
  cvrt$fileSkipRows       <- 1
  cvrt$fileSkipColumns    <- 1
  cvrt$LoadFile(meta_f)
  meta_id        <- cvrt$columnNames
  all_covariates <- cvrt$GetAllRowNames()
  stopifnot("Sex" %in% all_covariates)

  # ---- Gene positions for this cell type ----
  expr_genes <- gene$GetAllRowNames()
  genepos    <- dplyr::filter(genepos_all, gene %in% expr_genes) |>
    dplyr::group_by(gene, chr) |>
    dplyr::slice_min(s1, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  # ---- Build covariate configurations: Sex, Sex+PC1, ..., Sex+PC1..N ----
  # N (n_pcs_avail) is read back from how many PC rows this cell type's own
  # metadata file actually has, rather than hardcoded to 20 -- script 02 caps
  # PCs at min(20, n_samples - 1), so cell types with <=20 samples (e.g.
  # Mesenchymal_cells, Epithelial_cells) have fewer. See .living/decisions.md,
  # 2026-07-15.
  n_pcs_avail <- sum(grepl("^PC[0-9]+$", all_covariates))

  # Further cap by this cell type's own sample count: MatrixEQTL's covariate
  # step requires nrow(cvrt) < n_samples, and the actual per-gene cis test adds
  # intercept + SNP on top of cvrt. Cap PCs so the largest swept config (Sex +
  # PCs) still leaves >=1 residual df once those are added. Script 02's cap
  # (min(20, n_samples - 1)) only guarantees PCA validity, not regression
  # feasibility -- surfaced 2026-07-15 when Epithelial_cells (n=20) and
  # Mesenchymal_cells (n=19) crashed matrixEQTL with "number of covariates
  # exceeds number of samples" at the top of the sweep range. See
  # .living/decisions.md.
  n_samples   <- length(pheno_id)
  n_pcs_avail <- min(n_pcs_avail, max(n_samples - 4L, 0L))
  stopifnot("No PC covariates found in metadata file" = n_pcs_avail > 0)

  covariate_configs <- c(
    list(list(label = "Sex", covariates = "Sex")),
    lapply(seq_len(n_pcs_avail), function(k) {
      list(label = if (k == 1) "Sex and PC1" else paste0("Sex and PC1 to PC", k),
           covariates = c("Sex", paste0("PC", 1:k)))
    })
  )

  # ---- Sweep configurations (parallelized per chromosome) ----
  sig_counts     <- integer(length(covariate_configs))
  summary_lines  <- character(length(covariate_configs))
  raw_pvals_list <- vector("list", length(covariate_configs))

  for (cfg_i in seq_along(covariate_configs)) {
    covariate_label <- covariate_configs[[cfg_i]]$label
    covariate_names <- covariate_configs[[cfg_i]]$covariates

    keep_rows <- all_covariates %in% covariate_names

    cis_gene_pvals <- foreach(chr = 1:22,
                               .packages = c("MatrixEQTL", "Matrix"),
                               .combine  = rbind) %dopar% {

      snp_f <- file.path(GENO_DIR, paste0(chr, "_snp.txt"))
      snps  <- SlicedData$new()
      snps$fileDelimiter      <- "\t"
      snps$fileOmitCharacters <- "NA"
      snps$fileSkipRows       <- 1
      snps$fileSkipColumns    <- 1
      snps$fileSliceSize      <- 2000
      snps$LoadFile(snp_f)
      geno_id <- snps$columnNames

      pheno_base <- strip_celltype_batch(pheno_id, ct)

      # MAF filter (donor-level, > 5%)
      donor_base <- unique(pheno_base)
      donor_idx  <- match(donor_base, geno_id)
      if (anyNA(donor_idx)) {
        miss <- donor_base[is.na(donor_idx)]
        stop("Genotype columns missing for donors: ", paste(head(miss, 5), collapse = ", "))
      }
      snps_donor <- snps$Clone()
      snps_donor$ColumnSubsample(donor_idx)
      maf_mask <- unlist(lapply(seq_along(snps_donor), function(sl) {
        af <- rowMeans(snps_donor[[sl]], na.rm = TRUE) / 2
        pmin(af, 1 - af) > 0.05
      }))

      # Align covariates to phenotype samples
      idx_cvrt <- match(pheno_id, meta_id)
      if (anyNA(idx_cvrt)) stop("Covariate columns missing for some samples")
      cvrt_chr <- cvrt$Clone()
      cvrt_chr$ColumnSubsample(idx_cvrt)
      cvrt_chr$RowReorder(keep_rows)

      # Drop zero-variance covariates
      var_list <- lapply(seq_along(cvrt_chr), function(i) apply(cvrt_chr[[i]], 1, var))
      keep_cov <- unlist(var_list) > 0
      if (!all(keep_cov)) cvrt_chr$RowReorder(keep_cov)

      # Align genotypes to phenotype samples (replicate for duplicate individuals e.g. NA19114)
      idx_snp <- match(pheno_base, geno_id)
      snps_chr <- snps$Clone()
      snps_chr$ColumnSubsample(idx_snp)
      snps_chr$RowReorder(maf_mask)

      gene_chr <- gene$Clone()
      snps_chr$columnNames <- pheno_id
      gene_chr$columnNames <- pheno_id
      cvrt_chr$columnNames <- pheno_id

      snppos_f    <- file.path(SNPPOS_DIR, paste0(chr, "_snppos.txt"))
      snpspos     <- read.table(snppos_f, header = TRUE, stringsAsFactors = FALSE)
      snpspos_idx <- match(snps_chr$GetAllRowNames(), snpspos$snp)
      stopifnot("SNP position lookup failed for some cis-tested SNPs" = !anyNA(snpspos_idx))
      # ANALYSIS_OK[self-assignment]: reorder only, not a narrowing filter --
      # the stopifnot immediately above guarantees snpspos_idx has no NA, so
      # every row of snps_chr has a matching position row here.
      snpspos <- snpspos[snpspos_idx, ]

      me <- Matrix_eQTL_main(
        snps = snps_chr, gene = gene_chr, cvrt = cvrt_chr,
        useModel = modelLINEAR, verbose = FALSE,
        output_file_name.cis  = NULL,
        pvOutputThreshold.cis = 1,
        snpspos = snpspos, genepos = genepos, cisDist = CIS_DIST,
        pvalue.hist = FALSE, min.pv.by.genesnp = FALSE, noFDRsaveMemory = FALSE
      )

      eqtls <- me$cis$eqtls
      if (is.null(eqtls) || nrow(eqtls) == 0) {
        data.frame(gene = character(0), pvalue = numeric(0))
      } else {
        pcol <- grep("p.?value", colnames(eqtls), ignore.case = TRUE, value = TRUE)[1]
        gcol <- grep("^gene$", colnames(eqtls), ignore.case = TRUE, value = TRUE)[1]
        data.frame(gene = eqtls[[gcol]], pvalue = eqtls[[pcol]])
      }
    }

    if (is.null(cis_gene_pvals) || nrow(cis_gene_pvals) == 0) {
      n_sig <- 0L
      raw_pvals_list[[cfg_i]] <- data.frame(label = character(0), n_pcs = integer(0),
                                             gene = character(0), pvalue = numeric(0))
    } else {
      bh_fdr    <- p.adjust(cis_gene_pvals$pvalue, method = "BH")
      sig_genes <- unique(cis_gene_pvals$gene[bh_fdr <= BH_FDR_THRESHOLD & !is.na(bh_fdr)])
      n_sig     <- length(sig_genes)
      # Raw per-pair p-values (not gene-deduped, not BH-adjusted) -- this is
      # the same quantity p.adjust() above is computed from, so any future BH
      # threshold can be applied to it directly and reproduce this exact
      # counting method (pair-level BH, then dedupe by gene).
      raw_pvals_list[[cfg_i]] <- data.frame(
        label  = covariate_label,
        n_pcs  = length(covariate_names) - 1L,
        gene   = cis_gene_pvals$gene,
        pvalue = cis_gene_pvals$pvalue
      )
    }

    sig_counts[cfg_i]    <- n_sig
    summary_lines[cfg_i] <- paste0(covariate_label, ":\t", n_sig)
    message("  ", covariate_label, ": ", n_sig, " significant dGenes")
  }

  # Write sweep summary
  writeLines(summary_lines,
             file.path(OUT_DIR, paste0(ct, "_covariate_sweep.txt")))

  raw_pvals_f  <- file.path(OUT_DIR, paste0(ct, "_covariate_sweep_raw_pvals.rds"))
  raw_pvals_df <- dplyr::bind_rows(raw_pvals_list)
  saveRDS(raw_pvals_df, raw_pvals_f)
  message("  Written: ", raw_pvals_f, " (", nrow(raw_pvals_df), " gene-config p-value rows)")

  # ---- Auto-select optimal PC count (max significant dGenes) ----
  optimal_i <- which.max(sig_counts)
  optimal_k <- length(covariate_configs[[optimal_i]]$covariates) - 1L  # subtract Sex
  message("  Optimal: ", covariate_configs[[optimal_i]]$label,
          " (", sig_counts[optimal_i], " dGenes)")

  # ---- Write optimized metadata file ----
  meta_df        <- read.table(meta_f, header = TRUE, sep = "\t", check.names = FALSE)
  rows_keep      <- c("Sex", if (optimal_k > 0) paste0("PC", seq_len(optimal_k)))
  meta_ready     <- meta_df[meta_df$id %in% rows_keep, , drop = FALSE]
  meta_ready_idx <- match(rows_keep, meta_ready$id)
  stopifnot("Optimal covariate rows missing from metadata file" = !anyNA(meta_ready_idx))
  # ANALYSIS_OK[self-assignment]: reorder only, not a narrowing filter -- the
  # stopifnot immediately above guarantees meta_ready_idx has no NA, so every
  # name in rows_keep has a matching row here.
  meta_ready <- meta_ready[meta_ready_idx, ]

  # Archived auto-pick suggestion only -- NOT the live inputs/metadata/ file
  # 04_map_dqtls.R reads. 03b_apply_manual_pcs.R overwrites this same path
  # with the manually-chosen PC count once you've reviewed the sweep/plot.
  suggested_out <- file.path(OUT_DIR, paste0(ct, "_metadata_ready.tsv"))
  write.table(meta_ready, suggested_out, sep = "\t", quote = FALSE, row.names = FALSE)
  message("  Auto-pick suggestion written: ", suggested_out,
          " (review the sweep + PDF and run 03b_apply_manual_pcs.R before mapping)")
}

# ---- Record run provenance (config this sweep was run under) ----
# ANALYSIS_OK[bare-literal-fallback]: "unknown" is an intentional placeholder
# for a non-critical provenance field (run_info.txt only) -- git
# unavailability shouldn't block this sweep, and the message() logs the
# fallback explicitly rather than substituting it silently.
git_commit <- tryCatch(
  system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE),
  error = function(e) {
    message("Could not determine git commit (", conditionMessage(e), "); using 'unknown'")
    "unknown"
  }
)
if (length(git_commit) == 0) git_commit <- "unknown"

prep_config_f <- file.path(DQTL_BASE, "inputs/gene_sets", RUN_ID, "prep_config.txt")
prep_config   <- if (file.exists(prep_config_f)) {
  prep_lines <- readLines(prep_config_f)
  setNames(trimws(sub("^[^:]+:", "", prep_lines)), trimws(sub(":.*$", "", prep_lines)))
} else {
  message("inputs/gene_sets/", RUN_ID, "/prep_config.txt not found -- leaving cell_count_threshold/",
          "sd_cutoff_frac/n_genes_selected blank in run_info.txt")
  c(cell_count_threshold = "", sd_cutoff_frac = "", n_genes_selected = "")
}

writeLines(
  c(paste0("run_id: ", RUN_ID),
    paste0("created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste0("git_commit: ", git_commit),
    paste0("cell_count_threshold: ", prep_config[["cell_count_threshold"]]),
    paste0("sd_cutoff_frac: ", prep_config[["sd_cutoff_frac"]]),
    paste0("n_cell_types: ", length(CELL_TYPES)),
    paste0("cell_types: ", paste(CELL_TYPES, collapse = ",")),
    paste0("n_genes_selected: ", prep_config[["n_genes_selected"]]),
    paste0("bh_fdr_threshold: ", BH_FDR_THRESHOLD)),
  file.path(OUT_DIR, "run_info.txt")
)

message("\nDone. Sweep results: ", OUT_DIR)
