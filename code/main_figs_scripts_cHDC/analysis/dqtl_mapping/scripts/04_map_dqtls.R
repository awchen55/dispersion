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
RESULTS_BASE <- file.path(DQTL_BASE, "outputs/dqtl_results")
COVOPT_DIR   <- file.path(DQTL_BASE, "outputs/covariate_optimization")

RUN_ID   <- Sys.getenv("DQTL_RUN_ID", unset = format(Sys.time(), "%Y%m%d_%H%M%S"))
OUT_BASE <- file.path(RESULTS_BASE, RUN_ID)
message("Run ID: ", RUN_ID, " -> ", OUT_BASE)

# As of 2026-07-31, inputs/dispersion_matrices is versioned per RUN_ID (must
# match whatever 01_prepare_dispersion_inputs.R was run with) -- see that
# script's header and .living/decisions.md, 2026-07-31.
DISP_DIR <- file.path(DQTL_BASE, "inputs/dispersion_matrices", RUN_ID)
stopifnot("inputs/dispersion_matrices/{RUN_ID}/ not found -- run 01_prepare_dispersion_inputs.R (and 02, 03, 03b) with the same DQTL_RUN_ID first" =
            dir.exists(DISP_DIR))

# CELL_TYPES/MAF_THRESHOLD are env-var overridable to match 01's config (see
# DQTL_MAPPING.md "Run configuration"); defaults match 01's defaults.
CELL_TYPES <- strsplit(
  Sys.getenv("DQTL_CELL_TYPES",
             unset = "Activated_fibroblasts,Cardiac_fibroblasts,Cardiomyocytes"),
  ","
)[[1]]
CIS_WINDOW_KB <- as.numeric(Sys.getenv("DQTL_CIS_WINDOW_KB", unset = "250"))
CIS_DIST      <- CIS_WINDOW_KB * 1e3
MAF_THRESHOLD <- as.numeric(Sys.getenv("DQTL_MAF_THRESHOLD", unset = "0.05"))

# NOTE: duplicated in 03_optimize_covariates.R — keep both copies in sync
strip_celltype_batch <- function(x, ct) {
  x <- sub(paste0("_", ct, "$"), "", x)
  sub("_B[0-9]+$", "", x)
}

genepos_all <- read.csv(GENEPOS_FILE, stringsAsFactors = FALSE) |>
  dplyr::select(gene, chr, s1 = start, s2 = start)  # s2 = start: TSS-centred +/-CIS_WINDOW_KB window (intentional)

# ---- Manual PC-selection provenance check (see header comment above) ----
# NOTE: manual_pc_selection.txt parsing duplicated in 03b_apply_manual_pcs.R
# -- keep both copies in sync.
manual_f <- file.path(COVOPT_DIR, RUN_ID, "manual_pc_selection.txt")
stopifnot(
  "manual_pc_selection.txt missing for this RUN_ID -- review outputs/covariate_optimization/{RUN_ID}/covariate_sweep_by_celltype.pdf and run 03b_apply_manual_pcs.R before mapping dQTLs" =
    file.exists(manual_f)
)
manual_lines <- readLines(manual_f)
# ANALYSIS_OK[narrowing-filter]: drops only blank lines from a hand-edited
# file (formatting convenience); any remaining non-blank line that isn't
# "{cell_type}\t{k}" still fails loudly at the length==2 stopifnot() directly
# below, so no malformed row is silently dropped.
manual_lines <- manual_lines[nzchar(trimws(manual_lines))]
manual_parts <- strsplit(manual_lines, "\t")
# ANALYSIS_OK[magic-number]: "2" is the fixed field count of the
# "{cell_type}\t{k}" line format (not a tunable/scientific threshold).
stopifnot("Every manual_pc_selection.txt line must be '{cell_type}\\t{k}'" =
            all(vapply(manual_parts, length, integer(1)) == 2))
manual_pcs <- setNames(
  as.integer(vapply(manual_parts, `[`, character(1), 2)),
  vapply(manual_parts, `[`, character(1), 1)
)
stopifnot("manual_pc_selection.txt is missing an entry for a cell type in DQTL_CELL_TYPES" =
            all(CELL_TYPES %in% names(manual_pcs)))

for (ct in CELL_TYPES) {
  message("\n=== ", ct, " ===")

  out_dir <- file.path(OUT_BASE, ct)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  disp_f <- file.path(DISP_DIR, paste0(ct, "_dispersion.tsv"))
  # Reads the RUN_ID-versioned snapshot (not inputs/metadata/, a shared path
  # that can already exist from an unrelated past RUN_ID) so a fresh RUN_ID
  # can never silently pick up a stale/wrong PC selection -- see
  # 03b_apply_manual_pcs.R and .living/decisions.md, 2026-07-31.
  meta_f <- file.path(COVOPT_DIR, RUN_ID, paste0(ct, "_metadata_ready.tsv"))
  stopifnot(file.exists(disp_f), file.exists(meta_f))

  # Confirm this snapshot's actual PC count matches the manual choice for ct
  # -- catches a stale snapshot (e.g. an old auto-pick, or a manual choice
  # from a previous review) that Snakemake's DAG alone would not detect.
  ready_ids  <- read.table(meta_f, header = TRUE, sep = "\t", check.names = FALSE)$id
  actual_k   <- sum(grepl("^PC[0-9]+$", ready_ids))
  if (actual_k != manual_pcs[[ct]]) {
    stop(ct, ": ", meta_f, " has ", actual_k, " PC(s), but manual_pc_selection.txt specifies ",
         manual_pcs[[ct]], " -- run 03b_apply_manual_pcs.R to regenerate it before mapping dQTLs")
  }

  # ---- Load dispersion phenotype ----
  gene <- SlicedData$new()
  gene$fileDelimiter      <- "\t"
  gene$fileOmitCharacters <- "NA"
  gene$fileSkipRows       <- 1
  gene$fileSkipColumns    <- 1
  gene$fileSliceSize      <- 2000
  gene$LoadFile(disp_f)
  pheno_id <- gene$columnNames

  # ---- Load optimized covariates ----
  cvrt <- SlicedData$new()
  cvrt$fileDelimiter      <- "\t"
  cvrt$fileOmitCharacters <- "NA"
  cvrt$fileSkipRows       <- 1
  cvrt$fileSkipColumns    <- 1
  cvrt$LoadFile(meta_f)
  meta_id        <- cvrt$columnNames
  all_covariates <- cvrt$GetAllRowNames()
  stopifnot("Sex" %in% all_covariates)
  message("  Covariates: ", paste(all_covariates, collapse = ", "))

  # ---- Gene positions for this cell type ----
  expr_genes <- gene$GetAllRowNames()
  genepos    <- dplyr::filter(genepos_all, gene %in% expr_genes) |>
    dplyr::group_by(gene, chr) |>
    dplyr::slice_min(s1, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  pheno_base <- strip_celltype_batch(pheno_id, ct)

  # ---- Map per chromosome (parallelized) ----
  foreach(chr = 1:22, .packages = c("MatrixEQTL", "Matrix")) %dopar% {

    snp_f <- file.path(GENO_DIR, paste0(chr, "_snp.txt"))
    snps  <- SlicedData$new()
    snps$fileDelimiter      <- "\t"
    snps$fileOmitCharacters <- "NA"
    snps$fileSkipRows       <- 1
    snps$fileSkipColumns    <- 1
    snps$fileSliceSize      <- 2000
    snps$LoadFile(snp_f)
    geno_id <- snps$columnNames

    # MAF filter (donor-level)
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
      pmin(af, 1 - af) > MAF_THRESHOLD
    }))

    # Align covariates to phenotype samples
    idx_cvrt <- match(pheno_id, meta_id)
    if (anyNA(idx_cvrt)) stop(ct, " chr", chr, ": covariate columns missing for some samples")
    cvrt_chr <- cvrt$Clone()
    cvrt_chr$ColumnSubsample(idx_cvrt)

    # Drop zero-variance covariates
    var_list <- lapply(seq_along(cvrt_chr), function(i) apply(cvrt_chr[[i]], 1, var))
    keep_cov <- unlist(var_list) > 0
    if (!all(keep_cov)) cvrt_chr$RowReorder(keep_cov)

    # Align genotypes to phenotype samples
    idx_snp  <- match(pheno_base, geno_id)
    snps_chr <- snps$Clone()
    snps_chr$ColumnSubsample(idx_snp)
    snps_chr$RowReorder(maf_mask)

    gene_chr <- gene$Clone()
    snps_chr$columnNames <- pheno_id
    gene_chr$columnNames <- pheno_id
    cvrt_chr$columnNames <- pheno_id

    stopifnot(identical(gene_chr$columnNames, snps_chr$columnNames),
              identical(gene_chr$columnNames, cvrt_chr$columnNames))

    snppos_f    <- file.path(SNPPOS_DIR, paste0(chr, "_snppos.txt"))
    snpspos     <- read.table(snppos_f, header = TRUE, stringsAsFactors = FALSE)
    snpspos_idx <- match(snps_chr$GetAllRowNames(), snpspos$snp)
    stopifnot("SNP position lookup failed for some cis-tested SNPs" = !anyNA(snpspos_idx))
    # ANALYSIS_OK[self-assignment]: reorder only, not a narrowing filter -- the
    # stopifnot immediately above guarantees snpspos_idx has no NA, so every
    # row of snps_chr has a matching position row here.
    snpspos <- snpspos[snpspos_idx, ]

    out_file <- file.path(out_dir, paste0("dQTL_", ct, "_chr", chr, ".txt"))

    Matrix_eQTL_main(
      snps = snps_chr, gene = gene_chr, cvrt = cvrt_chr,
      useModel = modelLINEAR, verbose = TRUE,
      output_file_name.cis  = out_file,
      pvOutputThreshold.cis = 1,
      snpspos = snpspos, genepos = genepos, cisDist = CIS_DIST,
      pvalue.hist = FALSE, min.pv.by.genesnp = FALSE, noFDRsaveMemory = FALSE
    )

    message("  chr", chr, " done -> ", out_file)
  }

  message(ct, ": all chromosomes complete")
}

# ---- Snapshot gene-set provenance for this run (co-locates it with this
# RUN_ID's mapping output; inputs/gene_sets/{RUN_ID}/ is itself already a
# permanent, non-overwritten record as of 2026-07-31 -- see
# 01_prepare_dispersion_inputs.R's header) ----
gene_sets_src <- file.path(DQTL_BASE, "inputs/gene_sets", RUN_ID)
stopifnot("inputs/gene_sets/{RUN_ID}/ not found -- run 01_prepare_dispersion_inputs.R with the same DQTL_RUN_ID first" =
            dir.exists(gene_sets_src))
gene_sets_dst <- file.path(OUT_BASE, "gene_sets")
dir.create(gene_sets_dst, recursive = TRUE, showWarnings = FALSE)
file.copy(list.files(gene_sets_src, full.names = TRUE), gene_sets_dst, overwrite = TRUE)

prep_config_f <- file.path(gene_sets_dst, "prep_config.txt")
stopifnot("inputs/gene_sets/{RUN_ID}/prep_config.txt missing -- re-run 01_prepare_dispersion_inputs.R" =
            file.exists(prep_config_f))
prep_lines  <- readLines(prep_config_f)
prep_config <- setNames(trimws(sub("^[^:]+:", "", prep_lines)),
                        trimws(sub(":.*$", "", prep_lines)))
stopifnot("CELL_TYPES here doesn't match cell_types recorded in inputs/gene_sets/{RUN_ID}/prep_config.txt -- re-run 01_prepare_dispersion_inputs.R with a matching DQTL_CELL_TYPES" =
            identical(prep_config[["cell_types"]], paste(CELL_TYPES, collapse = ",")))

# ---- Record run provenance ----
# ANALYSIS_OK[bare-literal-fallback]: "unknown" is an intentional placeholder
# for a non-critical provenance field (run_info.txt/registry only) -- git
# unavailability shouldn't block a matrixEQTL run, and the message() above
# logs the fallback explicitly rather than substituting it silently.
git_commit <- tryCatch(
  system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE),
  error = function(e) {
    message("Could not determine git commit (", conditionMessage(e), "); using 'unknown'")
    "unknown"
  }
)
if (length(git_commit) == 0) git_commit <- "unknown"
writeLines(
  c(paste0("run_id: ", RUN_ID),
    paste0("created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste0("git_commit: ", git_commit),
    paste0("cell_count_threshold: ", prep_config[["cell_count_threshold"]]),
    paste0("sd_cutoff_frac: ", prep_config[["sd_cutoff_frac"]]),
    paste0("n_cell_types: ", length(CELL_TYPES)),
    paste0("cell_types: ", paste(CELL_TYPES, collapse = ",")),
    paste0("n_genes_selected: ", prep_config[["n_genes_selected"]]),
    paste0("maf_threshold: ", MAF_THRESHOLD),
    paste0("cis_window_kb: ", CIS_WINDOW_KB)),
  file.path(OUT_BASE, "run_info.txt")
)

# ---- Append to cross-run registry for easy config comparison across runs ----
registry_f    <- file.path(RESULTS_BASE, "run_registry.csv")
registry_row  <- data.frame(
  run_id               = RUN_ID,
  created              = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
  git_commit           = git_commit,
  cell_count_threshold = prep_config[["cell_count_threshold"]],
  sd_cutoff_frac       = prep_config[["sd_cutoff_frac"]],
  n_cell_types         = length(CELL_TYPES),
  cell_types           = paste(CELL_TYPES, collapse = ";"),
  n_genes_selected     = prep_config[["n_genes_selected"]],
  maf_threshold        = MAF_THRESHOLD,
  cis_window_kb        = CIS_WINDOW_KB,
  stringsAsFactors = FALSE
)
write.table(registry_row, registry_f, sep = ",", row.names = FALSE,
            col.names = !file.exists(registry_f), append = file.exists(registry_f))
message("Appended run config to ", registry_f)

message("\nDone. Run ID: ", RUN_ID)
