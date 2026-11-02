stopifnot(file.exists("CLAUDE.md"))

library(tidyverse)
library(data.table)

DQTL_BASE    <- "analysis/dqtl_mapping"
RESULTS_BASE <- file.path(DQTL_BASE, "outputs/dqtl_results")
REFORM_BASE  <- file.path(DQTL_BASE, "outputs/reformatted_dqtls")

args <- commandArgs(trailingOnly = TRUE)
RUN_ID <- if (length(args) >= 1) args[1] else Sys.getenv("DQTL_RUN_ID", unset = "")
stopifnot("No RUN_ID given -- pass one as the first CLI arg or set DQTL_RUN_ID" =
            !identical(RUN_ID, ""))
message("Run ID: ", RUN_ID)

DQTL_DIR   <- file.path(RESULTS_BASE, RUN_ID)
REFORM_DIR <- file.path(REFORM_BASE, RUN_ID)
LIST_DIR   <- file.path(REFORM_DIR, "list_files")
stopifnot("matrixEQTL results not found for this RUN_ID" = dir.exists(DQTL_DIR))

dir.create(REFORM_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LIST_DIR,   recursive = TRUE, showWarnings = FALSE)

# CELL_TYPES for this run isn't hardcoded here -- it's read back from
# run_info.txt (written by 04_map_dqtls.R) so this script always reformats
# whatever cell types that RUN_ID actually mapped, even if DQTL_CELL_TYPES
# was swept to a non-default set (see DQTL_MAPPING.md "Run configuration").
run_info_f <- file.path(DQTL_DIR, "run_info.txt")
stopifnot("run_info.txt missing for this RUN_ID -- re-run 04_map_dqtls.R" =
            file.exists(run_info_f))
run_info_lines <- readLines(run_info_f)
run_info <- setNames(trimws(sub("^[^:]+:", "", run_info_lines)),
                     trimws(sub(":.*$", "", run_info_lines)))
CELL_TYPES <- strsplit(run_info[["cell_types"]], ",")[[1]]

combined_list <- list()

# ---- Step 1-2: Combine per-chromosome files and BH-adjust ----
for (ct in CELL_TYPES) {
  message("=== ", ct, " ===")
  ct_dir <- file.path(DQTL_DIR, ct)

  chr_files <- list.files(ct_dir,
                           pattern = paste0("^dQTL_", ct, "_chr(\\d{1,2})\\.txt$"),
                           full.names = TRUE)

  missing_chrs <- setdiff(as.character(1:22),
                           as.character(as.integer(sub(".*_chr(\\d{1,2})\\.txt$", "\\1",
                                                        basename(chr_files)))))
  if (length(missing_chrs) > 0) {
    stop(ct, ": missing chromosome files: ", paste(missing_chrs, collapse = ", "))
  }

  all_eqtls <- bind_rows(lapply(chr_files, read_tsv, col_types = cols()))
  all_eqtls <- rename(all_eqtls, pvalue = `p-value`, t_stat = `t-stat`)
  all_eqtls$bh_adj_pvalue <- p.adjust(all_eqtls$pvalue, method = "BH")

  # Write AllChroms file alongside the per-chr files
  allchroms_out <- file.path(ct_dir, paste0("dQTL_", ct, "_AllChroms.txt"))
  write_tsv(all_eqtls, allchroms_out)
  message("  AllChroms: ", nrow(all_eqtls), " rows, ",
          sum(all_eqtls$bh_adj_pvalue <= 0.05, na.rm = TRUE), " sig at FDR 0.05")

  combined_list[[ct]] <- all_eqtls
}

# ---- Step 3: Reformat to fastQTL format, with placeholder-filling ----
#
# Gene selection is independent per cell type (01_prepare_dispersion_inputs.R
# no longer unions gene sets across cell types -- see .living/decisions.md,
# 2026-07-15), so most (gene, snps) pairs are only ever tested in one cell
# type. To let mash still borrow cross-condition information for those genes,
# missing (gene, snps, cell_type) entries are filled with placeholder values
# rather than dropped -- reintroduces the pre-2026-06-30 placeholder strategy
# (see .living/decisions.md, 2026-06-23), but with SE derived from (beta,
# pvalue) via mashr's own p2z() conversion instead of a fixed 1e16 for every
# placeholder row:
#   1. beta:   NA -> 0
#   2. pvalue: NA -> 1
#   3. se := beta / mashr:::p2z(pvalue, beta)  (for real rows: the SE implied
#      by the p-value's equivalent z-score, consistent with mash's own
#      z = beta/se assumption, rather than matrixEQTL's raw t-statistic SE;
#      for placeholder rows -- beta=0, pvalue=1 -- z = p2z(1, 0) = 0, so
#      0/0 = NaN, caught by step 4)
#   4. se: NA/NaN -> 1e16
stopifnot("mashr package required for placeholder SE computation (mashr:::p2z)" =
            requireNamespace("mashr", quietly = TRUE))

per_ct_long <- lapply(combined_list, function(df) {
  dt <- data.table::data.table(
    gene   = df$gene,
    snps   = df$SNP,
    beta   = ifelse(is.finite(df$beta),   df$beta,   NA_real_),
    pvalue = ifelse(is.finite(df$pvalue), df$pvalue, NA_real_)
  )
  dt[is.finite(beta) & is.finite(pvalue)]
})

long_all  <- data.table::rbindlist(per_ct_long, idcol = "cell_type")
beta_wide <- data.table::dcast(long_all, gene + snps ~ cell_type, value.var = "beta")
pval_wide <- data.table::dcast(long_all, gene + snps ~ cell_type, value.var = "pvalue")
stopifnot("dcast did not produce a column for every cell type" =
            all(CELL_TYPES %in% names(beta_wide)))

n_tested <- rowSums(!is.na(as.matrix(beta_wide[, ..CELL_TYPES])))
message("Gene-SNP pairs: ", nrow(beta_wide), " total (union across ", length(CELL_TYPES),
        " cell types); ", sum(n_tested == length(CELL_TYPES)),
        " tested in all cell types; placeholder-filled for the rest")

beta_mat <- as.matrix(beta_wide[, ..CELL_TYPES])
beta_mat[is.na(beta_mat)] <- 0
pval_mat <- as.matrix(pval_wide[, ..CELL_TYPES])
pval_mat[is.na(pval_mat)] <- 1
se_mat   <- beta_mat / mashr:::p2z(pval_mat, beta_mat)
se_mat[is.na(se_mat)] <- 1e16

reformatted_list <- setNames(lapply(CELL_TYPES, function(ct) {
  data.table::data.table(
    gene = beta_wide$gene, snps = beta_wide$snps,
    beta = beta_mat[, ct], se = se_mat[, ct], pvalue = pval_mat[, ct]
  )[order(gene)]
}), CELL_TYPES)

# ---- Step 4: Write .out.gz files and list file ----
out_files <- character(length(CELL_TYPES))

for (i in seq_along(CELL_TYPES)) {
  ct       <- CELL_TYPES[i]
  out_file <- file.path(REFORM_DIR, paste0("dQTL_", ct, "_AllChroms.out.gz"))
  data.table::fwrite(reformatted_list[[ct]], out_file,
                     sep = "\t", quote = FALSE, compress = "gzip")
  out_files[i] <- out_file
  message("Written: ", out_file)
}

# List file for fastqtl2mash (relative paths from list_files/ directory)
relative_paths <- file.path("..", basename(out_files))
writeLines(relative_paths, file.path(LIST_DIR, "all_celltypes.list"))

# ---- Record run provenance ----
run_info_src <- file.path(DQTL_DIR, "run_info.txt")
if (file.exists(run_info_src)) {
  file.copy(run_info_src, file.path(REFORM_DIR, "run_info.txt"), overwrite = TRUE)
}

message("\nList file written. Run fastqtl2mash next (see run.sh).")
message("Run ID: ", RUN_ID)
