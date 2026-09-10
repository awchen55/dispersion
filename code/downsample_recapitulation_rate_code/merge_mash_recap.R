#!/usr/bin/env Rscript
#
# For one downsampled allotetraploid memento DD run (one cell type, one
# cell-count / replicate), rebuild the full 21-cell-type tetraploid
# logFC/SE/pval matrices (production CSVs for every cell type except the
# downsampled target), refit MASH on them exactly as
# de_dd_snakemake_mt_rb_mean_filtered/scripts/tetraploid_dd_mash.R does, and
# compute the recapitulation rate for the target cell type against the
# (unchanged) production diploid DD results.
#
# Diploid data is never touched here -- only the target cell type's
# allotetraploid arm is downsampled -- so recapitulation rate as a function
# of allotetraploid cell count can be read directly off the diploid-DD gene
# set, held fixed across every point on the curve.

suppressPackageStartupMessages({
  library(mashr)
  library(optparse)
})

opt_list <- list(
  make_option("--target-ct", type = "character"),
  make_option("--target-csv", type = "character"),
  make_option("--n-cells", type = "integer"),
  make_option("--replicate", type = "integer"),
  make_option("--out-csv", type = "character"),
  make_option("--prod-dir", type = "character",
              default = "/project/gilad/awchen55/mechanism/data/dd_de_results/mt_rb_mean_filtered_q_thresholds/q0.2"),
  make_option("--pca-components", type = "integer", default = 5),
  make_option("--lfsr-threshold", type = "double", default = 0.05),
  make_option("--save-mash", action = "store_true", default = FALSE)
)
opt <- parse_args(OptionParser(option_list = opt_list))

celltypes_tetraploid <- c(
  "Bronchiolar_and_alveolar_epithelial_cells", "Cardiac_mesoderm_cells",
  "Cardiomyocytes", "Ciliated_epithelial_cells", "ENS_glia",
  "Early_mesoderm_cells", "Epicardial_fat_cells", "Glioblast", "Goblet_cells",
  "Hepatoblasts", "Intestinal_epithelial_cells", "Metanephric_cells",
  "Neuroblast", "Neuron", "Neuronal_IPC", "Radial_glia", "Retinal_neurons",
  "Squamous_epithelial_cells", "Stromal_cells", "Vascular_endothelial_cells",
  "iPSCs"
)

target_ct2 <- gsub("_", ".", opt[["target-ct"]])

# ---------------------------------------------------------------------------
# 1. Rebuild merged logFC / SE / pval matrices (tetraploid_dd_merge.R logic),
#    substituting the target cell type's production CSV with the downsampled
#    run's CSV.
# ---------------------------------------------------------------------------
merge_full <- function(existing, new_df) {
  if (is.null(existing)) return(new_df)
  merged <- merge(existing, new_df, by = "row.names", all = TRUE)
  rownames(merged) <- merged$Row.names
  merged$Row.names <- NULL
  merged
}

logfc_df <- se_df <- pval_df <- NULL

for (ct in celltypes_tetraploid) {
  ct2 <- gsub("_", ".", ct)
  csv_path <- if (ct == opt[["target-ct"]]) {
    opt[["target-csv"]]
  } else {
    file.path(opt[["prod-dir"]], "memento",
              paste0("tetraploid_memento_dd_results_", ct, ".csv"))
  }
  res <- read.csv(csv_path)

  make_df <- function(vals, prefix) {
    df <- data.frame(vals, row.names = res$gene)
    colnames(df) <- paste0(prefix, ct2)
    df
  }

  logfc_df <- merge_full(logfc_df, make_df(res$dv_coef, "dd.logFC.allotetraploid."))
  se_df    <- merge_full(se_df,    make_df(res$dv_se,   "dd.SE.allotetraploid."))
  pval_df  <- merge_full(pval_df,  make_df(res$dv_pval, "dd.pvalue.allotetraploid."))
}

# ---------------------------------------------------------------------------
# 2. Refit MASH on the rebuilt matrices (tetraploid_dd_mash.R logic).
# ---------------------------------------------------------------------------
gene_list <- read.csv(file.path(opt[["prod-dir"]], "dd",
                                 "diploid_de_dd_results_filtered.csv"))$gene

betahat <- as.matrix(logfc_df)
pval    <- as.matrix(pval_df)

keep    <- intersect(gene_list, rownames(betahat))
betahat <- betahat[keep, , drop = FALSE]
pval    <- pval[keep, , drop = FALSE]

betahat[is.na(betahat)] <- 0
pval[is.na(pval)]       <- 1

pval[pval == 0] <- 1e-300
shat <- betahat / mashr:::p2z(pval, betahat)
shat[is.na(shat) | !is.finite(shat)] <- 1e16

data_mash <- mash_set_data(betahat, shat)
m.1by1    <- mash_1by1(data_mash)
strong    <- get_significant_results(m.1by1, opt[["lfsr-threshold"]])

U.pca <- cov_pca(data_mash, opt[["pca-components"]], subset = strong)
U.ed  <- cov_ed(data_mash, U.pca, subset = strong)
m.ed  <- mash(data_mash, Ulist = c(U.ed, U.pca))

if (opt[["save-mash"]]) {
  saveRDS(m.ed, sub("\\.csv$", "_mash.rds", opt[["out-csv"]]))
}

# ---------------------------------------------------------------------------
# 3. Recapitulation rate for the target cell type vs. the fixed diploid
#    DD results (production; unchanged by this downsampling exercise).
# ---------------------------------------------------------------------------
allo_lfsr_col  <- paste0("dd.logFC.allotetraploid.", target_ct2)  # get_lfsr keeps betahat colnames
allo_logfc_col <- paste0("dd.logFC.allotetraploid.", target_ct2)

lfsr_allo  <- get_lfsr(m.ed)[, allo_lfsr_col]
logfc_allo <- get_pm(m.ed)[, allo_logfc_col]
allo_genes <- rownames(get_lfsr(m.ed))

dip_logfc_all <- read.csv(
  file.path(opt[["prod-dir"]], "dd", "mash",
            "diploid_memento_limma_all_celltypes_filtered_logFC.csv"),
  row.names = 1, check.names = FALSE
)
dip_mash <- readRDS(
  file.path(opt[["prod-dir"]], "dd", "mash",
            "diploid_memento_limma_pval_filtered_mash_results.rds")
)
dip_lfsr_all <- get_lfsr(dip_mash)
colnames(dip_lfsr_all) <- gsub("dd\\.logFC\\.diploid\\.", "dd.lfsr.diploid.",
                                colnames(dip_lfsr_all))

dip_logfc_col <- paste0("dd.logFC.diploid.", target_ct2)
dip_lfsr_col  <- paste0("dd.lfsr.diploid.", target_ct2)

common_genes <- Reduce(intersect, list(
  allo_genes, rownames(dip_logfc_all), rownames(dip_lfsr_all)
))

dip_lfsr  <- dip_lfsr_all[common_genes, dip_lfsr_col]
dip_logfc <- dip_logfc_all[common_genes, dip_logfc_col]
lfsr_allo  <- lfsr_allo[common_genes]
logfc_allo <- logfc_allo[common_genes]

dd_flag_dip <- !is.na(dip_lfsr) & dip_lfsr < opt[["lfsr-threshold"]]
n_dd_diploid <- sum(dd_flag_dip, na.rm = TRUE)

recap_flag <- dd_flag_dip &
  !is.na(lfsr_allo) & lfsr_allo < opt[["lfsr-threshold"]] &
  sign(dip_logfc) == sign(logfc_allo)
n_recap <- sum(recap_flag, na.rm = TRUE)

out <- data.frame(
  celltype     = opt[["target-ct"]],
  n_cells      = opt[["n-cells"]],
  replicate    = opt[["replicate"]],
  n_dd_diploid = n_dd_diploid,
  n_recap      = n_recap,
  pct_recap    = n_recap / n_dd_diploid
)

write.csv(out, opt[["out-csv"]], row.names = FALSE)
message(sprintf(
  "[%s] n_cells=%d rep=%d: %d/%d recapitulated (%.1f%%)",
  opt[["target-ct"]], opt[["n-cells"]], opt[["replicate"]],
  n_recap, n_dd_diploid, 100 * n_recap / n_dd_diploid
))
