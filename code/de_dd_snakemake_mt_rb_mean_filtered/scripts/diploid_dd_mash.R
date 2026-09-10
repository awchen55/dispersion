library(mashr)

logfc_path <- snakemake@input[["filtered_logfc"]]
se_path    <- snakemake@input[["filtered_se"]]
pval_path  <- snakemake@input[["filtered_pval"]]
out_rds    <- snakemake@output[["mash_rds"]]
n_pca      <- snakemake@params[["pca_components"]]
lfsr_thr   <- snakemake@params[["lfsr_threshold"]]

read_mat <- function(path) {
  df <- read.csv(path, row.names = 1, check.names = FALSE)
  as.matrix(df)
}

betahat <- read_mat(logfc_path)
sehat   <- read_mat(se_path)
pval    <- read_mat(pval_path)

# convert p-values to SE via z-score (matching existing Rmd logic)
shat <- betahat / mashr:::p2z(pval, betahat)
shat[is.na(shat)] <- 1e16

data_mash <- mash_set_data(betahat, shat)
m.1by1    <- mash_1by1(data_mash)
strong    <- get_significant_results(m.1by1, lfsr_thr)

U.pca <- cov_pca(data_mash, n_pca, subset = strong)
U.ed  <- cov_ed(data_mash, U.pca, subset = strong)
m.ed  <- mash(data_mash, Ulist = c(U.ed, U.pca))

saveRDS(m.ed, out_rds)
