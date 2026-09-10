library(mashr)

logfc_path     <- snakemake@input[["logfc"]]
se_path        <- snakemake@input[["se"]]
pval_path      <- snakemake@input[["pval"]]
gene_list_path <- snakemake@input[["gene_list"]]
out_rds        <- snakemake@output[["mash_rds"]]
n_pca          <- snakemake@params[["pca_components"]]
lfsr_thr       <- snakemake@params[["lfsr_threshold"]]

read_mat <- function(path) {
  df <- read.csv(path, row.names = 1, check.names = FALSE)
  as.matrix(df)
}

betahat <- read_mat(logfc_path)
sehat   <- read_mat(se_path)
pval    <- read_mat(pval_path)

# subset to genes in the diploid-derived filtered gene list
gene_list <- read.csv(gene_list_path)$gene
keep      <- intersect(gene_list, rownames(betahat))
betahat   <- betahat[keep, , drop = FALSE]
pval      <- pval[keep, , drop = FALSE]

# impute NAs arising from outer-join merge (CTs where gene was not tested)
betahat[is.na(betahat)] <- 0
pval[is.na(pval)]       <- 1

# convert p-values to SE via z-score (matching diploid MASH logic)
pval[pval == 0] <- 1e-300
shat <- betahat / mashr:::p2z(pval, betahat)
shat[is.na(shat) | !is.finite(shat)] <- 1e16

data_mash <- mash_set_data(betahat, shat)
m.1by1    <- mash_1by1(data_mash)
strong    <- get_significant_results(m.1by1, lfsr_thr)

U.pca <- cov_pca(data_mash, n_pca, subset = strong)
U.ed  <- cov_ed(data_mash, U.pca, subset = strong)
m.ed  <- mash(data_mash, Ulist = c(U.ed, U.pca))

saveRDS(m.ed, out_rds)
