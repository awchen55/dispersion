library(variancePartition)

fit_paths <- snakemake@input[["fits"]]
celltypes <- snakemake@params[["celltypes"]]

out_logfc    <- snakemake@output[["logfc"]]
out_se       <- snakemake@output[["se"]]
out_pval     <- snakemake@output[["pval"]]
out_adj_pval <- snakemake@output[["adj_pval"]]

# Build one data.frame per cell type, then full-outer-join
merge_full <- function(existing, new_df) {
  if (is.null(existing)) return(new_df)
  merge(existing, new_df, by = "gene", all = TRUE)
}

logfc_df <- se_df <- pval_df <- adj_df <- NULL

for (i in seq_along(celltypes)) {
  ct  <- celltypes[i]
  ct2 <- gsub("_", ".", ct)
  fit <- readRDS(fit_paths[i])
  res <- variancePartition::topTable(fit, coef = "L1", number = Inf)
  res$SE <- res$logFC / res$t
  genes  <- rownames(res)

  logfc_df <- merge_full(logfc_df,
    setNames(data.frame(genes, res$logFC,  stringsAsFactors = FALSE),
             c("gene", paste0("de.logFC.diploid.", ct2))))
  se_df    <- merge_full(se_df,
    setNames(data.frame(genes, res$SE,      stringsAsFactors = FALSE),
             c("gene", paste0("de.SE.diploid.", ct2))))
  pval_df  <- merge_full(pval_df,
    setNames(data.frame(genes, res$P.Value, stringsAsFactors = FALSE),
             c("gene", paste0("de.pvalue.diploid.", ct2))))
  adj_df   <- merge_full(adj_df,
    setNames(data.frame(genes, res$adj.P.Val, stringsAsFactors = FALSE),
             c("gene", paste0("de.adj.pvalue.diploid.", ct2))))
}

write_mat <- function(df, path) {
  rownames(df) <- df$gene; df$gene <- NULL
  write.csv(df, path, row.names = TRUE)
}
write_mat(logfc_df, out_logfc)
write_mat(se_df,    out_se)
write_mat(pval_df,  out_pval)
write_mat(adj_df,   out_adj_pval)
