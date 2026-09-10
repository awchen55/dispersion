library(variancePartition)

fit_paths <- snakemake@input[["fits"]]
celltypes <- snakemake@params[["celltypes"]]

out_logfc <- snakemake@output[["logfc"]]
out_se    <- snakemake@output[["se"]]
out_pval  <- snakemake@output[["pval"]]

merge_full <- function(existing, new_df) {
  if (is.null(existing)) return(new_df)
  merged <- merge(existing, new_df, by = "row.names", all = TRUE)
  rownames(merged) <- merged$Row.names
  merged$Row.names <- NULL
  merged
}

logfc_df <- se_df <- pval_df <- NULL

for (i in seq_along(celltypes)) {
  ct  <- celltypes[i]
  ct2 <- gsub("_", ".", ct)
  fit <- readRDS(fit_paths[i])
  res <- topTable(fit, coef = "specieshuman - specieschimp", number = Inf)
  res$SE <- res$logFC / res$t

  make_df <- function(vals, prefix) {
    df <- data.frame(vals, row.names = rownames(res))
    colnames(df) <- paste0(prefix, ct2)
    df
  }

  logfc_df <- merge_full(logfc_df, make_df(res$logFC,   "dd.logFC.diploid."))
  se_df    <- merge_full(se_df,    make_df(res$SE,       "dd.SE.diploid."))
  pval_df  <- merge_full(pval_df,  make_df(res$P.Value,  "dd.pvalue.diploid."))
}

write.csv(logfc_df, out_logfc, row.names = TRUE)
write.csv(se_df,    out_se,    row.names = TRUE)
write.csv(pval_df,  out_pval,  row.names = TRUE)
