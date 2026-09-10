library(dplyr)

dip_logfc_path   <- snakemake@input[["dip_dd_logfc"]]
dip_se_path      <- snakemake@input[["dip_dd_se"]]
dip_pval_path    <- snakemake@input[["dip_dd_pval"]]
tet_logfc_path   <- snakemake@input[["tet_dd_logfc"]]

out_logfc   <- snakemake@output[["filtered_logfc"]]
out_se      <- snakemake@output[["filtered_se"]]
out_pval    <- snakemake@output[["filtered_pval"]]
out_genes   <- snakemake@output[["gene_list"]]

threshold <- snakemake@params[["filter_threshold"]]

read_mat <- function(path) {
  df <- read.csv(path, row.names = 1, check.names = FALSE)
  as.matrix(df)
}

dip_logfc <- read_mat(dip_logfc_path)
dip_se    <- read_mat(dip_se_path)
dip_pval  <- read_mat(dip_pval_path)
tet_logfc <- read_mat(tet_logfc_path)

# --- Step 1: tetraploid filter (applied first) ---
# A gene is "tested" in a cell type if its logFC != 0.
# Keep genes tested in >= threshold fraction of tetraploid cell types.
n_tet_ct      <- ncol(tet_logfc)
min_tet       <- floor(threshold * n_tet_ct)
tet_tested    <- rowSums(tet_logfc != 0, na.rm = TRUE)
genes_tet_70  <- rownames(tet_logfc)[tet_tested >= min_tet]

message(sprintf(
  "Tetraploid filter (>= %.0f%% of %d CTs): %d / %d genes pass",
  threshold * 100, n_tet_ct, length(genes_tet_70), nrow(tet_logfc)
))

# --- Step 2: diploid filter ---
# Restrict to tetraploid-passing genes first, then apply threshold.
dip_logfc_sub <- dip_logfc[rownames(dip_logfc) %in% genes_tet_70, , drop = FALSE]
n_dip_ct      <- ncol(dip_logfc_sub)
min_dip       <- floor(threshold * n_dip_ct)
dip_tested    <- rowSums(dip_logfc_sub != 0, na.rm = TRUE)
final_genes   <- rownames(dip_logfc_sub)[dip_tested >= min_dip]

message(sprintf(
  "Diploid filter (>= %.0f%% of %d CTs, from tetraploid-passing set): %d / %d genes pass",
  threshold * 100, n_dip_ct, length(final_genes), nrow(dip_logfc_sub)
))
message(sprintf("Final gene list: %d genes", length(final_genes)))

# --- Subset and impute diploid matrices for MASH ---
subset_impute <- function(mat, genes, na_val) {
  m <- mat[genes, , drop = FALSE]
  m[is.na(m)] <- na_val
  m
}

logfc_out <- subset_impute(dip_logfc, final_genes, 0)
se_out    <- subset_impute(dip_se,    final_genes, 1e16)
pval_out  <- subset_impute(dip_pval,  final_genes, 1)

write.csv(logfc_out, out_logfc, row.names = TRUE)
write.csv(se_out,    out_se,    row.names = TRUE)
write.csv(pval_out,  out_pval,  row.names = TRUE)
write.csv(data.frame(gene = final_genes), out_genes, row.names = FALSE)
