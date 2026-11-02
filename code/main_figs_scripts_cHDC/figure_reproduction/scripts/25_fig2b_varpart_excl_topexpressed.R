#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(variancePartition)
  library(ggplot2)
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
fig_repro_dir <- file.path(repo_dir, "figure_reproduction")
data_dir <- file.path(fig_repro_dir, "data")
cluster_dir <- file.path(repo_dir, "cellranger_cluster-mode_trial")
analysis_dir <- fig_repro_dir
out_dir <- file.path(analysis_dir, "outputs")
fig2_dir <- file.path(out_dir, "fig2_replication")
dir.create(fig2_dir, recursive = TRUE, showWarnings = FALSE)

dispersion_csv <- file.path(data_dir, "res_var_excl_cluster24.csv")
# Epithelial_cells rows recomputed on sub-cluster-24-excluded cells (see header).
aggr_metadata_csv <- file.path(data_dir, "aggregated_metadata_excl_cluster24.csv")
retained_matrix_csv <- file.path(out_dir, "retained_gene_matrix.csv")
published_vp_csv <- file.path(cluster_dir, "analysis", "datasets", "limma",
                               "var_partition_new_pred_id.csv")

ribo_pattern <- "^RP"
individuals <- c("NA19099", "NA19093", "NA18522")
cell_types <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)
n_expected_cols <- length(individuals) * length(cell_types)  # 18

# ------------------------------- LOAD DISPERSION MATRIX -------------------------
if (!file.exists(dispersion_csv)) stop("Missing dispersion matrix: ", dispersion_csv)
res_var_df <- read.csv(dispersion_csv, row.names = 1, check.names = FALSE)
if (ncol(res_var_df) != n_expected_cols) {
  stop("Dispersion matrix has ", ncol(res_var_df), " cols; expected ", n_expected_cols)
}
res_var_mat <- as.matrix(res_var_df)
# ANALYSIS_OK[filtering]: drop ^RP ribosomal genes (repo convention, matches script 01/03/37); retained count logged below
res_var_mat <- res_var_mat[!grepl(ribo_pattern, rownames(res_var_mat)), , drop = FALSE]
message(sprintf("[load] %d non-ribo genes x %d samples", nrow(res_var_mat), ncol(res_var_mat)))

if (!file.exists(retained_matrix_csv)) {
  stop("Missing retained_gene_matrix.csv -- run 01_exclude_topexpressed_genes.R first")
}
retained_ind <- as.matrix(read.csv(retained_matrix_csv, row.names = 1, check.names = FALSE))
retained_ind <- retained_ind[, cell_types, drop = FALSE]

# Union-retained gene set: genes excluded in NO cell type (same rule as scripts 03/37).
union_retained_genes <- rownames(retained_ind)[rowSums(retained_ind == 1) == length(cell_types)]
message(sprintf("[excl] union-retained genes (excluded in no cell type): %d of %d non-ribo",
                length(union_retained_genes), nrow(res_var_mat)))

res_var_mat <- res_var_mat[union_retained_genes, , drop = FALSE]

# ------------------------------- SAMPLE METADATA --------------------------------
if (!file.exists(aggr_metadata_csv)) stop("Missing aggregated metadata: ", aggr_metadata_csv)
aggr_means_df <- read.csv(aggr_metadata_csv)

sample_cols <- colnames(res_var_mat)
aggr_means_df <- aggr_means_df[match(sample_cols, aggr_means_df$sample_id), ]
if (any(is.na(aggr_means_df$sample_id))) {
  stop("aggregated metadata missing an entry for sample(s): ",
       paste(sample_cols[is.na(aggr_means_df$sample_id)], collapse = ", "))
}
stopifnot(identical(aggr_means_df$sample_id, sample_cols))

metadata <- data.frame(
  cell_type = factor(aggr_means_df$new.predicted.id, levels = cell_types),
  cell_line = factor(aggr_means_df$vireo.individual, levels = individuals),
  cell_count = scale(aggr_means_df$cell_num),
  nFeatures = scale(aggr_means_df$nFeature_RNA),
  perc_mt = scale(aggr_means_df$percent.mt),
  row.names = sample_cols
)

# ------------------------------- VARIANCE PARTITION MODEL -----------------------
# Formula matches the published ground-truth CSV's columns (cell_line, cell_type,
# cell_count, nFeatures, perc_mt, Residuals) -- see header note on the formula
# discrepancy between the two reference notebooks.
form <- ~ cell_count + nFeatures + perc_mt + (1 | cell_line) + (1 | cell_type)

vp_res_var <- fitExtractVarPartModel(res_var_mat, form, metadata)
write.csv(vp_res_var, file.path(fig2_dir, "var_partition_retained.csv"), row.names = TRUE)

retained_median_ct <- median(vp_res_var$cell_type)
message(sprintf("[varpart] retained-set median cell_type variance = %.4f", retained_median_ct))

if (file.exists(published_vp_csv)) {
  published_vp <- read.csv(published_vp_csv, row.names = 1)
  published_median_ct <- median(published_vp$cell_type)
  message(sprintf("[varpart] published (full universe) median cell_type variance = %.4f",
                  published_median_ct))
} else {
  published_median_ct <- NA_real_
  message("[varpart] NOTE: published ", basename(published_vp_csv),
          " not found; skipping comparison")
}

# ------------------------------- PLOT (byte-for-byte copy of source Rmd spec) --
vp_plot_res_var <- plotVarPart(sortCols(vp_res_var)) +
  theme_classic(base_family = "Arial") +
  theme(
    text = element_text(size = 6),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    legend.position = "right",
    panel.grid = element_blank(),
    plot.margin = margin(2, 2, 2, 2, unit = "mm")
  )

out_pdf <- file.path(fig2_dir, "Fig2b_VarPart_ResVar_excl_topexpressed.pdf")
ggsave(
  filename = out_pdf,
  plot = vp_plot_res_var,
  device = cairo_pdf,
  width = 3,
  height = 2.5,
  units = "in"
)

# ------------------------------- REPORTABLE VALUES ------------------------------
write.csv(
  data.frame(
    key = c("fig2b_varpart_retained_n_genes",
            "fig2b_varpart_retained_median_celltype_variance",
            "fig2b_varpart_published_median_celltype_variance"),
    value = c(length(union_retained_genes), round(retained_median_ct, 4),
              if (is.na(published_median_ct)) NA else round(published_median_ct, 4)),
    provenance = "figure_reproduction/scripts/25_fig2b_varpart_excl_topexpressed.R"
  ),
  file.path(out_dir, "reportable_values_25.csv"),
  row.names = FALSE
)

message("[done] wrote ", out_pdf)
