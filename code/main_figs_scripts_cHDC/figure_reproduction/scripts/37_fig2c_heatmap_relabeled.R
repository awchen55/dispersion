#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(viridisLite)
  library(grid)
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
fig_repro_dir <- file.path(repo_dir, "figure_reproduction")
data_dir <- file.path(fig_repro_dir, "data")
cluster_dir <- file.path(repo_dir, "cellranger_cluster-mode_trial")
analysis_dir <- fig_repro_dir
out_dir <- file.path(analysis_dir, "outputs")
fig2_dir <- file.path(out_dir, "fig2_replication")

dispersion_csv <- file.path(data_dir, "res_var_excl_cluster24.csv")
# ANALYSIS_OK[external-path]: unaffected by the Epithelial sub-cluster-24
# exclusion (see header note above) -- left at its original location.
cell_counts_csv <- file.path(cluster_dir, "analysis", "datasets",
                              "cell_counts_new_pred_id.csv")
retained_matrix_csv <- file.path(out_dir, "retained_gene_matrix.csv")

ribo_pattern <- "^RP"
individuals <- c("NA19099", "NA19093", "NA18522")
cell_types <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)
n_expected_cols <- length(individuals) * length(cell_types)  # 18

# One word per line (after underscore->space).
wrap_by_word <- function(x) {
  parts <- strsplit(x, " ")
  vapply(parts, paste, character(1), collapse = "\n")
}
cell_type_labels <- wrap_by_word(gsub("_", " ", cell_types))

# ------------------------------- LOAD -----------------------------------------
if (!file.exists(dispersion_csv)) stop("Missing dispersion matrix: ", dispersion_csv)
res_var_df <- read.csv(dispersion_csv, row.names = 1, check.names = FALSE)
if (ncol(res_var_df) != n_expected_cols) {
  stop("Dispersion matrix has ", ncol(res_var_df), " cols; expected ", n_expected_cols)
}
res_var_mat <- as.matrix(res_var_df)
# ANALYSIS_OK[filtering]: drop ^RP ribosomal genes (repo convention, matches script 01/03/25); retained count logged below
res_var_mat <- res_var_mat[!grepl(ribo_pattern, rownames(res_var_mat)), , drop = FALSE]

if (!file.exists(retained_matrix_csv)) {
  stop("Missing retained_gene_matrix.csv -- run 01_exclude_topexpressed_genes.R first")
}
retained_ind <- as.matrix(read.csv(retained_matrix_csv, row.names = 1, check.names = FALSE))
retained_ind <- retained_ind[, cell_types, drop = FALSE]
union_retained_genes <- rownames(retained_ind)[rowSums(retained_ind == 1) == length(cell_types)]
res_var_mat <- res_var_mat[union_retained_genes, , drop = FALSE]

# ------------------------------- CORRELATION -----------------------------------
cormat <- cor(res_var_mat, method = "pearson")

# ------------------------------- SAMPLE METADATA --------------------------------
sample_cols <- colnames(res_var_mat)
individual <- sub("_.*$", "", sample_cols)
cell_type  <- sub("^[^_]*_", "", sample_cols)
stopifnot(all(individual %in% individuals), all(cell_type %in% cell_types))
individual <- factor(individual, levels = individuals)
# Factor built with the ORIGINAL (underscore) level names first (so matching
# against cell_type_cols below stays simple), then levels renamed to the
# wrapped/cleaned display text.
cell_type <- factor(cell_type, levels = cell_types)
cell_type_cols <- setNames(RColorBrewer::brewer.pal(12, "Paired")[seq_along(cell_types)], cell_types)
levels(cell_type) <- cell_type_labels
names(cell_type_cols) <- cell_type_labels

if (!file.exists(cell_counts_csv)) stop("Missing cell counts file: ", cell_counts_csv)
cell_counts_raw <- read.csv(cell_counts_csv, row.names = 1)
cell_counts <- cell_counts_raw[match(sample_cols, cell_counts_raw$sample_id), ]
if (any(is.na(cell_counts$Freq))) {
  stop("cell_counts missing an entry for sample(s): ",
       paste(sample_cols[is.na(cell_counts$Freq)], collapse = ", "))
}

individual_cols <- setNames(RColorBrewer::brewer.pal(8, "Paired")[seq_along(individuals)], individuals)

# --- Annotation: names relabeled via annotation_label (internal identifiers
#     Individual/Cell_Type/Cell_Count kept, only display text changes) -------
ha <- columnAnnotation(
  Individual = individual,
  Cell_Type  = cell_type,
  Cell_Count = anno_barplot(
    cell_counts$Freq,
    gp = gpar(fill = "grey30", col = NA),
    border = FALSE
  ),
  col = list(
    Individual = individual_cols,
    Cell_Type  = cell_type_cols
  ),
  annotation_label = c(Individual = "Individual", Cell_Type = "Cell Type", Cell_Count = "Cell Count"),
  gp = gpar(fontsize = 8, fontfamily = "Arial"),
  annotation_name_gp = gpar(fontsize = 8, fontfamily = "Arial"),
  annotation_legend_param = list(
    title_gp  = gpar(fontsize = 8, fontfamily = "Arial"),
    labels_gp = gpar(fontsize = 8, fontfamily = "Arial"),
    direction = "vertical",
    ncol = 1
  )
)

# --- Heatmap (byte-for-byte copy of script 26's spec, except font sizes) -------
HCheat <- Heatmap(
  cormat,
  top_annotation = ha,
  name = "Correlation",
  col = circlize::colorRamp2(
    seq(0.8, 1, length.out = 11),
    inferno(11)
  ),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = TRUE,
  column_title_gp = gpar(fontsize = 9, fontfamily = "Arial"),
  row_title_gp = gpar(fontsize = 9, fontfamily = "Arial"),
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = 8, fontfamily = "Arial"),
    labels_gp = gpar(fontsize = 8, fontfamily = "Arial"),
    direction = "vertical"
  )
)

out_pdf <- file.path(fig2_dir, "Fig2a_ResVar_CorrelationHeatmap_relabeled.pdf")
cairo_pdf(out_pdf, width = 6.2, height = 5.12, family = "Arial")
# align_*_legend = "heatmap_center": vertically centers the legend stack
# against the full heatmap body, clear of the top annotation row labels.
draw(HCheat, heatmap_legend_side = "right", annotation_legend_side = "right",
     align_heatmap_legend = "heatmap_center", align_annotation_legend = "heatmap_center")
dev.off()
message("[done] wrote ", out_pdf)
