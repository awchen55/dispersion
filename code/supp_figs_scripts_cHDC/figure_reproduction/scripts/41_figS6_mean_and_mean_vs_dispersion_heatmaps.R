#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(viridisLite)
  library(RColorBrewer)
  library(grid)
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
fig_repro_dir <- file.path(repo_dir, "figure_reproduction")
data_dir <- file.path(fig_repro_dir, "data")
out_dir <- file.path(fig_repro_dir, "outputs")
s6_dir <- file.path(out_dir, "figS6")

mean_csv <- file.path(data_dir, "mean_expression_excl_cluster24.csv")
dispersion_csv <- file.path(data_dir, "res_var_excl_cluster24.csv")
retained_matrix_csv <- file.path(out_dir, "retained_gene_matrix.csv")
# ANALYSIS_OK[external-path]: same original cell-count file as script 37 (sub-cluster 24 not recounted; flagged in header).
cell_counts_csv <- file.path(repo_dir, "cellranger_cluster-mode_trial", "analysis", "datasets",
                             "cell_counts_new_pred_id.csv")

ribo_pattern <- "^RP"
individuals <- c("NA19099", "NA19093", "NA18522")
cell_types <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)
n_expected_cols <- length(individuals) * length(cell_types)  # 18

# ANALYSIS_OK[magic-threshold]: colour-scale limits copied verbatim from the published S6 code (A 0.7-1, B -0.2-0.2); out-of-range values clamp, ranges logged below.
mean_col_range <- c(0.7, 1)
cross_col_range <- c(-0.2, 0.2)

fig_width_in <- 7.2
fig_height_in <- 2.6
# Horizontal split of the canvas: A (with its own legend) | B body + annotation
# names | B's legend column (Correlation, Individual, Cell type stacked).
panel_a_frac <- 0.42
legend_frac <- 0.17

# ------------------------------- LOAD -----------------------------------------
read_matrix <- function(path) {
  if (!file.exists(path)) stop("Missing ", path)
  m <- as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
  if (ncol(m) != n_expected_cols) stop(path, " has ", ncol(m), " cols; expected ", n_expected_cols)
  if (anyNA(m)) stop(path, " contains NA values")
  m
}
mean_mat <- read_matrix(mean_csv)
res_var_mat <- read_matrix(dispersion_csv)
if (!identical(colnames(mean_mat), colnames(res_var_mat))) stop("Mean and dispersion columns differ")
if (!identical(rownames(mean_mat), rownames(res_var_mat))) stop("Mean and dispersion genes differ")

# ANALYSIS_OK[filtering]: drop ^RP ribosomal genes (repo convention, as in 01/03/37); count logged below
keep <- !grepl(ribo_pattern, rownames(res_var_mat))
if (!file.exists(retained_matrix_csv)) stop("Missing ", retained_matrix_csv, " -- run 01 first")
retained_ind <- as.matrix(read.csv(retained_matrix_csv, row.names = 1, check.names = FALSE))
# ANALYSIS_OK[column-reorder]: select cell-type columns in canonical order; no rows dropped
retained_ind <- retained_ind[, cell_types, drop = FALSE]
union_retained_genes <- rownames(retained_ind)[rowSums(retained_ind == 1) == length(cell_types)]
genes <- intersect(rownames(res_var_mat)[keep], union_retained_genes)
if (length(genes) != length(union_retained_genes)) {
  stop("union-retained genes missing from the matrices: ", length(union_retained_genes) - length(genes))
}
message(sprintf("[load] %d genes -> %d non-ribo -> %d union-retained", nrow(res_var_mat), sum(keep), length(genes)))
mean_mat <- mean_mat[genes, , drop = FALSE]
res_var_mat <- res_var_mat[genes, , drop = FALSE]

# ------------------------------- CORRELATIONS ---------------------------------
cormat <- cor(mean_mat)  # Pearson, as in the source
cross_cor_mat <- cor(mean_mat, res_var_mat, method = "spearman")
dir.create(s6_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(cormat, file.path(s6_dir, "mean_correlation_matrix.csv"))
write.csv(cross_cor_mat, file.path(s6_dir, "mean_vs_dispersion_spearman_matrix.csv"))
off_diag <- cormat[upper.tri(cormat)]
same_sample <- diag(cross_cor_mat)
message(sprintf("[A] mean Pearson r off-diagonal range %.3f-%.3f", min(off_diag), max(off_diag)))
message(sprintf("[B] mean-vs-dispersion Spearman range %.3f-%.3f; same-sample range %.3f-%.3f",
                min(cross_cor_mat), max(cross_cor_mat), min(same_sample), max(same_sample)))

# ------------------------------- SAMPLE METADATA ------------------------------
sample_cols <- colnames(res_var_mat)
individual <- factor(sub("_.*$", "", sample_cols), levels = individuals)
cell_type_raw <- sub("^[^_]*_", "", sample_cols)
stopifnot(!anyNA(individual), all(cell_type_raw %in% cell_types))

# Display labels as in the published image: underscores removed, the longest
# name wrapped after its second word.
cell_type_labels <- gsub("_", " ", cell_types)
cell_type_labels <- sub("^(\\S+ \\S+) ", "\\1\n", cell_type_labels)
cell_type <- factor(cell_type_raw, levels = cell_types)
cell_type_cols <- setNames(brewer.pal(12, "Paired")[seq_along(cell_types)], cell_types)
levels(cell_type) <- cell_type_labels
names(cell_type_cols) <- cell_type_labels
individual_cols <- setNames(brewer.pal(8, "Paired")[seq_along(individuals)], individuals)

if (!file.exists(cell_counts_csv)) stop("Missing cell counts file: ", cell_counts_csv)
cell_counts_raw <- read.csv(cell_counts_csv, row.names = 1)
cell_counts <- cell_counts_raw$Freq[match(sample_cols, cell_counts_raw$sample_id)]
if (anyNA(cell_counts)) stop("cell_counts missing sample(s): ", paste(sample_cols[is.na(cell_counts)], collapse = ", "))

# No fontfamily in gpar: cairo_pdf(family = "Arial") already makes Arial the
# default, and naming it here makes ComplexHeatmap measure text with missing
# PostScript metrics, which mis-sizes annotation names (they overlapped the
# legends).
font_annot <- gpar(fontsize = 7)
font_legend_title <- gpar(fontsize = 6)
font_legend_labels <- gpar(fontsize = 6)
legend_key <- unit(3, "mm")

make_annotation <- function(show_names) {
  columnAnnotation(
    Individual = individual,
    Cell_Type = cell_type,
    Cell_Count = anno_barplot(cell_counts, gp = gpar(fill = "grey50", col = "black"),
                              axis_param = list(gp = font_legend_labels)),
    col = list(Individual = individual_cols, Cell_Type = cell_type_cols),
    annotation_label = c(Individual = "Individual", Cell_Type = "Cell type", Cell_Count = "Cell count"),
    show_annotation_name = show_names,
    annotation_name_gp = font_annot,
    annotation_legend_param = list(
      Individual = list(title = "Individual", title_gp = font_legend_title, labels_gp = font_legend_labels,
                        grid_height = legend_key, grid_width = legend_key),
      Cell_Type = list(title = "Cell type", title_gp = font_legend_title, labels_gp = font_legend_labels,
                       grid_height = legend_key, grid_width = legend_key)
    )
  )
}

make_heatmap <- function(mat, col_range, annotation, show_row_dend) {
  Heatmap(
    mat,
    top_annotation = annotation,
    name = "Correlation",
    col = colorRamp2(seq(col_range[1], col_range[2], length.out = 11), inferno(11)),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_row_dend = show_row_dend,
    show_column_dend = TRUE,
    heatmap_legend_param = list(title_gp = font_legend_title, labels_gp = font_legend_labels,
                                direction = "vertical", grid_width = legend_key,
                                legend_height = unit(18, "mm"))
  )
}

# ------------------------------- COMPOSITE ------------------------------------
# Heatmaps are built AFTER opening cairo_pdf: ComplexHeatmap measures text at
# construction, and the default device has no Arial metrics.
out_pdf <- file.path(s6_dir, "FigureS6_composite.pdf")
cairo_pdf(out_pdf, width = fig_width_in, height = fig_height_in, family = "Arial")
hm_a <- make_heatmap(cormat, mean_col_range, make_annotation(show_names = FALSE), show_row_dend = FALSE)
hm_b <- make_heatmap(cross_cor_mat, cross_col_range, make_annotation(show_names = TRUE), show_row_dend = FALSE)
grid.newpage()
panel_tag <- function(lbl) {
  grid.text(lbl, x = unit(1, "mm"), y = unit(1, "npc") - unit(1, "mm"), just = c("left", "top"),
            gp = gpar(fontsize = 12, fontface = "bold"))
}
pushViewport(viewport(x = 0, width = panel_a_frac, just = "left"))
draw(hm_a, newpage = FALSE, show_annotation_legend = FALSE, heatmap_legend_side = "right",
     padding = unit(c(2, 6, 2, 2), "mm"))
panel_tag("A")
popViewport()
# B is drawn without legends so its right-hand annotation names get the space;
# its legends are built explicitly and stacked in their own column.
pushViewport(viewport(x = panel_a_frac, width = 1 - panel_a_frac - legend_frac, just = "left"))
draw(hm_b, newpage = FALSE, show_heatmap_legend = FALSE, show_annotation_legend = FALSE,
     padding = unit(c(2, 6, 2, 2), "mm"))
panel_tag("B")
popViewport()
legends_b <- packLegend(
  Legend(title = "Correlation", col_fun = colorRamp2(seq(cross_col_range[1], cross_col_range[2], length.out = 11), inferno(11)),
         title_gp = font_legend_title, labels_gp = font_legend_labels,
         grid_width = legend_key, legend_height = unit(18, "mm")),
  Legend(title = "Individual", at = individuals, legend_gp = gpar(fill = individual_cols[individuals]),
         title_gp = font_legend_title, labels_gp = font_legend_labels,
         grid_height = legend_key, grid_width = legend_key),
  Legend(title = "Cell type", at = cell_type_labels, legend_gp = gpar(fill = cell_type_cols[cell_type_labels]),
         title_gp = font_legend_title, labels_gp = font_legend_labels,
         grid_height = legend_key, grid_width = legend_key),
  direction = "vertical", gap = unit(2, "mm"))
pushViewport(viewport(x = 1 - legend_frac, width = legend_frac, just = "left"))
draw(legends_b, x = unit(1, "mm"), just = "left")
popViewport()
dev.off()
message("[done] wrote ", out_pdf)

# ---- reportable values (same format as 39) -----------------------------------
reportable <- data.frame(
  key = c("figS6_n_genes", "figS6_mean_r_min", "figS6_mean_r_max",
          "figS6_mean_vs_disp_rho_min", "figS6_mean_vs_disp_rho_max",
          "figS6_mean_vs_disp_same_sample_rho_median"),
  value = c(length(genes), round(min(off_diag), 4), round(max(off_diag), 4),
            round(min(cross_cor_mat), 4), round(max(cross_cor_mat), 4), round(median(same_sample), 4)),
  provenance = "figure_reproduction/scripts/41_figS6_mean_and_mean_vs_dispersion_heatmaps.R")
write.csv(reportable, file.path(out_dir, "reportable_values_41.csv"), row.names = FALSE)
