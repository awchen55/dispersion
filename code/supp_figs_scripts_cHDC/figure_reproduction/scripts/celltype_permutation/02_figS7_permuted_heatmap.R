#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(viridisLite)  # inferno
  library(RColorBrewer)
  library(grid)
})

repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
perm_dir <- file.path(repo_dir, "figure_reproduction", "outputs", "celltype_permutation")
res_var_csv <- file.path(perm_dir, "memento", "res_var_permuted.csv")
fig_dir <- file.path(perm_dir, "figS7")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

all_individuals <- c("NA19099", "NA19093", "NA18522")
all_cell_types <- c(
  "Cardiomyocytes",
  "Mesenchymal_cells",
  "Cardiac_fibroblasts",
  "Epithelial_cells",
  "Activated_fibroblasts",
  "Foregut_endodermal_cells"
)

if (!file.exists(res_var_csv)) stop("Missing ", res_var_csv, " -- run 01 first")
perm_res_var <- read.csv(res_var_csv, row.names = 1, check.names = FALSE)
if (ncol(perm_res_var) != length(all_individuals) * length(all_cell_types)) {
  stop("Expected 18 columns, got ", ncol(perm_res_var))
}
if (anyNA(perm_res_var)) stop("NA values in permuted dispersion matrix")
message(sprintf("[load] %d genes x %d samples", nrow(perm_res_var), ncol(perm_res_var)))

split_parts <- strsplit(colnames(perm_res_var), "_", fixed = TRUE)

cormat <- cor(perm_res_var, method = "pearson")
write.csv(cormat, file.path(fig_dir, "permuted_correlation_matrix.csv"))

# Extract metadata
individual <- sapply(split_parts, `[`, 1)
cell_type <- sapply(split_parts, function(x) paste(x[-1], collapse = "_"))
stopifnot(all(individual %in% all_individuals), all(cell_type %in% all_cell_types))

# Force stable factor levels
individual <- factor(individual, levels = all_individuals)
cell_type <- factor(cell_type, levels = all_cell_types)

# Fixed color mappings
individual_palette <- RColorBrewer::brewer.pal(8, "Paired")
individual_cols <- setNames(individual_palette[seq_along(all_individuals)],
                            all_individuals)

cell_type_palette <- RColorBrewer::brewer.pal(12, "Paired")
cell_type_cols <- setNames(cell_type_palette[seq_along(all_cell_types)],
                           all_cell_types)

# Display labels as in the published image.
cell_type_display <- factor(cell_type, levels = all_cell_types)
levels(cell_type_display) <- gsub("_", " ", all_cell_types)
cell_type_display_cols <- setNames(cell_type_cols, gsub("_", " ", all_cell_types))

# --- Annotation ---
ha <- columnAnnotation(
  Individual = individual,
  Cell_type = cell_type_display,
  col = list(
    Individual = individual_cols,
    Cell_type = cell_type_display_cols
  ),
  annotation_label = c(Individual = "Individual", Cell_type = "Cell type"),
  gp = gpar(fontsize = 6, fontfamily = "Arial"),
  annotation_name_gp = gpar(fontsize = 6, fontfamily = "Arial"),
  annotation_legend_param = list(
    Individual = list(title = "Individual"),
    Cell_type = list(title = "Cell type"),
    title_gp = gpar(fontsize = 6, fontfamily = "Arial"),
    labels_gp = gpar(fontsize = 6, fontfamily = "Arial"),
    direction = "vertical",
    ncol = 1
  )
)

# --- Heatmap ---
# ANALYSIS_OK[magic-threshold]: colour-scale limits 0.8-1 copied verbatim from the published S7 chunk; values below 0.8 are clamped to the lowest colour (min logged below).
message(sprintf("[cor] off-diagonal range %.4f-%.4f",
                min(cormat[upper.tri(cormat)]), max(cormat[upper.tri(cormat)])))
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
  column_title_gp = gpar(fontsize = 7, fontfamily = "Arial"),
  row_title_gp = gpar(fontsize = 7, fontfamily = "Arial"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 6, fontfamily = "Arial"),
    labels_gp = gpar(fontsize = 6, fontfamily = "Arial"),
    direction = "vertical"
  )
)

cairo_pdf(file.path(fig_dir, "permut_cellid_CorrelationHeatmap.pdf"),
          width = 4,
          height = 3,
          family = "Arial")
draw(HCheat,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()

# --- Pair summary (added) ---
pairs <- which(upper.tri(cormat), arr.ind = TRUE)
same_ct <- cell_type[pairs[, 1]] == cell_type[pairs[, 2]]
same_ind <- individual[pairs[, 1]] == individual[pairs[, 2]]
r <- cormat[pairs]
pair_class <- ifelse(same_ct & !same_ind, "same_label_diff_individual",
              ifelse(!same_ct & same_ind, "diff_label_same_individual",
                     "diff_label_diff_individual"))
stopifnot(!any(same_ct & same_ind))  # each label x individual is one sample
summary_df <- do.call(rbind, lapply(split(r, pair_class), function(x) {
  data.frame(n_pairs = length(x), mean_r = mean(x), min_r = min(x), max_r = max(x))
}))
summary_df <- cbind(pair_class = rownames(summary_df), summary_df)
rownames(summary_df) <- NULL
write.csv(summary_df, file.path(fig_dir, "permuted_correlation_pair_summary.csv"),
          row.names = FALSE)
print(summary_df)

message("[done] outputs in ", fig_dir)
