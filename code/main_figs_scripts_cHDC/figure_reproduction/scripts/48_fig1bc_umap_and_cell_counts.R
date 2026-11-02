#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
seurat_rds <- file.path(repo_dir, "cellranger_cluster-mode_trial", "analysis", "datasets",
                        "100k_full_filtered_norm_seurat.rds")
out_dir <- file.path(repo_dir, "figure_reproduction", "outputs", "fig1_replication")
min_cells_per_sample <- 1000   # source: pass_threshold = Freq >= 1000 (computed, not plotted)

if (!file.exists(seurat_rds)) stop("Missing Seurat object: ", seurat_rds)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("[48] Reading ", seurat_rds)
brendan_obj <- readRDS(seurat_rds)
for (col in c("new.predicted.id", "vireo.individual")) {
  if (!col %in% colnames(brendan_obj@meta.data)) stop("Missing metadata column ", col)
}
message(sprintf("[48] %d cells", ncol(brendan_obj)))

# ---- Fig 1B: UMAP -------------------------------------------------------------
brendan_obj$new.predicted.id <- factor(brendan_obj$new.predicted.id)
brendan_obj$new.predicted.id <- droplevels(brendan_obj$new.predicted.id)
brendan_obj$clean.new.pred.id <- gsub("_", " ", brendan_obj$new.predicted.id)
brendan_obj$clean.new.pred.id <- factor(brendan_obj$clean.new.pred.id)

celltype_colors <- c(
  "Activated fibroblasts" = "#117733",
  "Cardiac fibroblasts" = "#4477AA",
  "Cardiomyocytes" = "#882255",
  "Endothelial cells" = "#C45F00",
  "Epicardial cells" = "#661100",
  "Epithelial cells" = "#AA4499",
  "Foregut endodermal cells" = "#CC6677",
  "Hematopoietic cells" = "#999933",
  "Hepatic endodermal cells" = "#44AA99",
  "Mesenchymal cells" = "#DDCC77",
  "Neural crest cells" = "#332288",
  "Neural progenitor cells" = "#AA7744",
  "Neuroectodermal cells" = "#6699CC",
  "Proliferating fibroblasts" = "gold1"
)
missing_cols <- setdiff(levels(brendan_obj$clean.new.pred.id), names(celltype_colors))
if (length(missing_cols) > 0) stop("No colour for cell types: ", paste(missing_cols, collapse = ", "))

fig1b <- DimPlot(
  brendan_obj,
  reduction = "umap",
  group.by = "clean.new.pred.id",
  cols = celltype_colors,
  label = FALSE,
  pt.size = 0.01,
  raster = FALSE
) +
  theme_classic(base_family = "Arial") +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 5),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 6),
    text = element_text(size = 6),
    plot.margin = margin(2, 2, 2, 2, unit = "mm"),
    legend.key.size = unit(2.5, "mm"),
    legend.spacing.y = unit(1, "mm")
  ) +
  labs(x = "UMAP 1", y = "UMAP 2", color = NULL)
ggsave(file.path(out_dir, "Fig1b_UMAP_CellTypes.pdf"), fig1b, device = cairo_pdf,
       width = 4, height = 2, units = "in")
message("[48] Wrote Fig1b_UMAP_CellTypes.pdf")

# ---- Fig 1C: cell counts --------------------------------------------------------
cell_counts <- as.data.frame(table(brendan_obj$new.predicted.id, brendan_obj$vireo.individual))
cell_counts <- cell_counts %>%
  mutate(pass_threshold = Freq >= min_cells_per_sample) %>%
  mutate(clean_ct = gsub("_", " ", Var1)) %>%
  mutate(sample_id = paste0(Var2, "_", Var1))
colnames(cell_counts) <- c("cell_type", "cell_line", "count", "pass_threshold", "clean_cell_type", "sample_id")
if (sum(cell_counts$count) != ncol(brendan_obj)) stop("Cell counts do not sum to the number of cells")
write.csv(cell_counts, file.path(out_dir, "fig1c_cell_counts.csv"), row.names = FALSE)

cell_line_colors <- c(
  "NA19099" = "#117733",
  "NA19093" = "#4477AA",
  "NA18522" = "#B22222"
)
if (!all(unique(as.character(cell_counts$cell_line)) %in% names(cell_line_colors))) stop("Unexpected individuals")
cell_counts$clean_cell_type <- factor(cell_counts$clean_cell_type)

p_cells_stacked <- ggplot(cell_counts, aes(x = clean_cell_type, y = count, fill = cell_line)) +
  geom_col(width = 0.8, linewidth = 0) +
  scale_fill_manual(values = cell_line_colors) +
  labs(x = "Cell type", y = "Number of cells", fill = "Individual") +
  theme_classic(base_family = "Arial") +
  theme(
    text = element_text(size = 6),
    axis.title = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.line = element_line(linewidth = 0.4),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    legend.position = "right",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.key.size = unit(6, "pt"),
    plot.margin = margin(2, 2, 2, 2, unit = "mm")
  )
ggsave(file.path(out_dir, "Fig1c_cell_count_stacked_barplot.pdf"), p_cells_stacked, device = cairo_pdf,
       width = 3.5, height = 3.2, units = "in")
message(sprintf("[48] Wrote Fig1c_cell_count_stacked_barplot.pdf (%d cell type x individual groups, %d pass >= %d cells)",
                nrow(cell_counts), sum(cell_counts$pass_threshold), min_cells_per_sample))
