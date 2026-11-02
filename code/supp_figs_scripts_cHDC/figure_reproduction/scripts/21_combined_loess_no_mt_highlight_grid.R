#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
fig_repro_dir <- file.path(repo_dir, "figure_reproduction")
gene_level_path <- file.path(fig_repro_dir, "outputs", "loess_mt_diagnostic",
                              "mt_gene_level_residuals.csv")
out_dir <- file.path(fig_repro_dir, "outputs", "loess_no_mt_highlight")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(gene_level_path)) {
  stop("Missing ", gene_level_path,
       " -- run 06_loess_mean_dispersion_diagnostic.R first.")
}

# ANALYSIS_OK[experimental-design]: row/column order -- individuals as rows,
# cell types as columns (matches the source script's explicit layout choice);
# matches 06's fixed 6-cell-type x 3-individual design either way.
individuals <- c("NA19099", "NA19093", "NA18522")
cell_types <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)
# ANALYSIS_OK[model-formula]: span=0.5, matching 06's loess_span -- the
# geom_smooth() call below re-fits per facet (individual x cell_type) rather
# than reusing 06's stored fitted_all curve, because ggplot's LOESS ribbon
# (se=TRUE) needs the fit object itself, not just point predictions; fitted
# values are numerically identical to 06's since the input data and span
# match exactly.
loess_span <- 0.5

gene_level_df <- read.csv(gene_level_path, stringsAsFactors = FALSE)
stopifnot(
  "mt_gene_level_residuals.csv must contain all 6 expected cell types" =
    setequal(unique(gene_level_df$cell_type), cell_types),
  "mt_gene_level_residuals.csv must contain all 3 expected individuals" =
    setequal(unique(gene_level_df$individual), individuals),
  "mt_gene_level_residuals.csv must contain no non-finite mean_expr/dispersion" =
    all(is.finite(gene_level_df$mean_expr)) && all(is.finite(gene_level_df$dispersion))
)

gene_level_df$cell_type <- factor(gsub("_", " ", gene_level_df$cell_type),
                                   levels = gsub("_", " ", cell_types))
gene_level_df$individual <- factor(gene_level_df$individual, levels = individuals)

p <- ggplot(gene_level_df, aes(x = mean_expr, y = dispersion)) +
  geom_point(alpha = 0.12, size = 0.3, colour = "grey30") +
  geom_smooth(method = "loess", span = loess_span, colour = "#4575b4",
              linewidth = 0.4, se = TRUE) +
  facet_grid(individual ~ cell_type) +
  labs(x = "Log mean expression", y = "Dispersion",
       title = "Dispersion vs. mean (ribo-dropped, -6-excluded)",
       subtitle = "Rows: individual. Columns: cell type.") +
  theme_classic(base_family = "Arial") +
  theme(
    text = element_text(size = 6), axis.text = element_text(size = 6),
    axis.title = element_text(size = 6), plot.title = element_text(size = 7),
    plot.subtitle = element_text(size = 6), strip.text = element_text(size = 6),
    legend.text = element_text(size = 6), legend.title = element_text(size = 6),
    axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
    plot.margin = margin(2, 2, 2, 2, unit = "mm"))

out_file <- file.path(out_dir, "combined_loess_no_mt_highlight_grid.pdf")
ggsave(out_file, p, width = length(cell_types) * 1.8, height = length(individuals) * 1.8,
       units = "in", device = cairo_pdf, limitsize = FALSE)
message("[done] combined grid saved: ", out_file)
