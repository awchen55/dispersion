#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ---- Named constants ---------------------------------------------------------
repo_dir     <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
analysis_dir <- file.path(repo_dir, "analysis/subcluster_dispersion_robustness")
config_label <- "mincells700_minind2"  # Brendan's explicit scoping, 2026-09-02
out_dir      <- file.path(analysis_dir, "outputs", config_label)

resolutions <- c("0.5", "1.0", "2.0")

# Matches script 03's heatmap_rho_range/hm_col_fun exactly (the shared scale
# Brendan noticed) -- NOT re-derived per panel here.
heatmap_rho_range <- c(0.8, 1)
hm_col_fun <- circlize::colorRamp2(
  c(heatmap_rho_range[1], mean(heatmap_rho_range), heatmap_rho_range[2]),
  c("blue", "white", "red")
)

source(file.path(analysis_dir, "scripts", "utils.R"))
log_msg <- make_logger("05_combined_heatmaps")

ht_list <- lapply(resolutions, function(res_label) {
  cor_csv <- file.path(out_dir, sprintf("across_subcluster_cor_res_%s.csv", res_label))
  if (!file.exists(cor_csv)) {
    stop(sprintf("Missing %s -- run 03_robustness_correlations.R for %s first.",
                  cor_csv, config_label))
  }
  cor_mat <- as.matrix(read.csv(cor_csv, row.names = 1, check.names = FALSE))
  log_msg("[res %s] %d x %d correlation matrix.", res_label, nrow(cor_mat), ncol(cor_mat))

  # Legend only on the last (rightmost) panel -- all three share hm_col_fun,
  # so one legend covers all of them; the other two would just duplicate it.
  ComplexHeatmap::Heatmap(
    cor_mat,
    name = "Spearman rho",
    col = hm_col_fun,
    column_title = sprintf("res %s (%d sub-clusters)", res_label, ncol(cor_mat)),
    row_names_gp = grid::gpar(fontsize = 6),
    column_names_gp = grid::gpar(fontsize = 6),
    show_heatmap_legend = (res_label == resolutions[[length(resolutions)]])
  )
})

out_file <- file.path(out_dir, "across_subcluster_cor_heatmap_combined.pdf")
pdf(out_file, width = 18, height = 6)
grid::grid.newpage()
grid::pushViewport(grid::viewport(
  layout = grid::grid.layout(nrow = 1, ncol = length(resolutions))))
for (i in seq_along(ht_list)) {
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = i))
  ComplexHeatmap::draw(ht_list[[i]], newpage = FALSE)
  grid::popViewport()
}
grid::popViewport()
dev.off()
log_msg("Saved: %s", out_file)
