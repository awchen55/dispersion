#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(variancePartition)
  library(ggplot2)
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
analysis_dir <- file.path(repo_dir, "figure_reproduction")
fig2_dir <- file.path(analysis_dir, "outputs", "fig2_replication")

vp_csv <- file.path(fig2_dir, "var_partition_retained.csv")
if (!file.exists(vp_csv)) stop("Missing ", vp_csv, " -- run 25_fig2b_varpart_excl_topexpressed.R first")
vp_res_var <- read.csv(vp_csv, row.names = 1, check.names = FALSE)

# General rule: strip underscores (replace with a space). "perc_mt" gets a
# full rewrite instead -- "% mito"/"chondrial reads" on 2 lines, splitting
# "mitochondrial" mid-word at "mito". "Residuals" is left untouched (no
# underscore, and sortCols() name-matches it for its always-last behavior).
label_map <- c(
  cell_line = "cell line",
  cell_type = "cell type",
  cell_count = "cell count",
  nFeatures = "nFeatures",
  perc_mt = "% mito\nchondrial reads",
  Residuals = "Residuals"
)
stopifnot("label_map must cover every column in vp_res_var exactly" =
            setequal(names(label_map), colnames(vp_res_var)))
colnames(vp_res_var) <- unname(label_map[colnames(vp_res_var)])

vp_plot_res_var <- plotVarPart(sortCols(vp_res_var), ylab = "Variance explained (%)") +
  theme_classic(base_family = "Arial") +
  theme(
    text = element_text(size = 6),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    legend.title = element_blank(),
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

out_pdf <- file.path(fig2_dir, "Fig2b_VarPart_ResVar_relabeled.pdf")
ggsave(
  filename = out_pdf,
  plot = vp_plot_res_var,
  device = cairo_pdf,
  width = 3,
  height = 2.5,
  units = "in"
)
message("[done] wrote ", out_pdf)
