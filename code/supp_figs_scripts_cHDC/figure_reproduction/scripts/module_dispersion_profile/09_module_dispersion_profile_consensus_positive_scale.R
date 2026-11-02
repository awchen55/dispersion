#!/usr/bin/env Rscript

stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

OUT_DIR <- "figure_reproduction/outputs/module_dispersion_profile/phase2_consensus"
SUMMARY_CSV <- file.path(OUT_DIR, "module_dispersion_summary_consensus.csv")
if (!file.exists(SUMMARY_CSV)) stop("Missing ", SUMMARY_CSV, " -- run script 02 first.")

CELL_TYPES <- c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                "Epithelial_cells", "Activated_fibroblasts",
                "Foregut_endodermal_cells")

celltype_colors <- c(
  "Activated fibroblasts"    = "#117733",
  "Cardiac fibroblasts"      = "#4477AA",
  "Cardiomyocytes"           = "#882255",
  "Epithelial cells"         = "#AA4499",
  "Foregut endodermal cells" = "#CC6677",
  "Mesenchymal cells"        = "#DDCC77"
)

log_msg <- function(...) message("[09_module_dispersion_profile_consensus_positive_scale] ", ...)

summary_df <- read.csv(SUMMARY_CSV, stringsAsFactors = FALSE)
stopifnot(
  "module_dispersion_summary_consensus.csv missing expected columns" =
    all(c("cell_type", "module", "module_size", "median_dispersion", "n_high_dispersion",
          "pct_high_dispersion", "module_rank") %in% colnames(summary_df))
)
stopifnot("expected all 6 cell types in summary_df" = setequal(unique(summary_df$cell_type), CELL_TYPES))

shift_rows <- list()
plot_df <- summary_df
plot_df$shifted_dispersion <- NA_real_
for (ct in CELL_TYPES) {
  ct_idx <- plot_df$cell_type == ct
  min_val <- min(plot_df$median_dispersion[ct_idx])
  plot_df$shifted_dispersion[ct_idx] <- plot_df$median_dispersion[ct_idx] - min_val
  log_msg(ct, ": shift constant = ", round(-min_val, 4),
          " (min median_dispersion ", round(min_val, 4), " -> 0)")
  shift_rows[[ct]] <- data.frame(cell_type = ct, min_median_dispersion = min_val,
                                  shift_constant = -min_val, stringsAsFactors = FALSE)
}
stopifnot("shift left negative values -- investigate" = all(plot_df$shifted_dispersion >= 0))

shift_df <- do.call(rbind, shift_rows); rownames(shift_df) <- NULL
shift_csv <- file.path(OUT_DIR, "module_dispersion_shift_constants_consensus.csv")
write.csv(shift_df, shift_csv, row.names = FALSE)
log_msg("Wrote: ", shift_csv)

plot_df$module_rank <- ave(plot_df$module_rank, plot_df$cell_type, FUN = function(x) rank(x, ties.method = "first"))
label_nudge <- 0.04 * diff(range(plot_df$shifted_dispersion, na.rm = TRUE))
plot_df$label_y <- plot_df$shifted_dispersion + label_nudge
plot_df$label_count <- as.character(plot_df$n_high_dispersion)
plot_df$label_pct <- sprintf("%.0f%%", plot_df$pct_high_dispersion)

LABEL_SIZE <- 1.7 * 1.5

make_positive_scale_panel <- function(ct_df, fill_color, label_col, panel_title) {
  ggplot(ct_df, aes(x = module_rank, y = shifted_dispersion)) +
    geom_col(width = 0.7, fill = fill_color, color = "grey20", linewidth = 0.15) +
    geom_text(aes(y = label_y, label = .data[[label_col]]), size = LABEL_SIZE, family = "Arial") +
    labs(x = "Module rank (ascending median dispersion)",
         y = "Median dispersion, shifted to a positive scale\n(per-cell-type min subtracted; relative differences unchanged)",
         title = panel_title) +
    theme_classic(base_family = "Arial") +
    theme(
      text = element_text(size = 13), axis.text = element_text(size = 11),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.title = element_text(size = 13),
      plot.title = element_text(size = 14, face = "bold"),
      axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
      plot.margin = margin(4, 4, 4, 4, unit = "mm"),
      panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.3)
    )
}

out_pdf <- file.path(OUT_DIR, "module_dispersion_profile_consensus_positive_scale.pdf")
cairo_pdf(out_pdf, width = 11, height = 15, onefile = TRUE)
for (ct in CELL_TYPES) {
  ct_label <- gsub("_", " ", ct)
  ct_df <- plot_df[plot_df$cell_type == ct, ]
  stopifnot("cell type missing from plot_df" = nrow(ct_df) > 0)
  fill_color <- celltype_colors[[ct_label]]

  p_count <- make_positive_scale_panel(ct_df, fill_color, "label_count",
                                        "Raw count of top-quintile-dispersion genes per module")
  p_pct   <- make_positive_scale_panel(ct_df, fill_color, "label_pct",
                                        "Proportion of top-quintile-dispersion genes per module")
  gridExtra::grid.arrange(
    p_count, p_pct, nrow = 2,
    top = grid::textGrob(
      paste0(ct_label, " -- consensus module dispersion profile (positive scale)"),
      gp = grid::gpar(fontsize = 17, fontface = "bold", fontfamily = "Arial")
    )
  )
}
dev.off()
log_msg("Wrote: ", out_pdf, " (", length(CELL_TYPES), " pages, one per cell type, count+proportion stacked)")

log_msg("DONE.")
