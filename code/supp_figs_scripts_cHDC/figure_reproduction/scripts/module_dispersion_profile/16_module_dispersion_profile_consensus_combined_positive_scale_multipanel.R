#!/usr/bin/env Rscript

stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})

OUT_DIR      <- "figure_reproduction/outputs/module_dispersion_profile/phase2_consensus"
SUMMARY_CSV  <- file.path(OUT_DIR, "module_dispersion_summary_consensus.csv")
CATEGORY_CSV <- file.path(OUT_DIR, "module_category_lookup_consensus.csv")
SHIFT_CSV    <- file.path(OUT_DIR, "module_dispersion_shift_constants_consensus.csv")
for (f in c(SUMMARY_CSV, CATEGORY_CSV, SHIFT_CSV)) {
  if (!file.exists(f)) stop("Missing ", f, " -- run scripts 02, 09, and 14 first.")
}

CELL_TYPES <- c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                "Epithelial_cells", "Activated_fibroblasts",
                "Foregut_endodermal_cells")

PAGE_WIDTH_IN  <- 11
PAGE_HEIGHT_IN <- 8.5
GRID_NCOL <- 3
GRID_NROW <- 2

MULTIPANEL_SCALE <- 1 / 3

AXIS_TITLE_SIZE <- 13 * MULTIPANEL_SCALE * 4
PLOT_TITLE_SIZE <- 16 * MULTIPANEL_SCALE * 2
AXIS_TEXT_Y_SIZE <- 11 * MULTIPANEL_SCALE * 2
LEGEND_TEXT_SIZE <- 11 * MULTIPANEL_SCALE * 2

log_msg <- function(...) message("[16_module_dispersion_profile_consensus_combined_positive_scale_multipanel] ", ...)

OTHER_MODULE_COLOR    <- "grey30"
CC_MODULE_COLOR       <- "#009E73"
CELLTYPE_MARKER_COLOR <- "#D55E00"
CHAMBER_MARKER_COLOR  <- "#CC79A7"

summary_df  <- read.csv(SUMMARY_CSV, stringsAsFactors = FALSE)
category_df <- read.csv(CATEGORY_CSV, stringsAsFactors = FALSE)
shift_df    <- read.csv(SHIFT_CSV, stringsAsFactors = FALSE)
stopifnot(
  all(c("cell_type", "module", "module_size", "median_dispersion", "module_rank",
        "n_high_dispersion", "pct_high_dispersion") %in% colnames(summary_df)),
  all(c("cell_type", "module", "category") %in% colnames(category_df)),
  all(c("cell_type", "shift_constant") %in% colnames(shift_df)),
  setequal(unique(summary_df$cell_type), CELL_TYPES),
  setequal(unique(shift_df$cell_type), CELL_TYPES)
)

plot_df <- merge(summary_df, shift_df[, c("cell_type", "shift_constant")], by = "cell_type")
stopifnot("merge lost or duplicated rows" = nrow(plot_df) == nrow(summary_df))
plot_df$shifted_dispersion <- plot_df$median_dispersion + plot_df$shift_constant
stopifnot("shift left negative values -- investigate" = all(plot_df$shifted_dispersion >= 0))

plot_df$module_rank <- ave(plot_df$module_rank, plot_df$cell_type, FUN = function(x) rank(x, ties.method = "first"))
plot_df <- merge(plot_df, category_df[, c("cell_type", "module", "category")],
                  by = c("cell_type", "module"), all.x = TRUE)
plot_df$category[is.na(plot_df$category)] <- "other"
label_nudge <- 0.04 * diff(range(plot_df$shifted_dispersion, na.rm = TRUE))
plot_df$count_label_y <- plot_df$shifted_dispersion + label_nudge
plot_df$count_label <- as.character(plot_df$n_high_dispersion)
proportion_str <- sprintf("%.2f", plot_df$pct_high_dispersion / 100)
plot_df$pct_label <- ifelse(substr(proportion_str, 1, 2) == "0.", substring(proportion_str, 2), proportion_str)

LABEL_SIZE <- 1.7 * 1.5 * MULTIPANEL_SCALE

fill_values <- c(other = OTHER_MODULE_COLOR, cell_cycle_module = CC_MODULE_COLOR,
                  celltype_marker = CELLTYPE_MARKER_COLOR, chamber_marker = CHAMBER_MARKER_COLOR)
fill_labels <- c(other = "Other module", cell_cycle_module = "Cell cycle module (F-060, large S/G2M)",
                  celltype_marker = "Cell-type marker gene", chamber_marker = "Chamber-specificity marker gene")

panels <- list()
for (i in seq_along(CELL_TYPES)) {
  ct <- CELL_TYPES[i]
  ct_col <- ((i - 1) %% GRID_NCOL) + 1
  ct_label <- gsub("_", " ", ct)
  ct_df <- plot_df[plot_df$cell_type == ct, ]
  stopifnot("cell type missing from plot_df" = nrow(ct_df) > 0)
  ct_df$category <- factor(ct_df$category,
                            levels = c("other", "cell_cycle_module", "celltype_marker", "chamber_marker"))
  ct_range <- range(ct_df$shifted_dispersion, na.rm = TRUE)
  ct_df$pct_label_y <- ct_range[1] - 0.14 * diff(ct_range)

  panels[[ct]] <- ggplot(ct_df, aes(x = module_rank, y = shifted_dispersion, fill = category)) +
    geom_col(width = 0.7, color = "grey20", linewidth = 0.15 * MULTIPANEL_SCALE) +
    geom_text(aes(y = count_label_y, label = count_label), size = LABEL_SIZE, family = "Arial") +
    geom_text(aes(y = pct_label_y, label = pct_label), size = LABEL_SIZE, family = "Arial",
              angle = 90, hjust = 0) +
    scale_fill_manual(values = fill_values, labels = fill_labels, name = NULL, drop = FALSE) +
    labs(x = NULL, y = NULL, title = ct_label) +
    theme_classic(base_family = "Arial") +
    theme(
      text = element_text(size = 13 * MULTIPANEL_SCALE),
      axis.text.y = if (ct_col == 1) element_text(size = AXIS_TEXT_Y_SIZE) else element_blank(),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      plot.title = element_text(size = PLOT_TITLE_SIZE, face = "bold"),
      legend.position = "none",
      axis.line = element_line(linewidth = 0.3 * MULTIPANEL_SCALE), axis.ticks = element_line(linewidth = 0.3 * MULTIPANEL_SCALE),
      plot.margin = margin(1, 1, 1, 1, unit = "mm"),
      panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.3 * MULTIPANEL_SCALE)
    )
}

legend_source <- panels[[CELL_TYPES[1]]] + theme(legend.position = "bottom", legend.text = element_text(size = LEGEND_TEXT_SIZE))
shared_legend <- get_legend(legend_source)

x_title <- ggdraw() + draw_label("Module rank (ascending median dispersion)", fontfamily = "Arial", size = AXIS_TITLE_SIZE)
y_title <- ggdraw() + draw_label("Median dispersion", fontfamily = "Arial", size = AXIS_TITLE_SIZE, angle = 90)

grid <- plot_grid(plotlist = panels, ncol = GRID_NCOL, nrow = GRID_NROW)
grid_with_y_title <- plot_grid(y_title, grid, ncol = 2, rel_widths = c(0.05, 1))
grid_with_axis_titles <- plot_grid(grid_with_y_title, x_title, ncol = 1, rel_heights = c(1, 0.06))
combined <- plot_grid(grid_with_axis_titles, shared_legend, ncol = 1, rel_heights = c(1, 0.06))

out_pdf <- file.path(OUT_DIR, "module_dispersion_profile_consensus_combined_positive_scale_multipanel.pdf")
ggsave(out_pdf, combined, device = cairo_pdf, width = PAGE_WIDTH_IN, height = PAGE_HEIGHT_IN, units = "in")
log_msg("Wrote: ", out_pdf, " (1 page, ", length(CELL_TYPES), " panels, ",
        GRID_NCOL, "x", GRID_NROW, " grid, ", PAGE_WIDTH_IN, "x", PAGE_HEIGHT_IN, "in)")

log_msg("DONE.")
