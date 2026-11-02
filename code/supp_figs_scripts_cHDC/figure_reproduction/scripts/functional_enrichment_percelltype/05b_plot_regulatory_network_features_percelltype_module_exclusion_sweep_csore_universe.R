#!/usr/bin/env Rscript

stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

OUT_DIR <- "figure_reproduction/outputs/functional_enrichment_percelltype"
RESULTS_CSV <- file.path(OUT_DIR, "regulatory_network_features_percelltype_module_exclusion_sweep_csore_universe_results.csv")
TATA_SWEEP_CSV <- file.path(OUT_DIR, "tata_enrichment_module_exclusion_sweep_csore_universe_results.csv")

for (f in c(RESULTS_CSV, TATA_SWEEP_CSV)) {
  if (!file.exists(f)) stop("Missing required input: ", f, " -- run script 05 first for the features CSV, script 04 first for the TATA sweep CSV.")
}

CELL_TYPES <- c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                "Epithelial_cells", "Activated_fibroblasts",
                "Foregut_endodermal_cells")
EXCLUSION_FRACTIONS <- seq(0.1, 0.7, by = 0.1)

log_msg <- function(...) message("[05b_plot_regulatory_network_features_percelltype_module_exclusion_sweep_csore_universe] ", ...)

celltype_colors <- c(
  "Activated fibroblasts"    = "#117733", "Cardiac fibroblasts" = "#4477AA",
  "Cardiomyocytes"           = "#882255", "Epithelial cells"    = "#AA4499",
  "Foregut endodermal cells" = "#CC6677", "Mesenchymal cells"   = "#DDCC77"
)

FEATURE_LABELS <- c(
  tss_count = "TSS count", enhancer_length = "Enhancer length (bp)",
  ppi_connectivity = "PPI connections", twn_hub_rank = "TWN shared-hub rank",
  intron_count = "Intron count"
)

results_df <- read.csv(RESULTS_CSV, stringsAsFactors = FALSE)
stopifnot("results CSV missing expected columns" =
            all(c("cell_type", "feature", "exclusion_fraction", "hl_shift", "ci_lower", "ci_upper") %in% colnames(results_df)))
plot_df <- results_df
plot_df$cell_type_label <- factor(gsub("_", " ", plot_df$cell_type), levels = gsub("_", " ", CELL_TYPES))

tata_df <- read.csv(TATA_SWEEP_CSV, stringsAsFactors = FALSE)
stopifnot("tata_enrichment_module_exclusion_sweep_csore_universe_results.csv missing expected columns" =
            all(c("cell_type", "exclusion_fraction", "odds_ratio", "ci_lower", "ci_upper") %in% colnames(tata_df)))
tata_df$cell_type_label <- factor(gsub("_", " ", tata_df$cell_type), levels = gsub("_", " ", CELL_TYPES))
tata_df$log2_or <- log2(tata_df$odds_ratio)
tata_df$log2_ci_lower <- log2(tata_df$ci_lower)
tata_df$log2_ci_upper <- log2(tata_df$ci_upper)

FONT_SCALE <- 1.5
X_LAB <- "Fraction of highest-dispersion consensus modules excluded"
HL_Y_LAB <- "Hodges-Lehmann shift"

make_panel <- function(df, y_col, ymin_col, ymax_col, panel_title, y_lab = "", show_x_text = TRUE) {
  p <- ggplot(df, aes(x = exclusion_fraction, y = .data[[y_col]], color = cell_type_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey40") +
    geom_errorbar(aes(ymin = .data[[ymin_col]], ymax = .data[[ymax_col]]), width = 0.01, linewidth = 0.3, alpha = 0.5) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.4) +
    scale_color_manual(values = celltype_colors, name = NULL) +
    scale_x_continuous(breaks = c(0, EXCLUSION_FRACTIONS), labels = scales::percent) +
    labs(x = NULL, y = y_lab, title = panel_title) +
    theme_classic(base_family = "Arial") +
    theme(
      text = element_text(size = 7 * FONT_SCALE), axis.title = element_text(size = 7 * FONT_SCALE),
      axis.text = element_text(size = 6), axis.line = element_line(linewidth = 0.4),
      legend.text = element_text(size = 7), legend.title = element_text(size = 7),
      legend.key.size = unit(7, "pt"),
      plot.title = element_text(size = 7 * FONT_SCALE, face = "bold"),
      plot.margin = margin(2, 2, 2, 2, unit = "mm"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4)
    )
  if (!show_x_text) p <- p + theme(axis.text.x = element_blank())
  p
}

TOP_ROW <- c("tss_count", "enhancer_length", "ppi_connectivity")
feature_panels <- lapply(names(FEATURE_LABELS), function(feat_name) {
  make_panel(plot_df[plot_df$feature == feat_name, ], "hl_shift", "ci_lower", "ci_upper",
             panel_title = FEATURE_LABELS[feat_name], show_x_text = !(feat_name %in% TOP_ROW))
})
names(feature_panels) <- names(FEATURE_LABELS)

tata_panel <- make_panel(tata_df, "log2_or", "log2_ci_lower", "log2_ci_upper",
                          panel_title = "TATA box enrichment",
                          y_lab = expression(log[2]("Odds ratio")))

y_label_grob <- wrap_elements(grid::textGrob(HL_Y_LAB, rot = 90, gp = grid::gpar(fontsize = 7 * FONT_SCALE, fontfamily = "Arial")))

panel_grid <- list(
  L = y_label_grob,
  A = feature_panels[["tss_count"]], B = feature_panels[["enhancer_length"]], C = feature_panels[["ppi_connectivity"]],
  D = feature_panels[["twn_hub_rank"]], E = feature_panels[["intron_count"]], G = tata_panel
)
grid_design <- "
LABC
LDEG
"
combined <- wrap_plots(panel_grid, design = grid_design, guides = "collect") +
  plot_layout(widths = c(0.05, 1, 1, 1)) +
  plot_annotation(
    title = "Regulatory/network feature persistence vs. top-x% consensus module exclusion",
    caption = X_LAB,
    theme = theme(
      plot.title = element_text(size = 8 * FONT_SCALE, face = "bold", family = "Arial"),
      plot.caption = element_text(size = 7 * FONT_SCALE, family = "Arial", hjust = 0.5)
    )
  ) &
  theme(legend.position = "right")

out_pdf <- file.path(OUT_DIR, "regulatory_network_features_percelltype_module_exclusion_sweep_csore_universe.pdf")
ggsave(out_pdf, combined, width = 9, height = 7, units = "in", device = cairo_pdf)
log_msg("Wrote: ", out_pdf, " (6 panels: 5 regulatory/network features + TATA box enrichment, bottom-right)")

log_msg("DONE.")
