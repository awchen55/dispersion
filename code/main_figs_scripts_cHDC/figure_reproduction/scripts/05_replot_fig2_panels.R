#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
fig_repro_dir <- file.path(repo_dir, "figure_reproduction")
out_dir <- file.path(fig_repro_dir, "outputs")
fig2_dir <- file.path(out_dir, "fig2_replication")

dd_sig_thresh <- 0.05  # unused directly here, kept for parity with 04's csv contents
n_top_terms_part3 <- 15  # matches 04's dotplot(..., showCategory = 15)

# =============================================================================
# PART 1 — Fig 2D: DE (mean) vs DD (dispersion) scatter, contrast K
# =============================================================================

# Byte-for-byte copy of 04's make_de_dd_scatter().
make_de_dd_scatter <- function(plot_df) {
  ggplot(plot_df, aes(x = de_logFC, y = dv_logFC, color = dv_sig)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey40") +
    geom_point(size = 0.8, alpha = 0.8) +
    scale_color_manual(
      values = c("FALSE" = "grey70", "TRUE" = "#B22222"),
      labels = c("FALSE" = "Non DD", "TRUE" = "DD")
    ) +
    labs(x = "DE log2 fold change", y = "DD log2 fold change", color = NULL) +
    theme_classic(base_family = "Arial") +
    theme(
      text = element_text(size = 6),
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 6),
      axis.line = element_line(linewidth = 0.4),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
      legend.position = "right",
      legend.text = element_text(size = 6),
      legend.key.size = unit(6, "pt"),
      plot.margin = margin(2, 2, 2, 2, unit = "mm")
    )
}

for (tag in c("full", "retained")) {
  csv_path <- file.path(fig2_dir, sprintf("fig2d_scatter_%s.csv", tag))
  if (!file.exists(csv_path)) stop("Missing ", csv_path, " -- run 04_replicate_fig2_panels.R first")
  df <- read.csv(csv_path)
  df$dv_sig <- factor(df$dv_sig, levels = c(FALSE, TRUE))
  p <- make_de_dd_scatter(df)
  ggsave(file.path(fig2_dir, sprintf("fig2d_scatter_%s.pdf", tag)), p,
         width = 3.5, height = 2, device = cairo_pdf)
  message(sprintf("[fig2d:%s] re-plotted from %s", tag, basename(csv_path)))
}

# =============================================================================
# PART 2 — Fig 2A: GO enrichment of consistently-low-dispersion genes
# =============================================================================

# Shared GO-dotplot style (Fig2A's own style; Part 3/Fig2E restyled to match
# it -- same theme_classic(base_family="Arial") 6pt convention as every other
# Fig2 panel).
#
# aspect.ratio = 2 REMOVED (matches source): the fixed ratio fights a variable
# label/legend budget across different GO-term lists.
#
# width=9 for both Part 2 (8 terms) and Part 3 (15 terms, longer labels):
# width=7 reproducibly collapses the panel in cairo_pdf, width=9 reproducibly
# renders every point/label cleanly (see source script's header for the
# verification detail).
make_go_dotplot <- function(plot_df, title) {
  ggplot(plot_df, aes(x = GeneRatio_num, y = Description)) +
    geom_point(aes(size = Count, color = log10_padj)) +
    scale_color_gradient(low = "grey85", high = "#8B0000",
                         name = "-log10(adjusted p)") +
    scale_size(range = c(2, 6)) +
    labs(title = title, x = "Gene ratio", y = NULL) +
    theme_classic(base_family = "Arial") +
    theme(
      text = element_text(size = 6),
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 6),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      plot.title = element_text(size = 6, hjust = 0.5),
      legend.key.height = unit(3, "mm"),
      legend.key.width  = unit(3, "mm"),
      plot.margin = margin(2, 2, 2, 2, unit = "mm")
    )
}

top8_csv <- file.path(fig2_dir, "fig2a_go_enrichment_retained_top8.csv")
if (!file.exists(top8_csv)) stop("Missing ", top8_csv, " -- run 04_replicate_fig2_panels.R first")
plot_df <- read.csv(top8_csv)
plot_df$Description <- factor(plot_df$ID, levels = rev(plot_df$ID))

p_go <- make_go_dotplot(plot_df, "GO enrichment: consistently low-dispersed genes")

ggsave(file.path(fig2_dir, "fig2a_go_enrichment_retained.pdf"), p_go,
       width = 9, height = 4, device = cairo_pdf)
message("[fig2a] re-plotted from ", basename(top8_csv))

# =============================================================================
# PART 3 — Fig 2E: Cardiomyocytes vs Foregut_endodermal_cells GO enrichment
# =============================================================================

for (tag in c("full", "retained")) {
  csv_path <- file.path(fig2_dir, sprintf("fig2_part3_card_vs_foregut_%s.csv", tag))
  if (!file.exists(csv_path)) stop("Missing ", csv_path, " -- run 04_replicate_fig2_panels.R first")
  res_df <- read.csv(csv_path)
  if (nrow(res_df) == 0) {
    message(sprintf("[part3:%s] no enriched terms in saved csv, skipping plot", tag))
    next
  }
  # ANALYSIS_OK[fixed-width-split]: GeneRatio is always clusterProfiler's
  # "k/n" string; str_split_fixed(..., 2) is guaranteed exactly 2 columns.
  res_df <- res_df %>%
    mutate(GeneRatio_num = as.numeric(str_split_fixed(GeneRatio, "/", 2)[, 1]) /
                            as.numeric(str_split_fixed(GeneRatio, "/", 2)[, 2]),
           log10_padj = -log10(p.adjust))
  # ANALYSIS_OK[plot-filter]: head(n_top_terms_part3) after ordering by
  # GeneRatio descending mirrors enrichplot::dotplot(orderBy="x")'s default
  # term selection/ordering; the saved csv already holds every term passing
  # the p.adjust/Count filter from 04, so this filter only affects what's
  # plotted, not the underlying enrichment result.
  plot_df <- res_df %>%
    arrange(desc(GeneRatio_num)) %>%
    head(n_top_terms_part3) %>%
    mutate(Description = factor(Description, levels = rev(Description)))

  p <- make_go_dotplot(plot_df, str_wrap(
    sprintf("Mean increased, variability decreased (%s)", tag), width = 30))

  ggsave(file.path(fig2_dir, sprintf("fig2_part3_card_vs_foregut_%s.pdf", tag)),
         p, width = 9, height = 4, device = cairo_pdf)
  message(sprintf("[part3:%s] re-plotted from %s", tag, basename(csv_path)))
}

message("[done] Fig 2 panels re-plotted in ", fig2_dir)
