#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
fig_repro_dir <- file.path(repo_dir, "figure_reproduction")
out_dir <- file.path(fig_repro_dir, "outputs", "mean_matched_regulatory_revision", "combined")
panel_rds <- file.path(out_dir, "panel_data.rds")
if (!file.exists(panel_rds)) stop("Missing ", panel_rds, " -- run 08_remaining_panels_matched.R first")
pd <- readRDS(panel_rds)

source(file.path(fig_repro_dir, "scripts", "utils.R"))  # make_low_high_boxplot()

# Same ladder as rnx 02 / 28 / 30: never print an exact p-value.
P_LADDER <- c(0.05, 0.01, 0.001, 1e-15)
# ANALYSIS_OK[within-file-helper]: called below in this script only; local by design
format_pvalue_threshold <- function(p) {
  if (p >= P_LADDER[1]) return("P >= 0.05")
  if (p >= P_LADDER[2]) return("P < 0.05")
  if (p >= P_LADDER[3]) return("P < 0.01")
  if (p >= P_LADDER[4]) return("P < 0.001")
  "P < 10^-15"
}
pvalue_text_size <- 2 + 1 / .pt  # rnx 02: 1pt above the panels' original 6pt
composite_theme <- theme(axis.text = element_text(size = 7),
                         axis.text.x = element_text(size = 5.5),
                         axis.title.y = element_text(size = 8),
                         plot.title = element_text(size = 7, hjust = 0.5))
two_line_labels <- scale_x_discrete(labels = c("1" = "Low-\ndispersion", "5" = "High-\ndispersion"))

# Low-vs-high boxplot with the Wilcoxon p-value as a threshold label (rnx 02).
# ANALYSIS_OK[within-file-helper]: called below in this script only; local by design
boxplot_panel <- function(df, y_var, y_label, title) {
  p_val <- wilcox.test(df[[y_var]] ~ df$rank)$p.value
  p <- make_low_high_boxplot(df, y_var, y_label)
  # ANALYSIS_OK[plot-layer-drop]: removes ggplot layers (ggpubr p-value text), not data rows; replaced by the annotate() below
  p$layers <- p$layers[!vapply(p$layers, function(l) inherits(l$stat, "StatCompareMeans"), logical(1))]
  p + annotate("text", x = 1.5, y = Inf, label = format_pvalue_threshold(p_val), vjust = 1.5,
               size = pvalue_text_size, family = "Arial", parse = TRUE) +
    two_line_labels + labs(title = title) + composite_theme +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))  # headroom so the P label clears the points
}

p_b <- boxplot_panel(pd$fig3b, "hub_rank", "Shared hub genes ranking", "TWN shared hubs")
p_c <- boxplot_panel(pd$fig3d, "int_count", "Protein connections per gene", "PPI connectivity")
p_d <- boxplot_panel(pd$fig3c_fib, "n_genes", "Connections per gene", "TWN, Cardiac fibroblasts")
p_e <- boxplot_panel(pd$fig3c_card, "n_genes", "Connections per gene", "TWN, Cardiomyocytes")
p_f <- boxplot_panel(pd$fig4a, "DE_count", "DE count per gene", "Differential expression")

# --- G: TATA Fisher enrichment (rnx 02 panel G styling; 05's matched TATA table)
pg <- pd$fig4b
pg$p_label <- vapply(pg$p_value, format_pvalue_threshold, character(1))
p_g <- ggplot(pg, aes(x = gene_set, y = log2_or)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey40") +
  geom_errorbar(aes(ymin = log2_ci_lower, ymax = log2_ci_upper), width = 0.15, linewidth = 0.4) +
  geom_point(size = 1.6) +
  geom_text(aes(label = p_label), y = Inf, vjust = 1.5, size = pvalue_text_size,
            family = "Arial", parse = TRUE) +
  scale_x_discrete(labels = c("Low-dispersion" = "Low-\ndispersion", "High-dispersion" = "High-\ndispersion")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = NULL, y = expression(log[2]("Odds ratio")), title = "TATA-box promoters") +
  theme_classic(base_family = "Arial") +
  theme(text = element_text(size = 6), axis.line = element_line(linewidth = 0.4),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        plot.margin = margin(2, 2, 2, 2, unit = "mm")) +
  composite_theme

p_h <- boxplot_panel(pd$fig4c, "tss_count", "TSS count", "TSS count")
p_i <- boxplot_panel(pd$fig4d, "enh_length", "Accum. enhancer length\n(per tissue/cell type)",
                     "Enhancer length")
p_j <- boxplot_panel(pd$fig4e, "intron_count", "Intron count", "Intron count")

composite <- p_b + p_c + p_d + p_e + p_f + p_g + p_h + p_i + p_j + plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "bold", family = "Arial"))

out_pdf <- file.path(out_dir, "regulatory_network_composite_matched.pdf")
ggsave(out_pdf, composite, device = cairo_pdf, width = 7.2, height = 6.4, units = "in")
ggsave(sub("\\.pdf$", ".png", out_pdf), composite, device = ragg::agg_png,
       width = 7.2, height = 6.4, units = "in", dpi = 300, bg = "white")
message("[done] wrote ", out_pdf)
