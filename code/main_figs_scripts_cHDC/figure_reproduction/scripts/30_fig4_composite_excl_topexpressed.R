#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggpubr)  # saved panels use stat_compare_means(), evaluated lazily at print/ggsave time
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
analysis_dir <- file.path(repo_dir, "figure_reproduction")
fig4_dir <- file.path(analysis_dir, "outputs", "fig4_replication")

load_panel <- function(name) {
  path <- file.path(fig4_dir, name)
  if (!file.exists(path)) stop("Missing ", path, " -- run 29_fig4_panels_excl_topexpressed.R first")
  readRDS(path)
}

p_a <- load_panel("fig4a_panel.rds")
p_b <- load_panel("fig4b_panel.rds")
p_c_tss <- load_panel("fig4c_panel.rds")
p_c_enh <- load_panel("fig4d_panel.rds")
p_c_intron <- load_panel("fig4e_panel.rds")

# ------------------------------- P-VALUE THRESHOLD LABELS -----------------------
# Identical convention to script 28's Fig3 composite: never report an exact
# p-value -- collapse to a fixed significance-threshold ladder, with a flat
# "P < 10^-15" floor for anything smaller (not the value's own exponent).
format_pvalue_threshold <- function(p) {
  if (p >= 0.05) return("P >= 0.05")  # plotmath-safe (no bare "=")
  if (p >= 0.01) return("P < 0.05")
  if (p >= 0.001) return("P < 0.01")
  if (p >= 1e-15) return("P < 0.001")
  "P < 10^-15"
}

drop_layers_by_class <- function(p, slot, class_name) {
  keep <- !vapply(p$layers, function(l) inherits(l[[slot]], class_name), logical(1))
  p$layers <- p$layers[keep]
  p
}

# geom_text/annotate's `size` is in mm, not points (1/.pt converts pt -> mm).
# Same 1pt bump over the original size=2 (~6pt) as script 28's Fig3 panels.
pvalue_text_size <- 2 + 1 / .pt

# Fig4A and Fig4C's 3 sub-panels all use ggpubr::stat_compare_means() for
# their p-value text -- remove it, add a manual threshold-formatted label at
# the same top-center position (x = 1.5, midpoint of the 2 factor levels).
add_pvalue_bucket <- function(p, y_var, group_var) {
  d <- p$data
  p_val <- wilcox.test(d[[y_var]] ~ d[[group_var]])$p.value
  p <- drop_layers_by_class(p, "stat", "StatCompareMeans")
  p + annotate("text", x = 1.5, y = Inf, label = format_pvalue_threshold(p_val),
               vjust = 1.5, size = pvalue_text_size, family = "Arial", parse = TRUE)
}

p_a <- add_pvalue_bucket(p_a, "DE_count", "rank")
p_c_tss <- add_pvalue_bucket(p_c_tss, "tss_count", "rank")
p_c_enh <- add_pvalue_bucket(p_c_enh, "enh_length", "rank")
p_c_intron <- add_pvalue_bucket(p_c_intron, "intron_count", "rank")

# Fig4B's 2 p-values (low-dispersion Fisher test, high-dispersion Fisher
# test) are baked into its own data (a `p_label` column feeding a
# geom_text() layer), like Fig3A -- recompute that column with the
# threshold formatter, then rebuild the geom_text layer with parse = TRUE.
pb_data <- p_b$data
pb_data$p_label <- vapply(pb_data$p_val, format_pvalue_threshold, character(1))
p_b <- drop_layers_by_class(p_b, "geom", "GeomText")
p_b$data <- pb_data
# y positioned at Inf (not the original data-based log2_ci_upper + vjust=-0.6)
# -- the high-dispersion errorbar already reaches near the panel's top, and
# at the bumped font size the old data-based position pushed the label
# outside the plot area (clipped). Inf + vjust anchors it to the panel edge
# instead, immune to font size.
p_b <- p_b + geom_text(aes(x = gene_set, label = p_label), y = Inf, vjust = 1.5,
                        size = pvalue_text_size, family = "Arial", parse = TRUE)

# ------------------------------- FONT SIZE BUMPS ---------------------------------
# Identical bump to script 28's Fig3 composite: axis tick labels 6->7pt,
# y-axis titles 6->8pt. None of Fig4's panels have an x-axis title or a
# sub-panel plot.title (unlike Fig3A/Fig3C) -- those bumps have no target
# here, so they're omitted.
axis_text_bump <- theme(axis.text = element_text(size = 7))
y_title_bump <- theme(axis.title.y = element_text(size = 8))

p_a <- p_a + axis_text_bump + y_title_bump
p_b <- p_b + axis_text_bump + y_title_bump
p_c_tss <- p_c_tss + axis_text_bump + y_title_bump
p_c_enh <- p_c_enh + axis_text_bump + y_title_bump
p_c_intron <- p_c_intron + axis_text_bump + y_title_bump

# Published "C" is a single 3-panel row (TSS | enhancer length | intron count).
p_c_combined <- p_c_tss + p_c_enh + p_c_intron + plot_layout(ncol = 3)
p_c_unit <- wrap_elements(full = p_c_combined)

# Tag-letter size: A/B's tags inherit each panel's own small `text` theme
# (cascading to plot.tag); C's wrap_elements()-wrapped tag instead falls
# back to patchwork's own larger default -- same mismatch as Fig3, same fix
# (the `&` operator forces one plot.tag size onto every constituent plot).
fig4_composite <- (p_a + p_b) / p_c_unit +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

out_pdf <- file.path(fig4_dir, "Figure4_composite_excl_topexpressed.pdf")
# width matches script 28's Fig3 composite exactly.
ggsave(out_pdf, fig4_composite, device = cairo_pdf, width = 7.2, height = 6.4, units = "in")
message("[done] wrote ", out_pdf)
