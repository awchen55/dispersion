#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggpubr)  # saved panels use stat_compare_means(), evaluated lazily at print/ggsave time
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
analysis_dir <- file.path(repo_dir, "figure_reproduction")
fig3_dir <- file.path(analysis_dir, "outputs", "fig3_replication")

load_panel <- function(name) {
  path <- file.path(fig3_dir, name)
  if (!file.exists(path)) stop("Missing ", path, " -- run 27_fig3_panels_excl_topexpressed.R first")
  readRDS(path)
}

p_a <- load_panel("fig3a_panel.rds")
p_b <- load_panel("fig3b_panel.rds")
p_c_fib <- load_panel("fig3c_fib_panel.rds")
p_c_card <- load_panel("fig3c_card_panel.rds")
p_d <- load_panel("fig3d_panel.rds")

# ------------------------------- P-VALUE THRESHOLD LABELS -----------------------
# Convention: never report an exact p-value on a panel -- collapse to a fixed
# significance-threshold ladder instead: P < 0.05, P < 0.01, P < 0.001, and
# P < 10^-15 as the floor for anything smaller (a fixed cap, NOT the value's
# own order of magnitude).
format_pvalue_threshold <- function(p) {
  if (p >= 0.05) return("P >= 0.05")  # not hit by any current Fig3 panel; plotmath-safe (no bare "=")
  if (p >= 0.01) return("P < 0.05")
  if (p >= 0.001) return("P < 0.01")
  if (p >= 1e-15) return("P < 0.001")
  "P < 10^-15"
}

# Drop all layers of a given stat/geom subclass (e.g. ggpubr's dynamically-
# computed stat_compare_means text) so a manually-labeled replacement can be
# added in its place.
drop_layers_by_class <- function(p, slot, class_name) {
  keep <- !vapply(p$layers, function(l) inherits(l[[slot]], class_name), logical(1))
  p$layers <- p$layers[keep]
  p
}

# Fig3B/C/D's boxplots all use ggpubr::stat_compare_means() for their p-value
# text -- remove it and add a manual, threshold-formatted label at the same
# top-center position (x = 1.5, the midpoint of the 2 factor levels).
# geom_text/annotate's `size` is in mm, not points (ggplot2's `.pt` constant
# converts pt -> mm: 1 / .pt). Original p-value text was size=2 (~6pt).
pvalue_text_size <- 2 + 1 / .pt

add_pvalue_bucket <- function(p, y_var, group_var) {
  d <- p$data
  p_val <- wilcox.test(d[[y_var]] ~ d[[group_var]])$p.value
  p <- drop_layers_by_class(p, "stat", "StatCompareMeans")
  p + annotate("text", x = 1.5, y = Inf, label = format_pvalue_threshold(p_val),
               vjust = 1.5, size = pvalue_text_size, family = "Arial", parse = TRUE)
}

p_b <- add_pvalue_bucket(p_b, "hub_rank", "rank")
p_c_fib <- add_pvalue_bucket(p_c_fib, "n_hub_genes", "rank_category")
p_c_card <- add_pvalue_bucket(p_c_card, "n_hub_genes", "rank_category")
p_d <- add_pvalue_bucket(p_d, "int_count", "rank")

# Fig3A's p-value text is baked into its own data (the `label` column,
# consumed by a geom_text() layer) rather than a stat_compare_means layer --
# recompute that column with the threshold formatter, then rebuild the
# geom_text layer with parse = TRUE so "10^-N" renders as a real superscript.
pa_data <- p_a$data
first_quintile <- levels(pa_data$quintile)[1]
pa_data$label <- ifelse(
  pa_data$quintile == first_quintile & pa_data$p.value < 0.05,
  vapply(pa_data$p.value, format_pvalue_threshold, character(1)),
  NA_character_
)
p_a <- drop_layers_by_class(p_a, "geom", "GeomText")
p_a$data <- pa_data
p_a <- p_a + geom_text(aes(label = label), position = position_nudge(y = 1),
                        size = pvalue_text_size, family = "Arial", parse = TRUE)

# Composite-only font bump: axis tick labels and y-axis titles were too small
# at the panels' original 6pt (set for the standalone panel PDFs, which are
# NOT touched here). axis.title.x = element_text() only reaches panel A,
# since B/C/D have no x-axis title.
axis_text_bump <- theme(axis.text = element_text(size = 7))
y_title_bump <- theme(axis.title.y = element_text(size = 8))
x_title_bump <- theme(axis.title.x = element_text(size = 8))
subplot_title_bump <- theme(plot.title = element_text(size = 8))

# Panel A's own title ("Enrichment of eGenes across dispersion quintiles")
# is redundant with the manuscript's figure caption -- dropped. labs(title =
# NULL) removes it entirely (not just blanked -- no reserved title space
# left behind).
p_a <- p_a + axis_text_bump + y_title_bump + x_title_bump + labs(title = NULL)
p_b <- p_b + axis_text_bump + y_title_bump
p_c_fib <- p_c_fib + axis_text_bump + y_title_bump + subplot_title_bump
p_c_card <- p_c_card + axis_text_bump + y_title_bump + subplot_title_bump
p_d <- p_d + axis_text_bump + y_title_bump

# C is itself a 2-panel row (Cardiac_fibroblasts | Cardiomyocytes), matching
# the published figure's sub-layout -- combined into one patchwork first,
# then wrap_elements() so the outer grid's auto-tagging treats it as a
# single panel ("C") instead of descending into its two leaf plots.
p_c_combined <- p_c_fib + p_c_card + plot_layout(ncol = 2)
p_c_unit <- wrap_elements(full = p_c_combined)

# Bottom row: C holds 2 sub-panels (Cardiac_fibroblasts + Cardiomyocytes,
# each with its own axis labels) vs. D's 1, so C needs roughly double D's
# width. widths = c(2, 1) gives each of C's two sub-panels about the same
# width as D, matching A/B's top-row proportions.
bottom_row <- p_c_unit + p_d + plot_layout(widths = c(2, 1))

# A/B/D's tag letters inherit each panel's own small `text` theme (6pt,
# cascading down to plot.tag by ggplot2's theme hierarchy), while C's tag
# (on the wrap_elements() wrapper, not a plain ggplot panel) falls back to
# patchwork's own larger default plot.tag size instead -- setting plot.tag
# explicitly on the outer plot_annotation() theme overrides every panel's
# inherited size uniformly, so all four tags render at the same (larger) size.
fig3_composite <- (p_a + p_b) / bottom_row +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

out_pdf <- file.path(fig3_dir, "Figure3_composite_excl_topexpressed.pdf")
ggsave(out_pdf, fig3_composite, device = cairo_pdf, width = 7.2, height = 6.4, units = "in")
message("[done] wrote ", out_pdf)
