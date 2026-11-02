#!/usr/bin/env Rscript

# Column names ({IND}_{CT}) for one cell type across the three individuals.
ct_columns <- function(cell_type) paste0(individuals, "_", cell_type)

# Wilcoxon p-value for x vs. a baseline distribution, NA if x is too small to
# test meaningfully. Factored out of scripts 08/09 to remove the same R025
# drift hazard (was defined identically in both).
wilcox_p_vs_baseline <- function(x, baseline, min_n) {
  if (length(x) >= min_n) wilcox.test(x, baseline)$p.value else NA_real_
}

# Low-vs-high consistently-dispersed-gene boxplot, byte-for-byte copy of the
# styling shared by main_figs_1-4_related_supp.Rmd's Fig3B/C/D and Fig4A/C/D/E
# (geom_boxplot + geom_jitter + stat_compare_means wilcox + steelblue/coral
# fill + theme_classic Arial 6pt). `df` needs columns `rank` (factor levels
# c("1","5")) and `y_var`; `y_label` is the y-axis title. Factored out here
# (scripts 27/29) to avoid an 8-way R025 drift hazard across Fig3B/C(x2)/D and
# Fig4A/C/D/E, which are otherwise identical apart from the y aesthetic/label.
suppressPackageStartupMessages(library(ggpubr))
make_low_high_boxplot <- function(df, y_var, y_label) {
  ggplot(df, aes(x = rank, y = .data[[y_var]], fill = rank)) +
    geom_boxplot(width = 0.6, linewidth = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 0.6, alpha = 0.4) +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 2) +
    scale_fill_manual(values = c("1" = "steelblue", "5" = "coral")) +
    scale_x_discrete(labels = c("1" = "Low-dispersion\ngenes", "5" = "High-dispersion\ngenes")) +
    labs(x = NULL, y = y_label) +
    theme_classic(base_family = "Arial") +
    theme(
      text = element_text(size = 6),
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 6),
      axis.line = element_line(linewidth = 0.4),
      axis.text.x = element_text(hjust = 0.5),
      legend.position = "none",
      plot.margin = margin(2, 2, 2, 2, unit = "mm"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4)
    )
}
