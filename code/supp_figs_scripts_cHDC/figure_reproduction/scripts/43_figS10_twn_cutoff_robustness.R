#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(patchwork)
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
cluster_dir <- file.path(repo_dir, "cellranger_cluster-mode_trial")
analysis_dir <- file.path(repo_dir, "figure_reproduction")
out_dir <- file.path(analysis_dir, "outputs")
filt_mat_dir <- file.path(out_dir, "filtered_matrices")
figS10_dir <- file.path(out_dir, "figS10")
dir.create(figS10_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(analysis_dir, "scripts", "utils.R"))  # loads ggpubr

twn_s2_xlsx <- file.path(cluster_dir, "analysis", "datasets", "external_data",
                          "saha_coexp_net_TWN_2017", "Supplemental_Data_S2.xlsx")

# Cell type -> TWN S2 sheet (same sheet indices as script 27 / source Rmd).
targets <- list(Cardiac_fibroblasts = list(sheet_index = 30L, sheet_prefix = "Fibroblasts"),
                Cardiomyocytes      = list(sheet_index = 34L, sheet_prefix = "Heart"))
percentages <- c(0.1, 0.2)  # manuscript S10: 10% and 20% tails (main Fig3C: 30%)

# ------------------------------- HELPERS -------------------------------------
# Same low/high definition and TWN join as script 27's
# ct_specific_twn_percentage_retained(); returns the plotting data frame.
twn_low_high_df <- function(cell_type, sheet, percentage) {
  disp_path <- file.path(filt_mat_dir, paste0(cell_type, "_dispersion_retained.csv"))
  if (!file.exists(disp_path)) {
    stop("Missing ", disp_path, " -- run 01_exclude_topexpressed_genes.R first")
  }
  disp_mat <- as.matrix(read.csv(disp_path, row.names = 1, check.names = FALSE))
  med_disp <- apply(disp_mat, 1, median)
  ct_ranks <- rank(med_disp, ties.method = "average")

  n_genes <- length(ct_ranks)
  low_disp <- names(ct_ranks)[ct_ranks <= percentage * n_genes]
  high_disp <- names(ct_ranks)[ct_ranks >= n_genes - percentage * n_genes]

  te_te <- read_excel(twn_s2_xlsx, sheet = sheet, col_names = FALSE)
  colnames(te_te) <- c("hub_name", "n_genes", "genes")

  df <- data.frame(genes = c(low_disp, high_disp),
                   rank_category = c(rep("low", length(low_disp)), rep("high", length(high_disp)))) %>%
    dplyr::filter(genes %in% te_te$hub_name) %>%
    dplyr::mutate(n_hub_genes = as.integer(te_te$n_genes[match(genes, te_te$hub_name)]))
  df$rank_category <- factor(df$rank_category, levels = c("low", "high"), ordered = TRUE)
  message(sprintf("[figS10:%s:%.0f%%] %d retained genes ranked; %d/%d low, %d/%d high found as TWN hubs",
                  cell_type, 100 * percentage, n_genes,
                  sum(df$rank_category == "low"), length(low_disp),
                  sum(df$rank_category == "high"), length(high_disp)))
  if (any(is.na(df$n_hub_genes))) stop("NA hub connection count for ", cell_type)
  if (sum(df$rank_category == "low") == 0 || sum(df$rank_category == "high") == 0) {
    stop("Empty low or high group for ", cell_type, " at ", percentage)
  }
  df
}

# Plot styling copied from script 27's Fig3C panel; y_max shared across panels.
twn_panel <- function(df, title, y_max) {
  ggplot(df, aes(x = rank_category, y = n_hub_genes, fill = rank_category)) +
    geom_boxplot(width = 0.6, linewidth = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 0.6, alpha = 0.4) +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 2) +
    scale_fill_manual(values = c("low" = "steelblue", "high" = "coral")) +
    scale_x_discrete(labels = c("low" = "Low-dispersion\ngenes", "high" = "High-dispersion\ngenes")) +
    coord_cartesian(ylim = c(0, y_max)) +
    labs(x = NULL, y = "Connections per gene", title = title) +
    theme_classic(base_family = "Arial") +
    theme(
      text = element_text(size = 6),
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 6),
      plot.title = element_text(size = 6, hjust = 0.5),
      axis.line = element_line(linewidth = 0.4),
      axis.text.x = element_text(hjust = 0.5),
      legend.position = "none",
      plot.margin = margin(2, 2, 2, 2, unit = "mm"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4)
    )
}

# --------------------------------- MAIN --------------------------------------
if (!file.exists(twn_s2_xlsx)) stop("Missing TWN table: ", twn_s2_xlsx)
s2_sheets <- excel_sheets(twn_s2_xlsx)

panel_dfs <- list()
stats <- list()
for (cell_type in names(targets)) {
  tg <- targets[[cell_type]]
  sheet <- s2_sheets[tg$sheet_index]
  if (!startsWith(sheet, tg$sheet_prefix)) {
    stop(sprintf("S2 sheet %d is '%s', expected prefix '%s'", tg$sheet_index, sheet, tg$sheet_prefix))
  }
  for (pct in percentages) {
    df <- twn_low_high_df(cell_type, sheet, pct)
    wt <- wilcox.test(n_hub_genes ~ rank_category, data = df)
    stats[[length(stats) + 1]] <- data.frame(
      cell_type = cell_type, percentage = pct,
      n_low_hubs = sum(df$rank_category == "low"), n_high_hubs = sum(df$rank_category == "high"),
      median_low = median(df$n_hub_genes[df$rank_category == "low"]),
      median_high = median(df$n_hub_genes[df$rank_category == "high"]),
      wilcox_p = wt$p.value)
    panel_dfs[[sprintf("%s, %.0f%% tails", gsub("_", " ", cell_type), 100 * pct)]] <- df
  }
}
stats_df <- bind_rows(stats)
y_max <- max(vapply(panel_dfs, function(d) max(d$n_hub_genes), numeric(1)))
panels <- lapply(names(panel_dfs), function(t) twn_panel(panel_dfs[[t]], t, y_max))
print(stats_df)
write.csv(stats_df, file.path(figS10_dir, "figS10_wilcox.csv"), row.names = FALSE)

# Rows = cell types (A, B), columns = 10% | 20%; tag only the first panel of each row.
composite <- wrap_plots(panels, ncol = 2) +
  plot_annotation(tag_levels = list(c("A", "", "B", ""))) &
  theme(plot.tag = element_text(size = 8, face = "bold", family = "Arial"))

out_pdf <- file.path(figS10_dir, "FigureS10_composite.pdf")
ggsave(out_pdf, composite, device = cairo_pdf, width = 4.4, height = 5.08, units = "in")
ggsave(sub("\\.pdf$", ".png", out_pdf), composite, device = ragg::agg_png,
       width = 4.4, height = 5.08, units = "in", dpi = 300)

write.csv(
  data.frame(
    key = sprintf("figS10_%s_%.0fpct_wilcox_p", stats_df$cell_type, 100 * stats_df$percentage),
    value = signif(stats_df$wilcox_p, 4),
    provenance = "figure_reproduction/scripts/43_figS10_twn_cutoff_robustness.R"),
  file.path(out_dir, "reportable_values_43.csv"), row.names = FALSE)

message("[done] wrote ", out_pdf, " and .png")
