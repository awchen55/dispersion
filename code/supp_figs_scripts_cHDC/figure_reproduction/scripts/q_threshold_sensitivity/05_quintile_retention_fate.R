#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(tibble)
})

# ------------------------------- CONSTANTS -----------------------------------
# PATCHED for figure_reproduction (see PORT_CHANGELOG_qthreshold.md).
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
analysis_dir <- file.path(repo_dir, "figure_reproduction", "outputs", "q_threshold_sensitivity")
quintile_dir <- file.path(analysis_dir, "quintiles")
fig_dir <- file.path(analysis_dir, "figures")
tab_dir <- file.path(analysis_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# ANALYSIS_OK[implicit-file-selection]: matches script 04's identical
# constants -- same swept values, same q=0.25 baseline (the manuscript's own
# canonical capture_rate, also the most permissive swept value).
q_values <- c("0.05", "0.1", "0.2", "0.25")
baseline_q <- "0.25"
other_qs <- setdiff(q_values, baseline_q)

cell_types <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)
n_quantiles <- 5

# --------------------------------- LOAD --------------------------------------

quintile_dfs <- setNames(lapply(q_values, function(q) {
  f <- file.path(quintile_dir, paste0("q", q), "quintile_df.csv")
  if (!file.exists(f)) stop("Missing quintile table for q=", q, ": ", f)
  read.csv(f, row.names = 1, check.names = FALSE)
}), q_values)

# =============================================================================
# Retention/fate of baseline-q (0.25) Q1/Q5 quintile genes as q decreases,
# per cell type, then averaged across cell types (pooled, no cell-type axis
# in the final plot).
# =============================================================================
set_targets <- c("Consistently low (Q1)" = 1, "Consistently high (Q5)" = n_quantiles)

retention_breakdown <- map_dfr(names(set_targets), function(set_label) {
  target_q <- set_targets[[set_label]]
  map_dfr(cell_types, function(ct) {
    baseline_df <- quintile_dfs[[baseline_q]]
    baseline_genes <- rownames(baseline_df)[!is.na(baseline_df[[ct]]) &
                                               baseline_df[[ct]] == target_q]
    n_baseline <- length(baseline_genes)

    bind_rows(
      tibble(set = set_label, celltype = ct, q = baseline_q, n_baseline = n_baseline,
             n_retained = n_baseline, n_reclassified = 0L, n_dropped = 0L,
             retention_pct = 100, reclassified_pct = 0, dropped_pct = 0),
      map_dfr(other_qs, function(q) {
        q_df <- quintile_dfs[[q]]
        universe_q <- rownames(q_df)
        dropped <- setdiff(baseline_genes, universe_q)
        present <- intersect(baseline_genes, universe_q)
        retained <- present[!is.na(q_df[present, ct]) & q_df[present, ct] == target_q]
        reclassified <- setdiff(present, retained)
        tibble(set = set_label, celltype = ct, q = q, n_baseline = n_baseline,
               n_retained = length(retained), n_reclassified = length(reclassified),
               n_dropped = length(dropped),
               retention_pct = if (n_baseline > 0) 100 * length(retained) / n_baseline else NA_real_,
               reclassified_pct = if (n_baseline > 0) 100 * length(reclassified) / n_baseline else NA_real_,
               dropped_pct = if (n_baseline > 0) 100 * length(dropped) / n_baseline else NA_real_)
      }))
  })
})
write.csv(retention_breakdown, file.path(tab_dir, "quintile_retention_fate_by_celltype.csv"), row.names = FALSE)

overall_retention <- retention_breakdown |>
  group_by(set, q) |>
  summarise(mean_retention = mean(retention_pct, na.rm = TRUE),
            mean_reclassified = mean(reclassified_pct, na.rm = TRUE),
            mean_dropped = mean(dropped_pct, na.rm = TRUE), .groups = "drop")
write.csv(overall_retention, file.path(tab_dir, "quintile_retention_fate_overall.csv"), row.names = FALSE)
print(overall_retention)

status_colours <- c("Retained" = "#4575b4", "Reclassified" = "#fdae61", "Dropped" = "grey70")

overall_long <- overall_retention |>
  pivot_longer(cols = c(mean_retention, mean_reclassified, mean_dropped),
               names_to = "status", values_to = "pct") |>
  mutate(status = recode(status,
           mean_retention = "Retained", mean_reclassified = "Reclassified", mean_dropped = "Dropped"),
         status = factor(status, levels = c("Dropped", "Reclassified", "Retained")),
         q = factor(q, levels = q_values))

p5 <- ggplot(overall_long, aes(x = q, y = pct, fill = status)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = status_colours) +
  labs(title = sprintf("Mean fate of q=%s dispersion-quintile genes vs. q (averaged across %d cell types)",
                        baseline_q, length(cell_types)),
       x = "q", y = "% of baseline quintile-set genes", fill = NULL) +
  facet_wrap(~ set) +
  theme_bw() +
  theme(legend.position = "top")
ggsave(file.path(fig_dir, "08_quintile_retention_fate_overall.pdf"), p5, width = 9, height = 5)

lowest_q <- other_qs[[which.min(as.numeric(other_qs))]]
cat(sprintf("\nRetention at q=%s (lowest swept q): Q1 %.1f%%, Q5 %.1f%%\n",
            lowest_q,
            overall_retention$mean_retention[overall_retention$set == "Consistently low (Q1)" &
                                                overall_retention$q == lowest_q],
            overall_retention$mean_retention[overall_retention$set == "Consistently high (Q5)" &
                                                overall_retention$q == lowest_q]))

message("\n[done] figure: ", file.path(fig_dir, "08_quintile_retention_fate_overall.pdf"))
message("[done] tables: ", file.path(tab_dir, "quintile_retention_fate_by_celltype.csv"), ", ",
        file.path(tab_dir, "quintile_retention_fate_overall.csv"))
