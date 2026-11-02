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
script_dir <- file.path(repo_dir, "figure_reproduction", "scripts", "q_threshold_sensitivity")
refit_dir <- file.path(analysis_dir, "dream_refit")
quintile_dir <- file.path(analysis_dir, "quintiles")
fig_dir <- file.path(analysis_dir, "figures")
tab_dir <- file.path(analysis_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

q_values <- c("0.05", "0.1", "0.2", "0.25")
baseline_q <- "0.25"  # most permissive swept value; also the manuscript's own canonical q
other_qs <- setdiff(q_values, baseline_q)

cell_types <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)
n_quantiles <- 5

source(file.path(script_dir, "utils.R"))

# --------------------------------- LOAD --------------------------------------

gene_universe <- setNames(lapply(q_values, function(q) {
  f <- file.path(refit_dir, paste0("q", q), "gene_universe.txt")
  if (!file.exists(f)) stop("Missing gene universe for q=", q, ": ", f)
  readLines(f)
}), q_values)

dd_sets <- setNames(lapply(q_values, function(q) {
  f <- file.path(refit_dir, paste0("q", q), "dd_gene_sets.rds")
  if (!file.exists(f)) stop("Missing DD gene sets for q=", q, ": ", f)
  readRDS(f)
}), q_values)

contrasts <- names(dd_sets[[q_values[1]]])
for (q in q_values) {
  stopifnot(setequal(names(dd_sets[[q]]), contrasts))
}
message(sprintf("[load] %d contrasts x %d q values", length(contrasts), length(q_values)))

quintile_dfs <- setNames(lapply(q_values, function(q) {
  f <- file.path(quintile_dir, paste0("q", q), "quintile_df.csv")
  if (!file.exists(f)) stop("Missing quintile table for q=", q, ": ", f)
  read.csv(f, row.names = 1, check.names = FALSE)
}), q_values)

# =============================================================================
# 1. Gene universe size per q
# =============================================================================
universe_df <- tibble(
  q = q_values,
  genes = sapply(q_values, function(q) length(gene_universe[[q]]))
)
write.csv(universe_df, file.path(tab_dir, "universe_size_per_q.csv"), row.names = FALSE)

p1 <- ggplot(universe_df, aes(x = q, y = genes)) +
  geom_col(fill = "#4575b4", width = 0.6) +
  labs(title = "Gene universe size per q (post ^RP filter)",
       x = "q (capture_rate)", y = "Number of genes") +
  theme_bw()
ggsave(file.path(fig_dir, "01_universe_size_per_q.pdf"), p1, width = 6, height = 5)
print(universe_df)

# =============================================================================
# 2. DD gene counts vs q, per contrast
# =============================================================================
dd_counts <- expand_grid(q = q_values, contrast = contrasts) |>
  rowwise() |>
  mutate(n_dd = length(dd_sets[[q]][[contrast]])) |>
  ungroup()
write.csv(dd_counts, file.path(tab_dir, "dd_gene_count_vs_q.csv"), row.names = FALSE)

p2 <- ggplot(dd_counts, aes(x = as.numeric(q), y = n_dd, group = contrast, colour = contrast)) +
  geom_line() +
  geom_point(size = 1.5) +
  labs(title = "DD gene count vs. q, per contrast", x = "q", y = "Number of DD genes") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none") +
  facet_wrap(~contrast, scales = "free_y", ncol = 5)
ggsave(file.path(fig_dir, "02_dd_gene_count_vs_q.pdf"), p2, width = 14, height = 8)

# =============================================================================
# 3. Pairwise-q Jaccard overlap of DD calls, per contrast
# =============================================================================
q_pairs <- combn(q_values, 2, simplify = FALSE)

pairwise_overlap <- map_dfr(q_pairs, function(pair) {
  map_dfr(contrasts, function(ct) {
    a <- dd_sets[[pair[1]]][[ct]]
    b <- dd_sets[[pair[2]]][[ct]]
    n_both <- length(intersect(a, b))
    n_union <- length(union(a, b))
    tibble(
      q1 = pair[1], q2 = pair[2], contrast = ct,
      n_q1 = length(a), n_q2 = length(b), n_both = n_both,
      jaccard = if (n_union > 0) n_both / n_union else NA_real_)
  })
})
write.csv(pairwise_overlap, file.path(tab_dir, "dd_pairwise_jaccard_by_contrast.csv"), row.names = FALSE)

cat(sprintf("Median pairwise Jaccard (DD calls, all q-pairs x contrasts): %.3f\n",
            median(pairwise_overlap$jaccard, na.rm = TRUE)))
cat(sprintf("%% of (q-pair x contrast) combos with Jaccard > 0.5: %.1f%%\n",
            100 * mean(pairwise_overlap$jaccard > 0.5, na.rm = TRUE)))

avg_jaccard <- pairwise_overlap |>
  group_by(q1, q2) |>
  summarise(mean_jaccard = mean(jaccard, na.rm = TRUE), .groups = "drop")
jaccard_mat <- bind_rows(
  avg_jaccard,
  avg_jaccard |> rename(q1 = q2, q2 = q1),
  tibble(q1 = q_values, q2 = q_values, mean_jaccard = 1))
write.csv(jaccard_mat, file.path(tab_dir, "dd_jaccard_avg_matrix.csv"), row.names = FALSE)

p3 <- ggplot(jaccard_mat, aes(x = q1, y = q2, fill = mean_jaccard)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.0f%%", 100 * mean_jaccard)), size = 4) +
  scale_fill_gradient(low = "white", high = "#4575b4", limits = c(0, 1)) +
  labs(title = "DD calls -- mean %% identical across contrasts", x = "q", y = "q", fill = "Jaccard") +
  theme_bw()
ggsave(file.path(fig_dir, "03_dd_jaccard_heatmap.pdf"), p3, width = 6, height = 5)

# =============================================================================
# 6a. Core/consistent DD-call fraction per contrast
# =============================================================================
dd_core <- map_dfr(contrasts, function(ct) {
  all_sets <- lapply(q_values, function(q) dd_sets[[q]][[ct]])
  union_genes <- Reduce(union, all_sets)
  core_genes <- Reduce(intersect, all_sets)
  tibble(contrast = ct, n_union = length(union_genes), n_core = length(core_genes),
         pct_core = if (length(union_genes) > 0) 100 * length(core_genes) / length(union_genes) else NA_real_)
}) |> arrange(pct_core)
write.csv(dd_core, file.path(tab_dir, "dd_core_fraction_by_contrast.csv"), row.names = FALSE)

p6a <- ggplot(dd_core, aes(x = reorder(contrast, pct_core), y = pct_core)) +
  geom_col(fill = "#4575b4", width = 0.75) +
  labs(title = "%% of ever-DD genes called DD at EVERY swept q, per contrast",
       x = NULL, y = "%% core (DD at every q)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
ggsave(file.path(fig_dir, "06a_dd_core_fraction_by_contrast.pdf"), p6a, width = 10, height = 5)

cat(sprintf("\nMedian DD-call core fraction across contrasts: %.1f%%\n", median(dd_core$pct_core, na.rm = TRUE)))

# =============================================================================
# 6b. Core/consistent dispersion-quintile fraction per cell type
#     (consistently-low = quintile 1, consistently-high = quintile 5)
# =============================================================================
quintile_sets <- function(target_quintile) {
  setNames(lapply(cell_types, function(ct) {
    setNames(lapply(q_values, function(q) {
      qdf <- quintile_dfs[[q]]
      rownames(qdf)[!is.na(qdf[[ct]]) & qdf[[ct]] == target_quintile]
    }), q_values)
  }), cell_types)
}
low_sets <- quintile_sets(1)
high_sets <- quintile_sets(n_quantiles)

quintile_core <- function(sets_by_ct, label) {
  map_dfr(cell_types, function(ct) {
    all_sets <- sets_by_ct[[ct]]
    union_genes <- Reduce(union, all_sets)
    core_genes <- Reduce(intersect, all_sets)
    tibble(celltype = ct, set = label, n_union = length(union_genes), n_core = length(core_genes),
           pct_core = if (length(union_genes) > 0) 100 * length(core_genes) / length(union_genes) else NA_real_)
  })
}
quintile_core_df <- bind_rows(
  quintile_core(low_sets, "Consistently low (Q1)"),
  quintile_core(high_sets, "Consistently high (Q5)"))
write.csv(quintile_core_df, file.path(tab_dir, "quintile_core_fraction_by_celltype.csv"), row.names = FALSE)

p6b <- ggplot(quintile_core_df, aes(x = reorder(celltype, pct_core), y = pct_core, fill = set)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Consistently low (Q1)" = "#4575b4", "Consistently high (Q5)" = "#d73027")) +
  labs(title = "%% of ever-extreme-quintile genes classified consistently across q, per cell type",
       x = NULL, y = "%% core (same quintile at every q)", fill = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "top")
ggsave(file.path(fig_dir, "06b_quintile_core_fraction_by_celltype.pdf"), p6b, width = 8, height = 5)

cat(sprintf("Median quintile core fraction (Q1): %.1f%%; (Q5): %.1f%%\n",
            median(quintile_core_df$pct_core[quintile_core_df$set == "Consistently low (Q1)"], na.rm = TRUE),
            median(quintile_core_df$pct_core[quintile_core_df$set == "Consistently high (Q5)"], na.rm = TRUE)))

# =============================================================================
# 7. Retention/fate of baseline-q (0.25) DD genes as q decreases, per contrast
# =============================================================================
retention_breakdown <- map_dfr(contrasts, function(ct) {
  baseline_genes <- dd_sets[[baseline_q]][[ct]]
  n_baseline <- length(baseline_genes)

  bind_rows(
    tibble(contrast = ct, q = baseline_q, n_baseline = n_baseline,
           n_retained = n_baseline, n_lost_not_sig = 0L, n_lost_dropped = 0L,
           retention_pct = 100, lost_not_sig_pct = 0, lost_dropped_pct = 0),
    map_dfr(other_qs, function(q) {
      universe_q <- gene_universe[[q]]
      retained <- intersect(baseline_genes, dd_sets[[q]][[ct]])
      dropped <- setdiff(baseline_genes, universe_q)
      lost_not_sig <- setdiff(setdiff(baseline_genes, retained), dropped)
      tibble(contrast = ct, q = q, n_baseline = n_baseline,
             n_retained = length(retained), n_lost_not_sig = length(lost_not_sig),
             n_lost_dropped = length(dropped),
             retention_pct = if (n_baseline > 0) 100 * length(retained) / n_baseline else NA_real_,
             lost_not_sig_pct = if (n_baseline > 0) 100 * length(lost_not_sig) / n_baseline else NA_real_,
             lost_dropped_pct = if (n_baseline > 0) 100 * length(dropped) / n_baseline else NA_real_)
    }))
})
write.csv(retention_breakdown, file.path(tab_dir, "retention_fate_by_contrast.csv"), row.names = FALSE)

status_colours <- c("Retained" = "#4575b4", "No Longer Significant" = "#fdae61", "Dropped" = "grey70")

ct_order <- retention_breakdown |> distinct(contrast, n_baseline) |> arrange(desc(n_baseline)) |> pull(contrast)
label_lookup <- retention_breakdown |> distinct(contrast, n_baseline)

# Split "compare_{CT1}_{CT2}" into a two-row "CT1\nvs CT2" facet label (the
# single-row "compare_CT1_CT2 (n=x)" label was too long and got cut off).
split_contrast_label <- function(ct, cell_types) {
  ct_no_prefix <- sub("^compare_", "", ct)
  ct_sorted <- cell_types[order(-nchar(cell_types))]
  for (ctype in ct_sorted) {
    prefix <- paste0(ctype, "_")
    if (startsWith(ct_no_prefix, prefix)) {
      second <- sub(paste0("^", prefix), "", ct_no_prefix)
      if (second %in% cell_types) return(sprintf("%s\nvs %s", ctype, second))
    }
  }
  stop("Could not parse contrast into two known cell types: ", ct)
}
label_map <- setNames(vapply(label_lookup$contrast, split_contrast_label, character(1), cell_types = cell_types),
                       label_lookup$contrast)

retention_long <- retention_breakdown |>
  select(contrast, q, retention_pct, lost_not_sig_pct, lost_dropped_pct) |>
  mutate(contrast = factor(unname(label_map[contrast]), levels = unname(label_map[ct_order]))) |>
  pivot_longer(cols = c(retention_pct, lost_not_sig_pct, lost_dropped_pct),
               names_to = "status", values_to = "pct") |>
  mutate(status = recode(status,
           retention_pct = "Retained", lost_not_sig_pct = "No Longer Significant", lost_dropped_pct = "Dropped"),
         status = factor(status, levels = c("Dropped", "No Longer Significant", "Retained")),
         q = factor(q, levels = q_values))

p7 <- ggplot(retention_long, aes(x = q, y = pct, fill = status)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = status_colours) +
  labs(x = "q", y = sprintf("%% of DD genes (q=%s baseline)", baseline_q), fill = NULL) +
  theme_bw(base_size = 7) +
  theme(legend.position = "top", legend.text = element_text(size = 6), legend.key.size = unit(0.3, "cm"),
        axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),
        strip.text = element_text(size = 5.5), plot.margin = margin(2, 2, 2, 2)) +
  facet_wrap(~contrast, ncol = 5)
# Sized to fit within half a standard (US Letter, 1in margins) page: 6.5in text
# width x ~4.5in (half of the ~9in text height).
ggsave(file.path(fig_dir, "07_retention_fate_by_contrast.pdf"), p7, width = 6.5, height = 4.5)

overall_retention <- retention_breakdown |>
  group_by(q) |>
  summarise(mean_retention = mean(retention_pct, na.rm = TRUE),
            mean_lost_not_sig = mean(lost_not_sig_pct, na.rm = TRUE),
            mean_lost_dropped = mean(lost_dropped_pct, na.rm = TRUE), .groups = "drop")
write.csv(overall_retention, file.path(tab_dir, "retention_fate_overall.csv"), row.names = FALSE)
print(overall_retention)

overall_long <- overall_retention |>
  pivot_longer(cols = c(mean_retention, mean_lost_not_sig, mean_lost_dropped),
               names_to = "status", values_to = "pct") |>
  mutate(status = recode(status,
           mean_retention = "Retained", mean_lost_not_sig = "No Longer Significant", mean_lost_dropped = "Dropped"),
         status = factor(status, levels = c("Dropped", "No Longer Significant", "Retained")),
         q = factor(q, levels = q_values))

p7b <- ggplot(overall_long, aes(x = q, y = pct, fill = status)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = status_colours) +
  labs(title = sprintf("Mean fate of q=%s DD genes vs. q (averaged across %d contrasts)", baseline_q, length(contrasts)),
       x = "q", y = "%% of baseline DD genes", fill = NULL) +
  theme_bw() +
  theme(legend.position = "top")
ggsave(file.path(fig_dir, "07b_retention_fate_overall.pdf"), p7b, width = 7, height = 5)

message("\n[done] figures in ", fig_dir, "; tables in ", tab_dir)
