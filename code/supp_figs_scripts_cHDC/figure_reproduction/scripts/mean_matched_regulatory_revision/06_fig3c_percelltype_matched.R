#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readxl)
})

args       <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
source(file.path(script_dir, "config.R"))
panel_dir <- ensure_dir(file.path(out_dir, "fig3_fig4_panels"))

TWN_S2_XLSX <- file.path(src_dir, "analysis/datasets/external_data/saha_coexp_net_TWN_2017/Supplemental_Data_S2.xlsx")

# Fig3C's own design (Rmd 1639/1642): two cell types, each against its own TWN
# sheet, same 30% low/high split for both. Named here, not buried inline.
FIG3C_TARGETS <- list(
  Cardiac_fibroblasts = list(sheet_index = 30L, sheet_name_prefix = "Fibroblasts"),
  Cardiomyocytes       = list(sheet_index = 34L, sheet_name_prefix = "Heart")
)
FIG3C_PERCENTAGE <- 0.3

# ---- Inputs -------------------------------------------------------------
gene_sets <- readRDS(file.path(out_dir, "gene_sets", "dispersion_gene_sets.rds"))
assert_that(!is.null(gene_sets$ct_ranks_df),
            "dispersion_gene_sets.rds has no ct_ranks_df -- re-run 01_build_dispersion_gene_sets.R (updated 2026-09-03) first")
ct_ranks_df <- gene_sets$ct_ranks_df

assert_that(file.exists(mean_matrix), paste("missing:", mean_matrix))
mean_df <- read.csv(mean_matrix, row.names = 1)
# ANALYSIS_OK[filtering]: drop the same ribo (^RP) genes 01/02 remove, keeping
# this script's own mean matrix on the same gene universe as ct_ranks_df;
# retained count logged immediately below.
mean_df <- mean_df[!rownames(mean_df) %in% gene_sets$ribo_genes_removed, , drop = FALSE]
log_step("loaded mean matrix after ribo filter (genes)", nrow(mean_df))
assert_that(all(rownames(ct_ranks_df) %in% rownames(mean_df)),
            "some ct_ranks_df genes are absent from the mean matrix")

assert_that(file.exists(TWN_S2_XLSX), paste("missing:", TWN_S2_XLSX))
s2_sheets <- excel_sheets(TWN_S2_XLSX)

# ---- Per-cell-type 2-group mean-matching (generalizes 02's match_by_mean) --
# For ONE cell type: rank genes by that cell type's own mean expression
# (median across its 3 individual columns), bin into n_mean_percentile_bins
# percentiles, and downsample low_disp/high_disp to the smaller count per bin.
# ANALYSIS_OK[within-file-helper]: called below per cell type per seed; local
# to this script, not shared, so no cross-file caller by design.
match_by_mean_one_celltype <- function(cell_type, low_disp, high_disp, seed) {
  cols <- paste0(individuals, "_", cell_type)
  assert_that(all(cols %in% colnames(mean_df)),
              sprintf("missing mean-expression columns for %s", cell_type))
  ct_mean <- apply(mean_df[, cols, drop = FALSE], 1, median)
  ct_mean_rank <- rank(ct_mean, ties.method = "average")

  bin_of <- cut(ct_mean_rank,
                breaks = quantile(ct_mean_rank, probs = seq(0, 1, length.out = n_mean_percentile_bins + 1),
                                  na.rm = TRUE),
                labels = FALSE, include.lowest = TRUE)
  names(bin_of) <- rownames(mean_df)

  # ANALYSIS_OK[seeded-sweep]: seed is an explicit argument, set once at entry,
  # matching 02's identical pattern -- called sequentially over match_seeds,
  # not in parallel, for the stability check.
  set.seed(seed)
  low_equal <- character(0); high_equal <- character(0)
  for (b in sort(unique(c(bin_of[low_disp], bin_of[high_disp])))) {
    low_genes  <- low_disp[bin_of[low_disp] == b]
    high_genes <- high_disp[bin_of[high_disp] == b]
    n_sample <- min(length(low_genes), length(high_genes))
    # ANALYSIS_OK[mean-matching]: sample() downsamples real genes to equalize
    # per-bin counts between the low/high groups (mean-matching), drawing from
    # the observed gene set; seeded above.
    low_keep  <- if (length(low_genes)  > n_sample) sample(low_genes,  n_sample) else low_genes
    high_keep <- if (length(high_genes) > n_sample) sample(high_genes, n_sample) else high_genes
    low_equal  <- c(low_equal, low_keep)
    high_equal <- c(high_equal, high_keep)
  }
  list(low = low_equal, high = high_equal)
}

# ANALYSIS_OK[within-file-helper]: called below per cell type; local to this
# script, not shared, so no cross-file caller by design.
assert_balanced_two_group <- function(matched, cell_type, seed) {
  assert_that(length(matched$low) == length(matched$high),
              sprintf("%s seed %d: matched low/high set sizes differ (%d vs %d) -- per-bin balance broken",
                      cell_type, seed, length(matched$low), length(matched$high)))
}

# ---- Per-cell-type: define low/high, mean-match, join TWN hub table, test --
for (cell_type in names(FIG3C_TARGETS)) {
  log_step(sprintf("=== Fig3C: %s ===", cell_type))
  target <- FIG3C_TARGETS[[cell_type]]
  assert_that(grepl(paste0("^", target$sheet_name_prefix), s2_sheets[target$sheet_index]),
              sprintf("Supplemental_Data_S2.xlsx sheet %d is '%s', expected prefix '%s' -- file layout changed",
                      target$sheet_index, s2_sheets[target$sheet_index], target$sheet_name_prefix))

  # Rank only this cell type's retained genes (NA = excluded by -6 here).
  ct_ranks <- ct_ranks_df[[cell_type]]
  names(ct_ranks) <- rownames(ct_ranks_df)
  ct_ranks <- ct_ranks[!is.na(ct_ranks)]
  n_genes_ct <- length(ct_ranks)
  assert_that(min(ct_ranks) >= 1 && max(ct_ranks) <= n_genes_ct,
              sprintf("%s ranks fall outside 1..%d (its retained genes)", cell_type, n_genes_ct))
  low_cut  <- FIG3C_PERCENTAGE * n_genes_ct
  high_cut <- n_genes_ct - FIG3C_PERCENTAGE * n_genes_ct
  low_disp  <- names(ct_ranks)[ct_ranks <= low_cut]
  high_disp <- names(ct_ranks)[ct_ranks >= high_cut]
  log_step(sprintf("  %s: %d genes total, low_disp=%d, high_disp=%d",
                   cell_type, n_genes_ct, length(low_disp), length(high_disp)))

  te_te <- read_excel(TWN_S2_XLSX, sheet = s2_sheets[target$sheet_index], col_names = FALSE)
  colnames(te_te) <- c("hub_name", "n_genes", "genes")
  te_te$n_genes <- as.integer(te_te$n_genes)

  matched_by_seed <- setNames(lapply(match_seeds, function(s) {
    m <- match_by_mean_one_celltype(cell_type, low_disp, high_disp, s)
    assert_balanced_two_group(m, cell_type, s)
    log_step(sprintf("  seed %d matched-set sizes: low=%d high=%d", s, length(m$low), length(m$high)))
    m
  }), as.character(match_seeds))

  hub_lookup <- function(genes) as.integer(te_te$n_genes[match(genes, te_te$hub_name)])

  primary <- matched_by_seed[[as.character(match_seed_primary)]]
  primary_labelled <- bind_rows(
    data.frame(gene = primary$low,  rank_category = "Low-dispersion"),
    data.frame(gene = primary$high, rank_category = "High-dispersion"))
  hub_frame <- primary_labelled %>%
    filter(gene %in% te_te$hub_name) %>%
    mutate(n_genes = hub_lookup(gene))
  log_step(sprintf("  matched genes that are TWN hubs in %s: %d", cell_type, nrow(hub_frame)))

  res <- wilcox_effect(hub_frame$n_genes, hub_frame$rank_category,
                        "Low-dispersion", "High-dispersion", "twn_connections_per_gene")
  write.csv(res, file.path(panel_dir, sprintf("fig3c_%s_twn_wilcox_matched.csv", cell_type)), row.names = FALSE)
  print(res)

  stability <- bind_rows(lapply(names(matched_by_seed), function(s) {
    m <- matched_by_seed[[s]]
    lab <- bind_rows(data.frame(gene = m$low,  rank_category = "Low-dispersion"),
                     data.frame(gene = m$high, rank_category = "High-dispersion"))
    v <- hub_lookup(lab$gene)
    keep <- !is.na(v)
    r <- wilcox_effect(v[keep], lab$rank_category[keep],
                       "Low-dispersion", "High-dispersion", "twn_connections_per_gene")
    r$seed <- s; r
  }))
  write.csv(stability, file.path(panel_dir, sprintf("fig3c_%s_seed_stability.csv", cell_type)), row.names = FALSE)
  assert_that(length(unique(sign(stability$hl_shift))) == 1,
              sprintf("Fig3C (%s) Low-vs-High effect changes SIGN across matching seeds -- not stable", cell_type))
  log_step(sprintf("  Fig3C (%s) Low-vs-High effect direction is stable across all matching seeds", cell_type))

  hub_frame$rank_category <- factor(hub_frame$rank_category, levels = c("Low-dispersion", "High-dispersion"))
  p <- ggplot(hub_frame, aes(x = rank_category, y = n_genes, fill = rank_category)) +
    geom_boxplot(width = 0.6, linewidth = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 0.6, alpha = 0.4) +
    # p-value annotation, same style/wording as 4B's forest-plot labels
    # (format_pvalue_label(), config.R) -- Brendan's request, 2026-09-03.
    annotate("text", x = 1.5, y = Inf, label = format_pvalue_label(res$p_value),
             vjust = 1.3, size = 2, family = "Arial") +
    scale_fill_manual(values = c("Low-dispersion" = "steelblue", "High-dispersion" = "coral")) +
    scale_x_discrete(labels = c("Low-dispersion" = "Low-dispersion\ngenes",
                                 "High-dispersion" = "High-dispersion\ngenes")) +
    labs(x = NULL, y = "Connections per gene",
         title = sprintf("Mean-matched TWN connectivity (%s)", gsub("_", " ", cell_type))) +
    theme_classic(base_family = "Arial") +
    theme(text = element_text(size = 6), axis.title = element_text(size = 6),
          axis.text = element_text(size = 6), axis.line = element_line(linewidth = 0.4),
          axis.text.x = element_text(hjust = 0.5), legend.position = "none",
          plot.margin = margin(2, 2, 2, 2, unit = "mm"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4))
  ggsave(file.path(panel_dir, sprintf("fig3c_%s_matched_boxplot.pdf", cell_type)), p,
         width = 3.2, height = 2.8, units = "in", device = cairo_pdf)
  saveRDS(p, file.path(panel_dir, sprintf("fig3c_%s_plot.rds", cell_type)))  # for 07's combined-figure patchwork
}

log_step(paste("Fig3C mean-matched results written to", panel_dir))
