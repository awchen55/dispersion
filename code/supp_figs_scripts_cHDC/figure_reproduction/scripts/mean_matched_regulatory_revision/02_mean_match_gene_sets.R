#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

args       <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
source(file.path(script_dir, "config.R"))
mm_dir <- ensure_dir(file.path(out_dir, "mean_matched"))

# ---- Inputs -----------------------------------------------------------------
gene_sets <- readRDS(file.path(out_dir, "gene_sets", "dispersion_gene_sets.rds"))
var_quintile_dfs        <- gene_sets$var_quintile_dfs
var_first_quint_genes   <- gene_sets$var_first_quint_genes
var_fifth_quint_genes   <- gene_sets$var_fifth_quint_genes
never_first_fifth_genes <- gene_sets$never_first_fifth_genes
ribo_genes_full         <- gene_sets$ribo_genes_removed

assert_that(file.exists(mean_matrix), paste("missing:", mean_matrix))
mean_df <- read.csv(mean_matrix, row.names = 1)
# ANALYSIS_OK[filtering]: drop the same ribo (^RP) genes removed from the dispersion
# matrix in 01 (ribo_genes_removed), keeping the two matrices on one gene set; the
# retained gene count is logged immediately below.
mean_df <- mean_df %>% filter(!rownames(.) %in% ribo_genes_full)   # Rmd 59-60
log_step("loaded mean matrix after ribo filter (genes)", nrow(mean_df))
# Matching is only defined for genes present in BOTH matrices.
assert_that(all(rownames(var_quintile_dfs) %in% rownames(mean_df)),
            "some quintile-table genes are absent from the mean matrix")

# ---- Global mean-expression rank + bins (Rmd 2903-2929) ---------------------
global_mean_ranks_df <- data.frame(
  gene        = rownames(mean_df),
  global_rank = rank(apply(mean_df, 1, median), ties.method = "average")
)

n_before_pd <- nrow(var_quintile_dfs)
plot_df <- data.frame(gene = rownames(var_quintile_dfs)) %>%
  left_join(global_mean_ranks_df, by = "gene")
stopifnot(nrow(plot_df) == n_before_pd)   # one-to-one gene join (unique gene keys both sides)
plot_df <- plot_df %>%
  mutate(first_quint_var  = gene %in% var_first_quint_genes,
         fifth_quint_var  = gene %in% var_fifth_quint_genes,
         never_first_fifth= gene %in% never_first_fifth_genes,
         other_genes      = !gene %in% var_first_quint_genes &
                            !gene %in% var_fifth_quint_genes &
                            !gene %in% never_first_fifth_genes)
assert_that(!anyNA(plot_df$global_rank), "gene missing a global mean rank")

# 5-bin (for the balance barplot) and 100-bin (for the actual matching) cuts.
plot_df$mean_quintiles <- cut(
  plot_df$global_rank,
  breaks = quantile(plot_df$global_rank, probs = seq(0, 1, length.out = n_quantiles + 1),
                    na.rm = TRUE),
  labels = FALSE, include.lowest = TRUE)
plot_df$mean_percentiles <- cut(
  plot_df$global_rank,
  breaks = quantile(plot_df$global_rank, probs = seq(0, 1, length.out = n_mean_percentile_bins + 1),
                    na.rm = TRUE),
  labels = FALSE, include.lowest = TRUE)

# ANALYSIS_OK[reshape-not-plot]: despite the legacy 'plot_' name (carried from the
# source Rmd), this is the long-form matching input, not a plotting-only view:
# each gene contributes one row per dispersion group it belongs to, and keeping
# is_member rows is the reshape that defines those memberships. It feeds
# match_by_mean() below; the downstream analysis population is exactly these rows.
plot_long <- plot_df %>%
  pivot_longer(cols = c(fifth_quint_var, first_quint_var, never_first_fifth, other_genes),
               names_to = "quintile_group", values_to = "is_member") %>%
  filter(is_member)

# ---- Matching function (Rmd 2941-2993) --------------------------------------
# For each mean-percentile bin, downsample first / fifth / never to the smallest
# of the three counts in that bin. Faithful to the source loop; wrapped so we
# can iterate seeds for the stability check.
# ANALYSIS_OK[within-file-helper]: called only in this script (seed loop below);
# a per-run matching helper, not shared, so no cross-file caller by design.
match_by_mean <- function(plot_long, seed) {
  # ANALYSIS_OK[seeded-sweep]: seed is an explicit argument, set once at entry, so
  # each matching run is reproducible; called sequentially over match_seeds (not
  # in parallel) precisely to sweep seeds for the stability check in 03.
  set.seed(seed)
  first_equal <- character(0)
  fifth_equal <- character(0)
  never_equal <- character(0)

  for (i in sort(unique(plot_long$mean_percentiles))) {
    first_genes <- plot_long$gene[plot_long$quintile_group == "first_quint_var"  & plot_long$mean_percentiles == i]
    fifth_genes <- plot_long$gene[plot_long$quintile_group == "fifth_quint_var"  & plot_long$mean_percentiles == i]
    never_genes <- plot_long$gene[plot_long$quintile_group == "never_first_fifth"& plot_long$mean_percentiles == i]

    n_sample <- min(length(first_genes), length(fifth_genes), length(never_genes))

    # sample() only when strictly larger than the target (avoids sample()'s
    # single-integer footgun and matches the source's if/else).
    # ANALYSIS_OK[mean-matching]: sample() downsamples real genes to equalize
    # per-bin counts across dispersion groups (mean-matching); it draws from the
    # observed gene set, not synthetic data, and is seeded above.
    first_keep <- if (length(first_genes) > n_sample) sample(first_genes, n_sample) else first_genes
    fifth_keep <- if (length(fifth_genes) > n_sample) sample(fifth_genes, n_sample) else fifth_genes
    never_keep <- if (length(never_genes) > n_sample) sample(never_genes, n_sample) else never_genes

    first_equal <- c(first_equal, first_keep)
    fifth_equal <- c(fifth_equal, fifth_keep)
    never_equal <- c(never_equal, never_keep)
  }
  list(first = first_equal, fifth = fifth_equal, never = never_equal)
}

# ---- Balance assertion ------------------------------------------------------
# Guarantee of the algorithm: within every occupied mean-percentile bin the
# three matched groups have identical counts. Fail loudly otherwise.
# ANALYSIS_OK[within-file-helper]: called only in this script (seed loop below);
# a per-run balance check, not shared, so no cross-file caller by design.
assert_balanced <- function(matched, plot_long, seed) {
  bin_of <- plot_long %>% dplyr::select(gene, mean_percentiles) %>% distinct()
  matched_long <- bind_rows(
    data.frame(gene = matched$first, group = "first"),
    data.frame(gene = matched$fifth, group = "fifth"),
    data.frame(gene = matched$never, group = "never")
  )
  n_before <- nrow(matched_long)
  matched_long <- matched_long %>% left_join(bin_of, by = "gene")
  stopifnot(nrow(matched_long) == n_before)  # one bin per gene (bin_of is distinct on gene)
  tab <- matched_long %>%
    count(mean_percentiles, group) %>%
    pivot_wider(names_from = group, values_from = n, values_fill = 0)
  unbalanced <- tab %>% filter(!(first == fifth & fifth == never))
  assert_that(nrow(unbalanced) == 0,
              sprintf("seed %d: %d mean-percentile bins are not balanced across groups",
                      seed, nrow(unbalanced)))
  log_step(sprintf("seed %d balance check passed across %d occupied bins",
                   seed, nrow(tab)))
}

# ---- Run matching for every seed --------------------------------------------
matched_by_seed <- list()
for (s in match_seeds) {
  m <- match_by_mean(plot_long, s)
  assert_balanced(m, plot_long, s)
  log_step(sprintf("seed %d matched-set sizes: first=%d fifth=%d never=%d",
                   s, length(m$first), length(m$fifth), length(m$never)))
  matched_by_seed[[as.character(s)]] <- m
}

# ---- Balance-check barplot for the primary seed (Rmd 2996-3024) -------------
primary <- matched_by_seed[[as.character(match_seed_primary)]]
gene_quintile_map <- plot_long %>% dplyr::select(gene, mean_quintiles) %>% distinct()
primary_long <- bind_rows(
  data.frame(gene = primary$first, group = "first"),
  data.frame(gene = primary$fifth, group = "fifth"),
  data.frame(gene = primary$never, group = "never")
)
n_before_dc <- nrow(primary_long)
primary_long <- primary_long %>% left_join(gene_quintile_map, by = "gene")
stopifnot(nrow(primary_long) == n_before_dc)  # one quintile per gene (map distinct on gene)
df_counts <- primary_long %>% count(mean_quintiles, group)

# Two-panel balance figure (2026-10-02 revision): left = mean-expression-quintile
# distribution of the FULL low (first) / high (fifth) dispersion sets BEFORE mean
# matching (plot_long, all members); right = the matched sets (primary seed).
# The "never" group is omitted from both panels' bars and legend (its matching
# above is unchanged -- it still sets n_sample per bin).
plot_groups <- c("fifth", "first")
df_counts_before <- plot_long %>%
  filter(quintile_group %in% c("fifth_quint_var", "first_quint_var")) %>%
  mutate(group = ifelse(quintile_group == "first_quint_var", "first", "fifth")) %>%
  count(mean_quintiles, group)
df_counts_after <- df_counts %>% filter(group %in% plot_groups)
log_step("balance plot before-matching counts (fifth + first genes)", sum(df_counts_before$n))
log_step("balance plot after-matching counts (fifth + first genes)", sum(df_counts_after$n))

balance_panel <- function(d, title) {
  ggplot(d, aes(x = factor(mean_quintiles), y = n, fill = group)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c(fifth = "steelblue", first = "firebrick"), breaks = plot_groups) +
    labs(x = "Mean-expression quintile", y = "Number of genes", fill = "Group", title = title) +
    theme_minimal()
}
balance_plot <- (balance_panel(df_counts_before, "Before mean matching") +
                 balance_panel(df_counts_after,
                               sprintf("Mean-matched (seed %d)", match_seed_primary))) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")

ggsave(file.path(mm_dir, "balance_check_barplot.pdf"),
       balance_plot, width = 9, height = 3.5, units = "in")

# ---- Save -------------------------------------------------------------------
saveRDS(list(matched_by_seed    = matched_by_seed,
             match_seed_primary = match_seed_primary,
             plot_long          = plot_long,       # carries mean-bin membership
             global_mean_ranks  = global_mean_ranks_df),
        file.path(mm_dir, "matched_gene_sets.rds"))
for (s in names(matched_by_seed)) {
  writeLines(matched_by_seed[[s]]$first, file.path(mm_dir, sprintf("matched_first_seed%s.txt", s)))
  writeLines(matched_by_seed[[s]]$fifth, file.path(mm_dir, sprintf("matched_fifth_seed%s.txt", s)))
  writeLines(matched_by_seed[[s]]$never, file.path(mm_dir, sprintf("matched_never_seed%s.txt", s)))
}
log_step(paste("wrote matched gene sets + balance plot to", mm_dir))
