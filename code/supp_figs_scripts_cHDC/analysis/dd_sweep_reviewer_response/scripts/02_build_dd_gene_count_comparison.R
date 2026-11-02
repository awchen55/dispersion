library(variancePartition)
library(limma)

source("scripts/00_config.R")

OUT_DIR <- "outputs"
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)

n_sig_genes <- function(fit, contrast) {
  tt <- topTable(fit, coef = contrast, number = Inf, adjust.method = "BH", sort.by = "none")
  sum(!is.na(tt$adj.P.Val) & tt$adj.P.Val < 0.05)
}

rows <- list()

for (thresh in THRESHOLDS) {
  brendan_fit <- readRDS(file.path(BRENDAN_DD_FIT_DIR, paste0("dd_fit_cell_count_", thresh, ".rds")))

  for (i in seq_len(nrow(SHARED_CONTRASTS))) {
    contrast   <- SHARED_CONTRASTS$contrast[i]
    erik_batch <- SHARED_CONTRASTS$erik_batch[i]
    erik_fit   <- readRDS(file.path(erik_dd_fit_dir(erik_batch), paste0("dd_fit_cell_count_", thresh, ".rds")))

    stopifnot(contrast %in% colnames(brendan_fit$coefficients))
    stopifnot(contrast %in% colnames(erik_fit$coefficients))

    n_brendan <- n_sig_genes(brendan_fit, contrast)
    n_erik    <- n_sig_genes(erik_fit, contrast)

    rows[[length(rows) + 1]] <- data.frame(
      cell_count_threshold = thresh,
      contrast             = contrast,
      erik_batch            = toupper(erik_batch),
      n_sig_dd_genes_brendan = n_brendan,
      n_sig_dd_genes_erik    = n_erik,
      more_in_brendan        = n_brendan - n_erik,
      fold_more_in_brendan   = ifelse(n_erik == 0, NA, round(n_brendan / n_erik, 2)),
      stringsAsFactors = FALSE
    )
  }
}

df <- do.call(rbind, rows)
df$contrast_label <- gsub("_", " ", sub("^compare_", "", df$contrast))
df <- df[order(-df$cell_count_threshold, df$contrast), ]

out_csv <- file.path(OUT_DIR, "dd_gene_count_comparison_brendan_vs_erik.csv")
write.csv(df, out_csv, row.names = FALSE)
message("Wrote: ", out_csv)

# Summary rolled up by threshold (sum across the 4 shared contrasts) --
# the single "how many more DD genes at each downsampling level" number.
summary_df <- aggregate(
  cbind(n_sig_dd_genes_brendan, n_sig_dd_genes_erik, more_in_brendan) ~ cell_count_threshold,
  data = df, FUN = sum
)
summary_df <- summary_df[order(-summary_df$cell_count_threshold), ]
summary_df$fold_more_in_brendan <- round(
  summary_df$n_sig_dd_genes_brendan / summary_df$n_sig_dd_genes_erik, 2
)

out_summary_csv <- file.path(OUT_DIR, "dd_gene_count_comparison_summary_by_threshold.csv")
write.csv(summary_df, out_summary_csv, row.names = FALSE)
message("Wrote: ", out_summary_csv)

print(summary_df)

# Brendan full-depth (no downsampling) reference point, per shared contrast --
# added per Brendan's request (2026-09-03) as one more x-axis column in the
# boxplot below. Uses script 01's already-fit no-random-effect reference
# (same model spec as every downsampled-sweep fit: ~ 0 + cell_type). Erik has
# no full-depth counterpart in this comparison, so this column carries only a
# Brendan box.
brendan_full_fit <- readRDS(file.path(
  "outputs", "brendan_full_reference_no_random_effect", "dd_fit_brendan_full_no_random_effect.rds"
))
full_rows <- list()
for (i in seq_len(nrow(SHARED_CONTRASTS))) {
  contrast <- SHARED_CONTRASTS$contrast[i]
  stopifnot(contrast %in% colnames(brendan_full_fit$coefficients))
  full_rows[[length(full_rows) + 1]] <- data.frame(
    contrast = contrast,
    contrast_label = gsub("_", " ", sub("^compare_", "", contrast)),
    n_sig_dd_genes_brendan = n_sig_genes(brendan_full_fit, contrast),
    stringsAsFactors = FALSE
  )
}
brendan_full_df <- do.call(rbind, full_rows)
out_full_csv <- file.path(OUT_DIR, "dd_gene_count_brendan_full.csv")
write.csv(brendan_full_df, out_full_csv, row.names = FALSE)
message("Wrote: ", out_full_csv)
print(brendan_full_df)

# Boxplot: significant DD genes PER SHARED CONTRAST (each of the 4 shared
# contrasts is one data point), not summed across contrasts -- per Brendan's
# request (2026-09-03), replacing the earlier summed bar chart above. X-axis
# adds a "Full" column (Brendan only, no downsampling) ahead of the 5
# downsampled thresholds.
library(ggplot2)
CELL_COUNT_LEVELS <- c(as.character(sort(THRESHOLDS, decreasing = FALSE)), "Full")
box_df <- rbind(
  data.frame(cell_count_group = "Full", cohort = "Brendan",
             n_sig_genes = brendan_full_df$n_sig_dd_genes_brendan, stringsAsFactors = FALSE),
  data.frame(cell_count_group = as.character(df$cell_count_threshold), cohort = "Brendan",
             n_sig_genes = df$n_sig_dd_genes_brendan, stringsAsFactors = FALSE),
  data.frame(cell_count_group = as.character(df$cell_count_threshold), cohort = "Erik (B2+B6)",
             n_sig_genes = df$n_sig_dd_genes_erik, stringsAsFactors = FALSE)
)
box_df$cell_count_group <- factor(box_df$cell_count_group, levels = CELL_COUNT_LEVELS)
# Legend labels per Brendan's explicit request, 2026-09-03 -- the Erik label
# is deliberately wrapped to two lines via an embedded "\n" (ggplot's default
# legend text renders embedded newlines as a line break, no extra package
# needed).
box_df$cohort <- factor(box_df$cohort, levels = c("Brendan", "Erik (B2+B6)"),
                         labels = c("n=3 cHDC cohort", "McIntire et al.\ncHDC cohort"))

p <- ggplot(box_df, aes(x = cell_count_group, y = n_sig_genes, fill = cohort)) +
  geom_boxplot(position = position_dodge2(preserve = "single"), outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
             size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = setNames(c("#882255", "#4477AA"), levels(box_df$cohort))) +
  labs(
    title = "Significant DD genes by downsampled cell count, per shared contrast",
    subtitle = "each point = 1 of the 4 shared contrasts",
    x = "Cells per sample", y = "Significant DD genes (per contrast)",
    fill = NULL
  ) +
  theme_bw()

out_fig <- file.path(OUT_DIR, "figures", "dd_gene_count_comparison_brendan_vs_erik.pdf")
ggsave(out_fig, p, width = 8, height = 5.5)
message("Wrote: ", out_fig)

# Fold-increase trend: median (across the 4 shared contrasts) of Brendan's
# per-contrast fold-more-DD-genes-than-Erik, at each downsampled threshold --
# per Brendan's request (2026-09-03). Reuses df's already-computed per-contrast
# fold_more_in_brendan (n_sig_dd_genes_brendan / n_sig_dd_genes_erik); median,
# not mean, so one contrast with an extreme fold value can't dominate. "Full"
# is excluded -- Erik has no full-depth counterpart, so no fold is defined there.
fold_median_df <- aggregate(fold_more_in_brendan ~ cell_count_threshold, data = df, FUN = median)
colnames(fold_median_df)[colnames(fold_median_df) == "fold_more_in_brendan"] <- "median_fold_more_in_brendan"
fold_median_df <- fold_median_df[order(fold_median_df$cell_count_threshold), ]

out_fold_csv <- file.path(OUT_DIR, "dd_fold_increase_median_by_threshold.csv")
write.csv(fold_median_df, out_fold_csv, row.names = FALSE)
message("Wrote: ", out_fold_csv)
print(fold_median_df)

fold_median_df$cell_count_threshold <- factor(fold_median_df$cell_count_threshold,
                                               levels = sort(THRESHOLDS, decreasing = FALSE))

p_fold <- ggplot(fold_median_df, aes(x = cell_count_threshold, y = median_fold_more_in_brendan, group = 1)) +
  geom_line(color = "#882255") +
  geom_point(color = "#882255", size = 2.5) +
  labs(
    title = "Fold increase in significant DD gene count, Brendan vs. Erik, by cell count",
    subtitle = "Median (across the 4 shared contrasts) of the per-contrast fold increase",
    x = "Cells per sample", y = "Median fold increase (n=3 cHDC cohort / McIntire et al. cHDC cohort)"
  ) +
  theme_bw()

out_fold_fig <- file.path(OUT_DIR, "figures", "dd_fold_increase_median_by_threshold.pdf")
ggsave(out_fold_fig, p_fold, width = 7, height = 5)
message("Wrote: ", out_fold_fig)
