#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

# ---- Named constants ---------------------------------------------------------
repo_dir     <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
analysis_dir <- file.path(repo_dir, "analysis/subcluster_dispersion_robustness")

# Both step-02 configs Brendan asked for (2026-09-02) -- matches script 03's
# SUBCLUSTER_CONFIG_LABEL values, both already run.
CONFIG_LABELS <- c("mincells1000_minind2", "mincells700_minind2")

# Distinct, colourblind-safe qualitative colour per resolution (ColorBrewer
# Dark2, first 3) -- no pre-existing per-resolution palette convention in
# .living/conventions.md (that file's celltype_colors is cell-type-specific),
# so this is a new but standard qualitative choice.
resolution_colors <- c("0.5" = "#1b9e77", "1.0" = "#d95f02", "2.0" = "#7570b3")

source(file.path(analysis_dir, "scripts", "utils.R"))
log_msg <- make_logger("04_combined_vs_original_rho")

for (config_label in CONFIG_LABELS) {
  out_dir <- file.path(analysis_dir, "outputs", config_label)
  vs_csv <- file.path(out_dir, "subcluster_vs_original_correlation.csv")
  if (!file.exists(vs_csv)) {
    stop(sprintf("Missing %s -- run 03_robustness_correlations.R with SUBCLUSTER_CONFIG_LABEL=%s first.",
                  vs_csv, config_label))
  }

  # colClasses forces "resolution" to stay character -- read.csv's default
  # type-sniffing would otherwise parse "0.5"/"1.0"/"2.0" as numeric and
  # silently drop the trailing zero from "1.0"/"2.0", breaking the exact
  # string match against resolution_colors' names below.
  vs_df <- read.csv(vs_csv, stringsAsFactors = FALSE,
                     colClasses = c(resolution = "character"))
  stopifnot(
    "subcluster_vs_original_correlation.csv must have the 5 expected columns" =
      identical(colnames(vs_df), c("resolution", "subcluster", "parent_ct", "n_genes", "rho")),
    "Every resolution value must have an assigned colour" =
      all(unique(vs_df$resolution) %in% names(resolution_colors))
  )
  vs_df$resolution <- factor(vs_df$resolution, levels = names(resolution_colors))

  log_msg("[%s] %d (cell type x resolution) rows across %d cell types, %d resolutions.",
          config_label, nrow(vs_df), length(unique(vs_df$parent_ct)),
          length(unique(vs_df$resolution)))
  present <- table(vs_df$parent_ct, vs_df$resolution)
  if (any(present == 0)) {
    missing_combos <- which(present == 0, arr.ind = TRUE)
    for (i in seq_len(nrow(missing_combos))) {
      log_msg("  [no sub-clusters -- box omitted] %s x res %s",
              rownames(present)[missing_combos[i, "row"]],
              colnames(present)[missing_combos[i, "col"]])
    }
  }

  p <- ggplot(vs_df, aes(x = parent_ct, y = rho, fill = resolution)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge2(preserve = "single")) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
               size = 1.2, alpha = 0.7) +
    scale_fill_manual(values = resolution_colors) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = sprintf("Sub-cluster vs original cell-type dispersion, all resolutions (%s)", config_label),
         x = "Parent cell type", y = "Spearman rho (sub-cluster vs original)", fill = "Resolution")

  # PDF (not PNG) per Brendan's request, 2026-09-02 -- vector output for
  # editing in Affinity Designer; cairo_pdf matches the repo's device
  # convention for vector figures elsewhere (.living/conventions.md).
  out_file <- file.path(out_dir, "subcluster_vs_original_rho_combined.pdf")
  ggsave(out_file, p, width = 10, height = 6, device = cairo_pdf)
  log_msg("  Saved: %s", out_file)
}

log_msg("DONE.")
