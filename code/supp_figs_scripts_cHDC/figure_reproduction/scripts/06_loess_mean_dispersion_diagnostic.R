#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
fig_repro_dir <- file.path(repo_dir, "figure_reproduction")
filt_mat_dir <- file.path(fig_repro_dir, "outputs", "filtered_matrices")
out_dir <- file.path(fig_repro_dir, "outputs", "loess_mt_diagnostic")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(fig_repro_dir, "scripts", "utils.R"))

# ANALYSIS_OK[experimental-design]: the 3 pooled iPSC lines are the fixed
# experiment (data/DATA_MANIFEST.md chdc-scrnaseq), matching 01's identical
# constant; asserted against the retained-matrix columns below.
individuals <- c("NA19099", "NA19093", "NA18522")
cell_types <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)

# ANALYSIS_OK[filtering]: ^MT- mitochondrial-gene identification pattern --
# used here to LABEL genes for comparison, not to exclude them (this script's
# purpose is testing whether MT genes should stay in); retained MT counts
# logged per cell type below.
mito_pattern <- "^MT-"
# ANALYSIS_OK[model-formula]: span=0.5 matches Alexander Chen's
# dispersion_diagnostics.Rmd Section 1 method -- kept identical here so the
# trend-flatness comparison is apples-to-apples.
loess_span <- 0.5
min_mt_for_test <- 2  # Wilcoxon needs >=2 obs per group to be meaningful

residual_rows <- list()
curve_shift_rows <- list()
mt_summary_rows <- list()
gene_level_rows <- list()  # full per-gene/per-individual table, for downstream follow-up diagnostics

for (ct in cell_types) {
  mean_path <- file.path(filt_mat_dir, paste0(ct, "_mean_retained.csv"))
  disp_path <- file.path(filt_mat_dir, paste0(ct, "_dispersion_retained.csv"))
  if (!file.exists(mean_path) || !file.exists(disp_path)) {
    stop("Missing filtered matrix for ", ct,
         " -- run 01_exclude_topexpressed_genes.R first: ", mean_path, " / ", disp_path)
  }
  mean_mat <- as.matrix(read.csv(mean_path, row.names = 1, check.names = FALSE))
  disp_mat <- as.matrix(read.csv(disp_path, row.names = 1, check.names = FALSE))

  stopifnot(
    "mean/dispersion retained matrices must have identical gene rows in identical order" =
      identical(rownames(mean_mat), rownames(disp_mat)),
    "retained matrices must have exactly the 3-individual columns for this cell type, in order" =
      identical(colnames(mean_mat), ct_columns(ct)) && identical(colnames(disp_mat), ct_columns(ct)),
    "retained matrices must contain no non-finite values" =
      all(is.finite(mean_mat)) && all(is.finite(disp_mat))
  )

  genes <- rownames(mean_mat)
  is_mt <- grepl(mito_pattern, genes)
  n_mt <- sum(is_mt)
  message(sprintf("[%s] %d retained genes, %d mitochondrial (%s)",
                  ct, length(genes), n_mt, paste(genes[is_mt], collapse = ", ")))

  long_df <- do.call(rbind, lapply(individuals, function(ind) {
    col <- paste0(ind, "_", ct)
    data.frame(gene = genes, individual = ind,
               mean_expr = mean_mat[, col], dispersion = disp_mat[, col],
               is_mt = is_mt, stringsAsFactors = FALSE)
  }))
  long_df$individual <- factor(long_df$individual, levels = individuals)

  # --- per-individual LOESS fit (MT-kept, the actual analysis config), plus
  # the MT-dropped sensitivity fit for the curve-shift comparison -----------
  panel_dfs <- list()
  for (ind in individuals) {
    ind_df <- long_df[long_df$individual == ind, ]
    fit_all <- loess(dispersion ~ mean_expr, data = ind_df, span = loess_span)
    ind_df$fitted_all <- predict(fit_all, newdata = ind_df)
    ind_df$resid_all <- ind_df$dispersion - ind_df$fitted_all

    noMT_df <- ind_df[!ind_df$is_mt, ]
    fit_noMT <- loess(dispersion ~ mean_expr, data = noMT_df, span = loess_span)
    grid_x <- seq(min(noMT_df$mean_expr), max(noMT_df$mean_expr), length.out = 200)
    pred_all <- predict(fit_all, newdata = data.frame(mean_expr = grid_x))
    pred_noMT <- predict(fit_noMT, newdata = data.frame(mean_expr = grid_x))
    valid <- is.finite(pred_all) & is.finite(pred_noMT)
    curve_shift_rows[[paste(ct, ind)]] <- data.frame(
      cell_type = ct, individual = ind, n_grid_points = sum(valid),
      curve_correlation = cor(pred_all[valid], pred_noMT[valid]),
      curve_max_abs_diff = max(abs(pred_all[valid] - pred_noMT[valid])),
      curve_mean_abs_diff = mean(abs(pred_all[valid] - pred_noMT[valid])),
      stringsAsFactors = FALSE)

    panel_dfs[[ind]] <- ind_df
  }
  long_df <- do.call(rbind, panel_dfs)

  gene_level_rows[[ct]] <- data.frame(
    cell_type = ct, gene = long_df$gene, individual = as.character(long_df$individual),
    mean_expr = long_df$mean_expr, dispersion = long_df$dispersion, is_mt = long_df$is_mt,
    fitted_all = long_df$fitted_all, resid_all = long_df$resid_all, stringsAsFactors = FALSE)

  # --- pooled (across 3 individuals) residual test: MT vs. non-MT ----------
  mt_resid <- long_df$resid_all[long_df$is_mt]
  nonmt_resid <- long_df$resid_all[!long_df$is_mt]
  if (length(mt_resid) >= min_mt_for_test) {
    wtest_p <- wilcox.test(mt_resid, nonmt_resid)$p.value
  } else {
    message(sprintf("  [warn] only %d MT residuals pooled for %s -- skipping Wilcoxon test",
                    length(mt_resid), ct))
    wtest_p <- NA_real_
  }
  residual_rows[[ct]] <- data.frame(
    cell_type = ct, n_mt_retained = n_mt,
    n_mt_residuals_pooled = length(mt_resid), n_nonmt_residuals_pooled = length(nonmt_resid),
    median_mt_residual = median(mt_resid), median_nonmt_residual = median(nonmt_resid),
    wilcox_p = wtest_p, stringsAsFactors = FALSE)

  mt_summary_rows[[ct]] <- data.frame(
    cell_type = ct, n_retained = length(genes), n_mt_retained = n_mt,
    mt_genes = paste(sort(genes[is_mt]), collapse = ";"), stringsAsFactors = FALSE)

  # --- plot: per-individual scatter + LOESS (MT-kept fit), MT highlighted --
  p <- ggplot(long_df, aes(x = mean_expr, y = dispersion)) +
    geom_point(data = long_df[!long_df$is_mt, ], alpha = 0.12, size = 0.3, colour = "grey30") +
    geom_smooth(method = "loess", span = loess_span, colour = "#4575b4",
                linewidth = 0.4, se = TRUE) +
    geom_point(data = long_df[long_df$is_mt, ], colour = "#d73027", size = 1.1) +
    facet_wrap(~ individual, ncol = 3) +
    labs(x = "Log mean expression", y = "Dispersion",
         title = paste0(gsub("_", " ", ct),
                        " -- dispersion vs. mean (ribo-dropped, -6-excluded, MT kept)"),
         subtitle = sprintf("%d retained genes (%d mitochondrial, red)", length(genes), n_mt)) +
    theme_classic(base_family = "Arial") +
    theme(
      text = element_text(size = 6), axis.text = element_text(size = 6),
      axis.title = element_text(size = 6), plot.title = element_text(size = 6),
      plot.subtitle = element_text(size = 6), strip.text = element_text(size = 6),
      legend.text = element_text(size = 6), legend.title = element_text(size = 6),
      axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
      plot.margin = margin(2, 2, 2, 2, unit = "mm"))

  out_file <- file.path(out_dir, paste0(ct, "_loess_mt_diagnostic.pdf"))
  ggsave(out_file, p, width = 3 * 1.8, height = 1.8 + 0.4, units = "in",
         device = cairo_pdf, limitsize = FALSE)
  message("  Saved: ", out_file)
}

residual_df <- do.call(rbind, residual_rows); rownames(residual_df) <- NULL
curve_shift_df <- do.call(rbind, curve_shift_rows); rownames(curve_shift_df) <- NULL
mt_summary_df <- do.call(rbind, mt_summary_rows); rownames(mt_summary_df) <- NULL
gene_level_df <- do.call(rbind, gene_level_rows); rownames(gene_level_df) <- NULL

write.csv(residual_df, file.path(out_dir, "mt_vs_nonmt_residual_test.csv"), row.names = FALSE)
write.csv(curve_shift_df, file.path(out_dir, "mt_kept_vs_dropped_curve_shift.csv"), row.names = FALSE)
write.csv(mt_summary_df, file.path(out_dir, "mt_gene_retention_summary.csv"), row.names = FALSE)
# Full per-gene/per-individual mean, dispersion, LOESS-fitted value, and
# residual -- this is what script 21 (Fig 1) reads directly.
write.csv(gene_level_df, file.path(out_dir, "mt_gene_level_residuals.csv"), row.names = FALSE)

message("\n[done] figures + tables in ", out_dir)
print(residual_df)
print(curve_shift_df)

# ---- reportable values (kept for parity with the source script; harmless
# in this folder even though figure_reproduction/ has no report-generator
# scaffolding) -----------------------------------------------------------
reportable <- do.call(rbind, lapply(seq_len(nrow(residual_df)), function(i) {
  row <- residual_df[i, ]
  data.frame(
    key = c(paste0("n_mt_retained_", row$cell_type),
            paste0("wilcox_p_mt_vs_nonmt_residual_", row$cell_type)),
    value = c(row$n_mt_retained, signif(row$wilcox_p, 4)),
    provenance = "figure_reproduction/scripts/06_loess_mean_dispersion_diagnostic.R",
    stringsAsFactors = FALSE)
}))
write.csv(reportable, file.path(fig_repro_dir, "outputs", "reportable_values_06.csv"), row.names = FALSE)
