#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)              # topTable
  library(variancePartition)  # S4 method dispatch for the dream/eBayes fit object
  library(ggplot2)
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
fig_repro_dir <- file.path(repo_dir, "figure_reproduction")
filt_mat_dir <- file.path(fig_repro_dir, "outputs", "filtered_matrices")
refit_dir <- file.path(fig_repro_dir, "outputs", "dd_refit")
out_dir <- file.path(fig_repro_dir, "outputs", "dd_trend_correction")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(fig_repro_dir, "scripts", "utils.R"))

# ANALYSIS_OK[experimental-design]: the 3 pooled iPSC lines are the fixed
# experiment (data/DATA_MANIFEST.md chdc-scrnaseq), matching 01's/03's
# identical constant; asserted against the filtered-matrix columns below.
individuals <- c("NA19099", "NA19093", "NA18522")
cell_types <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)
# ANALYSIS_OK[model-formula]: span=0.5 matches 06's identical parameter --
# kept the same so the trend estimates used here are consistent with the
# rest of this diagnostic family.
loess_span <- 0.5
dd_adjp_thresh <- 0.05  # matches 03's DD-call threshold, for the reference flag only

fit_path <- file.path(refit_dir, "dream_fit_retained.rds")
if (!file.exists(fit_path)) {
  stop("Missing ", fit_path, " -- run 03_rerun_dd_rank_analyses.R Part B first")
}
fit_retained <- readRDS(fit_path)

# 15 pairwise contrasts named compare_{B}_{A} -- GENERATED identically to
# script 03 (combn + same paste0 logic), not parsed from strings, since
# several cell-type names (Cardiac_fibroblasts, Activated_fibroblasts,
# Foregut_endodermal_cells) contain underscores and would make string-
# splitting on "_" ambiguous.
pairs <- combn(cell_types, 2, simplify = FALSE)
contrast_names <- vapply(pairs, function(x) paste0("compare_", x[2], "_", x[1]), character(1))
contrast_pair_a <- setNames(vapply(pairs, function(x) x[1], character(1)), contrast_names)
contrast_pair_b <- setNames(vapply(pairs, function(x) x[2], character(1)), contrast_names)
# NOTE: fit_retained$coefficients also retains the 6 underlying cell_type
# fixed-effect coefficients alongside the 15 contrast coefficients (a
# consequence of how makeContrastsDream/eBayes preserves the model design) --
# so we check the 15 generated names are a SUBSET, not a set-equal match.
stopifnot(
  "generated contrast names must all be present among dream_fit_retained's coefficient names" =
    all(contrast_names %in% colnames(fit_retained$coefficients)))
message(sprintf("[load] dream_fit_retained.rds loaded; %d genes x %d contrasts, verified against generated names",
                nrow(fit_retained$coefficients), length(contrast_names)))

# --------------------------- pooled per-cell-type trends ---------------------

pooled_fits_gaussian <- list()
pooled_fits_robust <- list()
gene_mean_by_ct <- list()  # per cell type: named vector, gene -> individual-averaged mean_expr

for (ct in cell_types) {
  mean_path <- file.path(filt_mat_dir, paste0(ct, "_mean_retained.csv"))
  disp_path <- file.path(filt_mat_dir, paste0(ct, "_dispersion_retained.csv"))
  if (!file.exists(mean_path) || !file.exists(disp_path)) {
    stop("Missing filtered matrix for ", ct, ": ", mean_path, " / ", disp_path)
  }
  mean_mat <- as.matrix(read.csv(mean_path, row.names = 1, check.names = FALSE))
  disp_mat <- as.matrix(read.csv(disp_path, row.names = 1, check.names = FALSE))
  stopifnot(
    "mean/dispersion retained matrices must have identical gene rows in identical order" =
      identical(rownames(mean_mat), rownames(disp_mat)),
    "retained matrices must have exactly the 3-individual columns for this cell type, in order" =
      identical(colnames(mean_mat), ct_columns(ct)) && identical(colnames(disp_mat), ct_columns(ct)),
    "retained matrices must contain no non-finite values" =
      all(is.finite(mean_mat)) && all(is.finite(disp_mat)))

  genes <- rownames(mean_mat)
  pooled_df <- do.call(rbind, lapply(individuals, function(ind) {
    col <- paste0(ind, "_", ct)
    data.frame(mean_expr = mean_mat[, col], dispersion = disp_mat[, col])
  }))
  pooled_fits_gaussian[[ct]] <- loess(dispersion ~ mean_expr, data = pooled_df, span = loess_span)
  pooled_fits_robust[[ct]] <- loess(dispersion ~ mean_expr, data = pooled_df, span = loess_span,
                                    family = "symmetric")
  message(sprintf("[%s] pooled trend fit on %d points (%d genes x 3 individuals)",
                  ct, nrow(pooled_df), length(genes)))

  gene_mean_by_ct[[ct]] <- setNames(rowMeans(mean_mat), genes)
}

# --------------------------- per-contrast topTable + correction --------------

union_genes <- rownames(fit_retained$coefficients)
message(sprintf("[union] %d genes in dream_fit_retained's tested universe", length(union_genes)))
for (ct in cell_types) {
  stopifnot("every union-retained gene must have a mean-expression value in every cell type" =
              all(union_genes %in% names(gene_mean_by_ct[[ct]])))
}

contrast_rows <- list()
for (cn in contrast_names) {
  ct_a <- contrast_pair_a[[cn]]
  ct_b <- contrast_pair_b[[cn]]

  tt <- topTable(fit_retained, coef = cn, number = Inf, sort.by = "none")
  stopifnot("topTable row names must exactly match the union gene universe, in the same order" =
              identical(rownames(tt), union_genes))

  mean_a <- gene_mean_by_ct[[ct_a]][union_genes]
  mean_b <- gene_mean_by_ct[[ct_b]][union_genes]

  pred_a_gauss <- predict(pooled_fits_gaussian[[ct_a]], newdata = data.frame(mean_expr = mean_a))
  pred_b_gauss <- predict(pooled_fits_gaussian[[ct_b]], newdata = data.frame(mean_expr = mean_b))
  pred_a_robust <- predict(pooled_fits_robust[[ct_a]], newdata = data.frame(mean_expr = mean_a))
  pred_b_robust <- predict(pooled_fits_robust[[ct_b]], newdata = data.frame(mean_expr = mean_b))

  expected_artifact_gaussian <- pred_b_gauss - pred_a_gauss
  expected_artifact_robust <- pred_b_robust - pred_a_robust

  contrast_rows[[cn]] <- data.frame(
    contrast = cn, gene = union_genes, cell_type_A = ct_a, cell_type_B = ct_b,
    mean_expr_A = mean_a, mean_expr_B = mean_b,
    observed_logFC = tt$logFC, adj.P.Val = tt$adj.P.Val,
    dd_call_current = !is.na(tt$adj.P.Val) & tt$adj.P.Val <= dd_adjp_thresh,
    expected_artifact_gaussian = expected_artifact_gaussian,
    expected_artifact_robust = expected_artifact_robust,
    trend_corrected_logfc_gaussian = tt$logFC - expected_artifact_gaussian,
    trend_corrected_logfc_robust = tt$logFC - expected_artifact_robust,
    stringsAsFactors = FALSE)
  message(sprintf("  [%s] %s vs %s: %d genes, %d current DD calls",
                  cn, ct_b, ct_a, length(union_genes), sum(contrast_rows[[cn]]$dd_call_current)))
}

effects_df <- do.call(rbind, contrast_rows); rownames(effects_df) <- NULL
stopifnot("output table must have exactly contrasts x genes rows" =
            nrow(effects_df) == length(contrast_names) * length(union_genes))
write.csv(effects_df, file.path(out_dir, "dd_trend_corrected_effects.csv"), row.names = FALSE)
message(sprintf("[done] wrote %d rows (%d contrasts x %d genes) to %s",
                nrow(effects_df), length(contrast_names), length(union_genes),
                file.path(out_dir, "dd_trend_corrected_effects.csv")))

# --------------------------------- summary ------------------------------------

summary_df <- do.call(rbind, lapply(contrast_names, function(cn) {
  sub <- effects_df[effects_df$contrast == cn, ]
  data.frame(
    contrast = cn, n_genes = nrow(sub), n_dd_call_current = sum(sub$dd_call_current),
    median_expected_artifact_gaussian = median(sub$expected_artifact_gaussian),
    iqr_expected_artifact_gaussian = IQR(sub$expected_artifact_gaussian),
    median_expected_artifact_robust = median(sub$expected_artifact_robust),
    iqr_expected_artifact_robust = IQR(sub$expected_artifact_robust),
    stringsAsFactors = FALSE)
}))
write.csv(summary_df, file.path(out_dir, "dd_trend_corrected_effects_summary_by_contrast.csv"),
          row.names = FALSE)
message("\n[summary by contrast]")
print(summary_df, digits = 3)

# ----------------------------------- plot -------------------------------------

p <- ggplot(effects_df, aes(x = expected_artifact_gaussian, y = observed_logFC,
                            colour = dd_call_current)) +
  geom_point(alpha = 0.15, size = 0.3) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50", linewidth = 0.3, linetype = "dashed") +
  scale_colour_manual(name = "Current DD call\n(adj.P.Val<=0.05)",
                      values = c("TRUE" = "#d73027", "FALSE" = "grey60")) +
  facet_wrap(~ contrast, ncol = 5) +
  labs(x = "Expected artifact (gaussian pooled trend)", y = "Observed logFC",
       title = "Observed DD effect vs. mean-expression-trend-predicted artifact, per contrast",
       subtitle = "Dashed line: observed = expected (i.e. entirely trend-explained). MT-kept, ribo-dropped, -6-excluded gene set.") +
  theme_classic(base_family = "Arial") +
  theme(
    text = element_text(size = 6), axis.text = element_text(size = 6),
    axis.title = element_text(size = 6), plot.title = element_text(size = 7),
    plot.subtitle = element_text(size = 6), strip.text = element_text(size = 6),
    legend.text = element_text(size = 6), legend.title = element_text(size = 6),
    legend.key.size = unit(6, "pt"),
    axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
    plot.margin = margin(2, 2, 2, 2, unit = "mm"))

out_file <- file.path(out_dir, "dd_observed_vs_expected_artifact.pdf")
ggsave(out_file, p, width = 14, height = 9, units = "in", device = cairo_pdf, limitsize = FALSE)
message("  Saved: ", out_file)

# ------------------------- low-ratio tail (follow-up) ------------------------
notable_artifact_thresh <- 0.1
low_ratio_thresh <- 2
notable <- effects_df[effects_df$dd_call_current &
                        abs(effects_df$expected_artifact_gaussian) >= notable_artifact_thresh, ]
notable$ratio_obs_to_artifact <- abs(notable$observed_logFC) / abs(notable$expected_artifact_gaussian)
low_ratio_tail <- notable[notable$ratio_obs_to_artifact < low_ratio_thresh, ]
n_low_ratio_before_sort <- nrow(low_ratio_tail)
low_ratio_tail <- low_ratio_tail[order(low_ratio_tail$ratio_obs_to_artifact), ]
stopifnot("row reorder must not change row count" = nrow(low_ratio_tail) == n_low_ratio_before_sort)
write.csv(low_ratio_tail, file.path(out_dir, "dd_low_ratio_tail.csv"), row.names = FALSE)
message(sprintf(
  "\n[low-ratio tail] %d of %d notable-artifact (|artifact|>=%.1f) DD calls have ratio < %d (artifact >= half the observed effect) -- written to %s",
  nrow(low_ratio_tail), nrow(notable), notable_artifact_thresh, low_ratio_thresh,
  file.path(out_dir, "dd_low_ratio_tail.csv")))

# ---- reportable values (kept for parity with the source script) -------------
reportable <- data.frame(
  key = c("dd_trend_correction_median_expected_artifact_gaussian",
          "dd_trend_correction_median_expected_artifact_robust",
          "dd_trend_correction_n_contrasts", "dd_trend_correction_n_genes_per_contrast",
          "dd_trend_correction_n_notable_artifact_dd_calls",
          "dd_trend_correction_n_low_ratio_tail"),
  value = c(signif(median(effects_df$expected_artifact_gaussian), 4),
            signif(median(effects_df$expected_artifact_robust), 4),
            length(contrast_names), length(union_genes),
            nrow(notable), nrow(low_ratio_tail)),
  provenance = "figure_reproduction/scripts/14_dd_trend_corrected_effects.R",
  stringsAsFactors = FALSE)
write.csv(reportable, file.path(fig_repro_dir, "outputs", "reportable_values_14.csv"), row.names = FALSE)

message("\n[done] all outputs in ", out_dir)
