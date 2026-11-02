stopifnot(file.exists("CLAUDE.md"))  # must run from project root

suppressPackageStartupMessages(library(ggplot2))

# Shared read_wide()/read_active_cell_types() live in utils.R (single
# canonical definition; see R025). read_wide() looks up OUT_DIR/CELL_TYPES
# from this script's environment at call time.
# PATCHED for figure_reproduction (see PORT_CHANGELOG_crosscohort.md).
source(file.path("figure_reproduction", "scripts", "cross_cohort_dispersion_corr", "utils.R"))

# ---------------------------------------------------------------------------
# CONSTANTS
# ---------------------------------------------------------------------------
# Named seed constant — bootstrap is reproducible (robust-analysis).
SEED <- 20260720L

CORR_METHOD_PRIMARY   <- "spearman"  # primary (rank-based)
CORR_METHOD_SECONDARY <- "pearson"
N_BOOT   <- 2000L                    # bootstrap resamples for the 95% CI
CI_PROBS <- c(0.025, 0.975)
MIN_PAIRS_FOR_CORR <- 3L             # minimum gene pairs to attempt a correlation

# Erik cell-count threshold — must match the value used in 01 (outputs are
# namespaced by it). See SENSITIVITY note in the header.
ERIK_CELL_COUNT_THRESHOLD <- as.integer(
  Sys.getenv("CCDC_CELL_COUNT_THRESHOLD", unset = "150")
)

# PATCHED for figure_reproduction: matches 01's OUT_ROOT.
OUT_DIR <- file.path("figure_reproduction", "outputs", "cross_cohort_dispersion_corr",
                     paste0("cc", ERIK_CELL_COUNT_THRESHOLD))
ALIGN_DIR <- file.path(OUT_DIR, "aligned")
FIG_DIR   <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Cell types that survived 01's per-CT Erik-individual-coverage filter at
# this threshold (see filter_cell_types_by_erik_coverage() in utils.R) —
# read back rather than hardcoded, so re-running at a different
# CCDC_CELL_COUNT_THRESHOLD needs no script edits.
CELL_TYPES <- read_active_cell_types(OUT_DIR)

# Display names / colours consistent with existing repo figures.
CT_DISPLAY <- c(
  Activated_fibroblasts    = "Activated fibroblasts",
  Cardiac_fibroblasts      = "Cardiac fibroblasts",
  Cardiomyocytes           = "Cardiomyocytes",
  Epithelial_cells         = "Epithelial cells",
  Foregut_endodermal_cells = "Foregut endodermal cells",
  Mesenchymal_cells        = "Mesenchymal cells"
)
celltype_colors <- c(
  "Activated fibroblasts"    = "#117733",
  "Cardiac fibroblasts"      = "#4477AA",
  "Cardiomyocytes"           = "#882255",
  "Epithelial cells"         = "#AA4499",
  "Foregut endodermal cells" = "#CC6677",
  "Mesenchymal cells"        = "#DDCC77"
)

message("=== 02_correlate_cohorts.R ===")
message("Erik cell-count threshold: ", ERIK_CELL_COUNT_THRESHOLD)
message("Reading aligned tables from: ", ALIGN_DIR)
if (!dir.exists(ALIGN_DIR)) {
  stop("Aligned-table dir not found: ", ALIGN_DIR,
       " — run 01_load_and_align_dispersion.R first (same threshold).")
}

# ---------------------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------------------

# Correlation with a bootstrap percentile CI over resampled gene PAIRS.
# Returns named list: estimate, ci_lo, ci_hi, boots (raw draws -- kept, not
# discarded, so the matched-vs-mismatched-null significance test below can
# reuse them directly rather than re-bootstrapping; added 2026-09-03 per
# Brendan's request to add statistics to matched_vs_mismatched.pdf).
# ANALYSIS_OK[within-file-helper]: called for both correlation methods in the
# per-cell-type loop below; R026 cross-file rule cannot see intra-file callers.
boot_corr <- function(x, y, method, n_boot, probs) {
  stopifnot(length(x) == length(y), length(x) >= MIN_PAIRS_FOR_CORR)
  est <- cor(x, y, method = method)
  n <- length(x)
  boots <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    boots[i] <- cor(x[idx], y[idx], method = method)
  }
  ci <- quantile(boots, probs = probs, names = FALSE, na.rm = TRUE)
  list(estimate = est, ci_lo = ci[1], ci_hi = ci[2], boots = boots)
}

# Bootstrap CI on the MEAN correlation across several mismatched-CT pairings
# (anchor_ct in mat_x vs each of partner_cts in mat_y). At each bootstrap
# iteration, genes are resampled independently WITHIN each pairing (their
# shared-gene sets differ pair to pair) and the resampled rhos are then
# averaged — this targets the sampling distribution of the reported MEAN
# statistic itself, not the pooled spread of individual pairings (pooling raw
# per-pairing bootstrap draws would conflate between-pairing variability with
# resampling noise and overstate the CI width).
# ANALYSIS_OK[within-file-helper]: called twice per cell type in the main
# loop below (cross-cohort null + within-Brendan baseline); R026 cross-file
# rule cannot see intra-file callers.
mean_mismatched_corr <- function(mat_x, mat_y, anchor_ct, partner_cts, method,
                                  n_boot, probs) {
  pairs <- lapply(partner_cts, function(pct) {
    x <- mat_x[, anchor_ct]
    y <- mat_y[, pct]
    # ANALYSIS_OK[finite-filter]: keep only finite dispersion values per
    # pairing (non-finite = too sparse to estimate); the resulting
    # shared-gene count is recorded per pairing in the `detail` data.frame
    # (n_genes column) returned below and written to
    # mismatched_pairwise_detail.csv.
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    shared <- intersect(names(x), names(y))
    stopifnot(length(shared) >= MIN_PAIRS_FOR_CORR)
    list(partner_ct = pct, x = x[shared], y = y[shared], n = length(shared))
  })
  point_rhos <- vapply(pairs, function(p) cor(p$x, p$y, method = method),
                       numeric(1))
  boots <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    iter_rhos <- vapply(pairs, function(p) {
      idx <- sample.int(p$n, p$n, replace = TRUE)
      cor(p$x[idx], p$y[idx], method = method)
    }, numeric(1))
    boots[i] <- mean(iter_rhos)
  }
  ci <- quantile(boots, probs = probs, names = FALSE, na.rm = TRUE)
  list(
    mean_rho = mean(point_rhos), ci_lo = ci[1], ci_hi = ci[2], boots = boots,
    n_pairings = length(pairs),
    detail = data.frame(
      anchor_cell_type = anchor_ct,
      partner_cell_type = vapply(pairs, `[[`, character(1), "partner_ct"),
      n_genes = vapply(pairs, `[[`, integer(1), "n"),
      spearman_rho = point_rhos,
      stringsAsFactors = FALSE
    )
  )
}

# ---------------------------------------------------------------------------
# LOAD ALIGNED MATCHED TABLES + WIDE MEDIAN TABLES (for the null control)
# ---------------------------------------------------------------------------
aligned <- list()
for (ct in CELL_TYPES) {
  f <- file.path(ALIGN_DIR, paste0(ct, ".csv"))
  stopifnot(file.exists(f))
  df <- read.csv(f, check.names = FALSE)
  stopifnot(all(c("gene", "brendan_disp", "erik_disp") %in% colnames(df)))
  stopifnot(!anyNA(df), !anyDuplicated(df$gene))
  aligned[[ct]] <- df
  message("  loaded aligned ", ct, ": ", nrow(df), " shared genes")
}

# read_wide() is sourced from utils.R (see top of file). Its canonical version
# additionally asserts the 6 CELL_TYPES columns are present and returns them in
# CELL_TYPES order — behavior-preserving here, since the same assertion follows
# immediately below and downstream access is by column name.
brendan_wide <- read_wide("brendan_median_dispersion.csv")
erik_wide    <- read_wide("erik_median_dispersion.csv")
stopifnot(all(CELL_TYPES %in% colnames(brendan_wide)))
stopifnot(all(CELL_TYPES %in% colnames(erik_wide)))

# ---------------------------------------------------------------------------
# CORRELATE — matched cell types (primary) + mismatched null
# ---------------------------------------------------------------------------
set.seed(SEED)  # reproducible bootstrap for the entire run

rows <- list()
scatter_df <- list()
null_detail_rows <- list()
baseline_detail_rows <- list()
for (ct in CELL_TYPES) {
  df <- aligned[[ct]]
  sp <- boot_corr(df$brendan_disp, df$erik_disp, CORR_METHOD_PRIMARY,
                  N_BOOT, CI_PROBS)
  pe <- boot_corr(df$brendan_disp, df$erik_disp, CORR_METHOD_SECONDARY,
                  N_BOOT, CI_PROBS)

  other_cts <- setdiff(CELL_TYPES, ct)

  # Cross-cohort mismatched-CT null: Brendan ct vs Erik's 5 OTHER cell types,
  # bootstrap mean rho across all 5 pairings (replaces a single arbitrary
  # cyclic-shift partner, which could not distinguish "no cell-type-specific
  # signal" from "which one partner happened to be picked").
  null_res <- mean_mismatched_corr(brendan_wide, erik_wide, ct, other_cts,
                                   CORR_METHOD_PRIMARY, N_BOOT, CI_PROBS)
  null_detail_rows[[ct]] <- null_res$detail

  # Within-Brendan cross-cell-type baseline: same procedure, but ct vs
  # Brendan's OWN 5 other cell types (single cohort/pipeline, no cross-cohort
  # difference in the mix). Reference point for interpreting the cross-cohort
  # null above.
  baseline_res <- mean_mismatched_corr(brendan_wide, brendan_wide, ct, other_cts,
                                       CORR_METHOD_PRIMARY, N_BOOT, CI_PROBS)
  baseline_detail_rows[[ct]] <- baseline_res$detail

  # Matched-vs-mismatched-null significance test -- Welch t-test between the
  # two bootstrap distributions (same procedure as 05_matched_vs_mismatched_
  # ttest.R, reusing THESE boots directly rather than re-bootstrapping since
  # boot_corr()/mean_mismatched_corr() now return them). Added 2026-09-03 per
  # Brendan's request to add statistics to matched_vs_mismatched.pdf, which
  # previously showed only CI-overlap-by-eye. Same caveat as 05: a bootstrap
  # resampling distribution is not an iid sample, so this is an informal/
  # descriptive significance summary, not a classical hypothesis test.
  tt <- t.test(sp$boots, null_res$boots)

  message(sprintf(
    paste0("  %s: spearman=%.3f [%.3f, %.3f] (n=%d) | ",
           "cross-cohort null (mean of %d)=%.3f [%.3f, %.3f] | ",
           "within-Brendan baseline (mean of %d)=%.3f [%.3f, %.3f] | ",
           "matched-vs-null t=%.2f p=%.3g"),
    ct, sp$estimate, sp$ci_lo, sp$ci_hi, nrow(df),
    null_res$n_pairings, null_res$mean_rho, null_res$ci_lo, null_res$ci_hi,
    baseline_res$n_pairings, baseline_res$mean_rho, baseline_res$ci_lo,
    baseline_res$ci_hi, tt$statistic, tt$p.value))

  rows[[ct]] <- data.frame(
    cell_type                          = ct,
    n_genes                            = nrow(df),
    spearman_rho                       = sp$estimate,
    spearman_ci_lo                     = sp$ci_lo,
    spearman_ci_hi                     = sp$ci_hi,
    pearson_r                          = pe$estimate,
    pearson_ci_lo                      = pe$ci_lo,
    pearson_ci_hi                      = pe$ci_hi,
    null_n_pairings                    = null_res$n_pairings,
    null_mean_spearman_rho             = null_res$mean_rho,
    null_ci_lo                         = null_res$ci_lo,
    null_ci_hi                         = null_res$ci_hi,
    brendan_baseline_n_pairings        = baseline_res$n_pairings,
    brendan_baseline_mean_spearman_rho = baseline_res$mean_rho,
    brendan_baseline_ci_lo             = baseline_res$ci_lo,
    brendan_baseline_ci_hi             = baseline_res$ci_hi,
    matched_vs_null_t_statistic        = unname(tt$statistic),
    matched_vs_null_df                 = unname(tt$parameter),
    matched_vs_null_p_value            = tt$p.value,
    stringsAsFactors = FALSE
  )
  d <- df
  d$cell_type <- ct
  scatter_df[[ct]] <- d
}

summary_df <- do.call(rbind, rows)
# BH correction across the per-cell-type matched-vs-null tests (same
# rationale as 05_matched_vs_mismatched_ttest.R's identical correction: report
# effect sizes/CIs, don't hide multiple-comparisons structure).
summary_df$matched_vs_null_p_value_bh <- p.adjust(summary_df$matched_vs_null_p_value, method = "BH")
summary_csv <- file.path(OUT_DIR, "cohort_correlation_summary.csv")
write.csv(summary_df, summary_csv, row.names = FALSE)
message("Wrote ", summary_csv)

# Per-pairing detail: every individual mismatched comparison feeding the
# aggregate means above, so those means are auditable rather than a black box.
null_detail_df <- do.call(rbind, null_detail_rows)
null_detail_df$comparison <- "cross_cohort_null"
baseline_detail_df <- do.call(rbind, baseline_detail_rows)
baseline_detail_df$comparison <- "within_brendan_baseline"
mismatched_detail <- rbind(null_detail_df, baseline_detail_df)
mismatched_detail_csv <- file.path(OUT_DIR, "mismatched_pairwise_detail.csv")
write.csv(mismatched_detail, mismatched_detail_csv, row.names = FALSE)
message("Wrote ", mismatched_detail_csv)

# ---------------------------------------------------------------------------
# SCATTER PLOTS — one panel per matched cell type
# ---------------------------------------------------------------------------
scatter_all <- do.call(rbind, scatter_df)
scatter_all$cell_type_display <- factor(
  CT_DISPLAY[scatter_all$cell_type],
  levels = names(celltype_colors)
)
lab_df <- summary_df
lab_df$cell_type_display <- factor(CT_DISPLAY[lab_df$cell_type],
                                   levels = names(celltype_colors))
lab_df$label <- sprintf("rho==%.2f", lab_df$spearman_rho)

p_scatter <- ggplot(scatter_all,
                    aes(x = brendan_disp, y = erik_disp,
                        colour = cell_type_display)) +
  geom_point(size = 0.4, alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.4, colour = "grey20") +
  geom_text(data = lab_df, parse = TRUE,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.15, vjust = 1.4, size = 3, colour = "black",
            inherit.aes = FALSE) +
  facet_wrap(~ cell_type_display, scales = "free") +
  scale_colour_manual(values = celltype_colors) +
  labs(x = "Discovery cohort dispersion (median across individuals)",
       y = "Erik cohort dispersion (median across individuals)") +
  theme_classic(base_size = 9) +
  theme(legend.position = "none",
        strip.text = element_text(size = 7),
        panel.border = element_rect(colour = "black", fill = NA,
                                    linewidth = 0.4))

scatter_pdf <- file.path(FIG_DIR, "cohort_scatter.pdf")
ggsave(scatter_pdf, p_scatter, width = 7.5, height = 5, device = cairo_pdf)
message("Wrote ", scatter_pdf)

# ---------------------------------------------------------------------------
# NULL CONTROL FIGURE — matched vs cross-cohort null vs within-Brendan
# baseline Spearman rho per cell type, each with its bootstrap 95% CI
# ---------------------------------------------------------------------------
null_long <- rbind(
  data.frame(cell_type = summary_df$cell_type, rho = summary_df$spearman_rho,
             ci_lo = summary_df$spearman_ci_lo, ci_hi = summary_df$spearman_ci_hi,
             kind = "matched CT (cross-cohort)"),
  data.frame(cell_type = summary_df$cell_type,
             rho = summary_df$null_mean_spearman_rho,
             ci_lo = summary_df$null_ci_lo, ci_hi = summary_df$null_ci_hi,
             kind = "mismatched CT (cross-cohort null)")
)
null_long$cell_type_display <- factor(CT_DISPLAY[null_long$cell_type],
                                      levels = names(celltype_colors))
null_long$kind <- factor(null_long$kind, levels = c(
  "matched CT (cross-cohort)",
  "mismatched CT (cross-cohort null)"
))

# Per-cell-type matched-vs-null p-value label (Welch t-test on the bootstrap
# distributions, computed above) -- Brendan's request, 2026-09-03. Positioned
# just above that cell type's taller CI whisker; y-axis expanded on the top
# so labels never get clipped by the panel border.
lab_df <- summary_df
lab_df$cell_type_display <- factor(CT_DISPLAY[lab_df$cell_type],
                                   levels = names(celltype_colors))
plot_range <- max(null_long$ci_hi) - min(null_long$ci_lo)
lab_df$y_pos <- pmax(lab_df$spearman_ci_hi, lab_df$null_ci_hi) + 0.09 * plot_range
# ANALYSIS_OK[display-precision]: 2e-16 is R's double-precision floor for a
# p-value this small to print meaningfully, not a significance threshold.
lab_df$p_label <- ifelse(lab_df$matched_vs_null_p_value < 2e-16, "p < 2e-16",
                         paste0("p = ", formatC(lab_df$matched_vs_null_p_value,
                                                 format = "e", digits = 2)))

p_null <- ggplot(null_long,
                 aes(x = cell_type_display, y = rho, fill = kind)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                position = position_dodge(width = 0.75), width = 0.2,
                linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_text(data = lab_df, aes(x = cell_type_display, y = y_pos, label = p_label),
            inherit.aes = FALSE, size = 2.2, family = "Arial") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.17))) +
  scale_fill_manual(values = c(
    "matched CT (cross-cohort)" = "#4477AA",
    "mismatched CT (cross-cohort null)" = "grey70"
  )) +
  labs(x = NULL, y = expression("Spearman " * rho),
       fill = NULL,
       caption = "p-values: Welch t-test, matched vs. mismatched-null bootstrap distributions (BH-adjusted values in cohort_correlation_summary.csv)") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top", legend.text = element_text(size = 7),
        plot.caption = element_text(size = 5, hjust = 0))

null_pdf <- file.path(FIG_DIR, "matched_vs_mismatched.pdf")
# 6.5 x 4.5in -- comfortably fits within a standard Letter/A4 page margin
# width, confirmed unchanged from before (adding the p-value labels/caption
# did not require growing the canvas -- expansion above gives them room
# within the existing size).
ggsave(null_pdf, p_null, width = 6.5, height = 4.5, device = cairo_pdf)
message("Wrote ", null_pdf)
message("=== 02 complete ===")
