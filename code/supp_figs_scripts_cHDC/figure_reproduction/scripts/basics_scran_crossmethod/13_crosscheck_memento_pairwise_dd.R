stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(variancePartition)  # topTable (via limma)
  library(ggplot2)
})

# ---------------------------------------------------------------------------
# NAMED CONSTANTS
# ---------------------------------------------------------------------------
CELL_TYPES <- c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                 "Epithelial_cells", "Activated_fibroblasts",
                 "Foregut_endodermal_cells")

DD_ADJP_THRESH <- 0.05   # Memento DD call, matching this repo's established threshold
                         # (dd_excl_topexpressed/scripts/03's dd_adjp_thresh)
MIN_SHARED_GENES   <- 1000L  # canary: pairs share ~11,600+ tested genes; far
                              # fewer would mean gene-name alignment broke
MIN_PAIRS_SPEARMAN <- 3L     # need >= 3 paired points for a meaningful Spearman rho

BOOT_SEED <- 2024L   # matches 03_compare_to_memento.R's seed for consistency
N_BOOT    <- 2000L
CI_LEVEL  <- 0.95

BASICS_DIR   <- "figure_reproduction/outputs/basics_scran_crossmethod/test_de"
COMPARATOR   <- Sys.getenv("DD_COMPARATOR", "matched_1000cell")
COMPARATOR_FITS <- c(
  matched_1000cell = "figure_reproduction/outputs/basics_scran_crossmethod/memento_1000subset/dream_fit_1000subset.rds",
  fullcell         = "figure_reproduction/outputs/dd_refit/dream_fit_full.rds")
stopifnot("DD_COMPARATOR must be matched_1000cell or fullcell" = COMPARATOR %in% names(COMPARATOR_FITS))
MEMENTO_FIT  <- COMPARATOR_FITS[[COMPARATOR]]
message("Memento comparator: ", COMPARATOR, " (", MEMENTO_FIT, ")")
OUT_DIR      <- file.path("figure_reproduction/outputs/basics_scran_crossmethod/crosscheck_memento", COMPARATOR)
FIG_DIR      <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

stopifnot("Memento dream/eBayes fit not found -- run script 16 (matched_1000cell) or hc_03 (fullcell) first" =
            file.exists(MEMENTO_FIT))

# ---------------------------------------------------------------------------
# Bootstrap Spearman rho with percentile CI (paired resampling of genes) --
# identical logic to 03_compare_to_memento.R's boot_spearman, reused here for
# consistency of methodology across this analysis's cross-method comparisons.
# ---------------------------------------------------------------------------
boot_spearman <- function(x, y, n_boot, ci_level) {
  stopifnot(length(x) == length(y))
  if (length(x) < MIN_PAIRS_SPEARMAN) return(list(rho = NA_real_, lo = NA_real_, hi = NA_real_))
  n <- length(x)
  # ANALYSIS_OK[warning-suppression]: suppresses cor()'s zero-variance warning
  # (a bootstrap resample can be constant -> NA rho); rank agreement is still
  # valid, and NA reps are dropped by na.rm in the quantile() below.
  rho_obs <- suppressWarnings(cor(x, y, method = "spearman"))
  reps <- vapply(seq_len(n_boot), function(i) {
    idx <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  }, numeric(1))
  a <- (1 - ci_level) / 2
  ci <- quantile(reps, c(a, 1 - a), na.rm = TRUE)
  list(rho = rho_obs, lo = unname(ci[1]), hi = unname(ci[2]))
}

message("Loading Memento dream/eBayes fit (1000-cell/individual matched subset)...")
memento_fit <- readRDS(MEMENTO_FIT)
stopifnot("expected an MArrayLM fit object" = inherits(memento_fit, "MArrayLM"))

set.seed(BOOT_SEED)

pairs <- combn(CELL_TYPES, 2, simplify = FALSE)
summary_rows  <- vector("list", length(pairs))
merged_all    <- vector("list", length(pairs))

for (i in seq_along(pairs)) {
  ct1 <- pairs[[i]][1]; ct2 <- pairs[[i]][2]
  pair_label <- paste0(ct1, "_vs_", ct2)
  message("\n--- ", pair_label, " ---")

  basics_csv <- file.path(BASICS_DIR, paste0(pair_label, "_resdisp.csv"))
  stopifnot("BASiCS resdisp CSV not found -- run 04 first" = file.exists(basics_csv))
  basics_df <- read.csv(basics_csv, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot("resdisp CSV missing expected columns" =
              all(c("GeneName", "ResDispDistance", "ResultDiffResDisp") %in% colnames(basics_df)))

  memento_contrast <- paste0("compare_", ct2, "_", ct1)
  stopifnot("Memento contrast not found in fit -- CELL_TYPES/pairing mismatch" =
              memento_contrast %in% colnames(memento_fit$coefficients))
  tt <- topTable(memento_fit, coef = memento_contrast, number = Inf, sort.by = "none")
  tt$GeneName <- rownames(tt)
  # Sign flip: Memento's raw logFC for "compare_{ct2}_{ct1}" is ct2 - ct1; negate
  # to align to BASiCS's ResDispDistance = ct1 - ct2 (see header note).
  tt$memento_logFC_1m2 <- -tt$logFC

  merged <- merge(basics_df[, c("GeneName", "ResDispDistance", "ResultDiffResDisp")],
                   tt[, c("GeneName", "memento_logFC_1m2", "adj.P.Val")],
                   by = "GeneName")
  message("  BASiCS genes=", nrow(basics_df), " Memento genes=", nrow(tt),
          " shared/tested=", nrow(merged))
  stopifnot("fewer than MIN_SHARED_GENES shared tested genes -- gene-name alignment likely broken" =
              nrow(merged) >= MIN_SHARED_GENES)

  sig_memento <- !is.na(merged$adj.P.Val) & merged$adj.P.Val <= DD_ADJP_THRESH
  sig_basics  <- merged$ResultDiffResDisp %in% c(paste0(ct1, "+"), paste0(ct2, "+"))
  sig_both    <- sig_memento & sig_basics
  n_memento <- sum(sig_memento); n_basics <- sum(sig_basics)
  n_both    <- sum(sig_both);    n_union  <- sum(sig_memento | sig_basics)

  pct_jaccard      <- 100 * n_both / n_union
  pct_of_basics    <- 100 * n_both / n_basics
  pct_of_memento   <- 100 * n_both / n_memento

  message(sprintf("  significant: Memento=%d BASiCS=%d both=%d union=%d",
                  n_memento, n_basics, n_both, n_union))
  message(sprintf("  overlap: Jaccard=%.1f%%  (both/BASiCS-sig)=%.1f%%  (both/Memento-sig)=%.1f%%",
                  pct_jaccard, pct_of_basics, pct_of_memento))

  bs_sig <- boot_spearman(merged$memento_logFC_1m2[sig_both],
                          merged$ResDispDistance[sig_both], N_BOOT, CI_LEVEL)
  # Third category (Brendan's follow-up request): correlation restricted to
  # genes Memento calls significant, REGARDLESS of the BASiCS call -- distinct
  # from sig_both (both methods agree) and all_tested (no significance filter
  # at all). Answers "for genes Memento says are DD, how well does BASiCS's
  # effect size track Memento's logFC, whether or not BASiCS also flagged it."
  bs_sig_memento <- boot_spearman(merged$memento_logFC_1m2[sig_memento],
                                  merged$ResDispDistance[sig_memento], N_BOOT, CI_LEVEL)
  bs_all <- boot_spearman(merged$memento_logFC_1m2, merged$ResDispDistance, N_BOOT, CI_LEVEL)
  message(sprintf("  logFC Spearman rho (sig-in-both, n=%d) = %.3f [%.3f, %.3f]",
                  n_both, bs_sig$rho, bs_sig$lo, bs_sig$hi))
  message(sprintf("  logFC Spearman rho (Memento-sig, n=%d) = %.3f [%.3f, %.3f]",
                  n_memento, bs_sig_memento$rho, bs_sig_memento$lo, bs_sig_memento$hi))
  message(sprintf("  logFC Spearman rho (all shared tested, n=%d) = %.3f [%.3f, %.3f]",
                  nrow(merged), bs_all$rho, bs_all$lo, bs_all$hi))

  summary_rows[[i]] <- data.frame(
    cell_type_1 = ct1, cell_type_2 = ct2,
    n_shared_tested = nrow(merged),
    n_sig_memento = n_memento, n_sig_basics = n_basics,
    n_sig_both = n_both, n_sig_union = n_union,
    pct_overlap_jaccard = round(pct_jaccard, 2),
    pct_of_basics_sig_also_memento = round(pct_of_basics, 2),
    pct_of_memento_sig_also_basics = round(pct_of_memento, 2),
    logfc_spearman_rho_sig_both = bs_sig$rho, logfc_spearman_ci_lo_sig_both = bs_sig$lo,
    logfc_spearman_ci_hi_sig_both = bs_sig$hi,
    logfc_spearman_rho_sig_memento = bs_sig_memento$rho,
    logfc_spearman_ci_lo_sig_memento = bs_sig_memento$lo,
    logfc_spearman_ci_hi_sig_memento = bs_sig_memento$hi,
    logfc_spearman_rho_all_tested = bs_all$rho, logfc_spearman_ci_lo_all_tested = bs_all$lo,
    logfc_spearman_ci_hi_all_tested = bs_all$hi,
    stringsAsFactors = FALSE
  )

  merged$cell_type_1 <- ct1
  merged$cell_type_2 <- ct2
  merged$sig_memento <- sig_memento
  merged$sig_basics  <- sig_basics
  merged$sig_both    <- sig_both
  write.csv(merged, file.path(OUT_DIR, paste0(pair_label, "_merged.csv")), row.names = FALSE)
  merged_all[[i]] <- merged
}

summary_df <- do.call(rbind, summary_rows)
stopifnot("pairwise crosscheck summary row count mismatch" = nrow(summary_df) == length(pairs))
write.csv(summary_df, file.path(OUT_DIR, "pairwise_crosscheck_summary.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Pooled summary across all 15 pairs (genes stacked per category) -- an
# overall headline number alongside the per-pair breakdown above.
# ---------------------------------------------------------------------------
merged_stacked <- do.call(rbind, merged_all)
pooled_sig_both    <- merged_stacked[merged_stacked$sig_both, ]
pooled_sig_memento <- merged_stacked[merged_stacked$sig_memento, ]
bs_pooled          <- boot_spearman(pooled_sig_both$memento_logFC_1m2, pooled_sig_both$ResDispDistance,
                                    N_BOOT, CI_LEVEL)
bs_pooled_memento  <- boot_spearman(pooled_sig_memento$memento_logFC_1m2, pooled_sig_memento$ResDispDistance,
                                    N_BOOT, CI_LEVEL)
pooled_summary <- data.frame(
  n_pairs = length(pairs),
  n_sig_both_pooled = nrow(pooled_sig_both),
  n_sig_memento_pooled = nrow(pooled_sig_memento),
  mean_pct_overlap_jaccard = round(mean(summary_df$pct_overlap_jaccard), 2),
  median_pct_overlap_jaccard = round(median(summary_df$pct_overlap_jaccard), 2),
  logfc_spearman_rho_pooled_sig_both = bs_pooled$rho,
  logfc_spearman_ci_lo_pooled_sig_both = bs_pooled$lo,
  logfc_spearman_ci_hi_pooled_sig_both = bs_pooled$hi,
  logfc_spearman_rho_pooled_sig_memento = bs_pooled_memento$rho,
  logfc_spearman_ci_lo_pooled_sig_memento = bs_pooled_memento$lo,
  logfc_spearman_ci_hi_pooled_sig_memento = bs_pooled_memento$hi,
  stringsAsFactors = FALSE
)
write.csv(pooled_summary, file.path(OUT_DIR, "pooled_crosscheck_summary.csv"), row.names = FALSE)

message("\n=== Pooled across all ", length(pairs), " pairs ===")
message(sprintf("  Mean pairwise Jaccard overlap = %.1f%% (median %.1f%%)",
                pooled_summary$mean_pct_overlap_jaccard, pooled_summary$median_pct_overlap_jaccard))
message(sprintf("  Pooled logFC Spearman rho (Memento-sig, n=%d) = %.3f [%.3f, %.3f]",
                pooled_summary$n_sig_memento_pooled, bs_pooled_memento$rho,
                bs_pooled_memento$lo, bs_pooled_memento$hi))
message(sprintf("  Pooled logFC Spearman rho (sig-in-both, n=%d) = %.3f [%.3f, %.3f]",
                pooled_summary$n_sig_both_pooled, bs_pooled$rho, bs_pooled$lo, bs_pooled$hi))

# ---------------------------------------------------------------------------
# Faceted scatter: Memento logFC vs BASiCS ResDispDistance, sig-in-both genes,
# one facet per pair.
# ---------------------------------------------------------------------------
# ANALYSIS_OK[plot-filter]: intentional -- this scatter is specifically about
# the sig-in-both population (the "shared DD genes" Brendan asked to
# correlate), not the full tested population; the analogous all-tested-genes
# rho is reported separately in the summary CSVs (logfc_spearman_rho_all_tested).
plot_df <- merged_stacked[merged_stacked$sig_both, ]
plot_df$pair_label <- paste0(plot_df$cell_type_1, "\nvs ", plot_df$cell_type_2)
p <- ggplot(plot_df, aes(memento_logFC_1m2, ResDispDistance)) +
  geom_point(alpha = 0.3, size = 0.6) +
  facet_wrap(~ pair_label, scales = "free") +
  labs(title = "Memento DD logFC vs. BASiCS epsilon ResDispDistance",
       subtitle = "Genes significant in BOTH methods per pair (cell_type_1 - cell_type_2 convention)",
       x = "Memento logFC (cell_type_1 - cell_type_2)",
       y = "BASiCS ResDispDistance (cell_type_1 - cell_type_2)") +
  theme_bw(base_size = 8)
ggsave(file.path(FIG_DIR, "pairwise_logfc_scatter.pdf"), p, width = 14, height = 10, device = cairo_pdf)

message("\nDone. Outputs written to ", OUT_DIR)
print(summary_df[, c("cell_type_1", "cell_type_2", "pct_overlap_jaccard",
                     "logfc_spearman_rho_sig_both", "logfc_spearman_rho_sig_memento")])

# ---- reportable values (R-side register_value shim; see dd_excl_topexpressed/
# scripts/06_loess_mean_dispersion_diagnostic.R's identical-pattern precedent) --
reportable <- do.call(rbind, lapply(seq_len(nrow(summary_df)), function(i) {
  row <- summary_df[i, ]
  pair_key <- paste0(row$cell_type_1, "_vs_", row$cell_type_2)
  data.frame(
    key = c(paste0("pct_overlap_jaccard_", pair_key),
            paste0("logfc_spearman_rho_sig_both_", pair_key),
            paste0("logfc_spearman_rho_sig_memento_", pair_key)),
    value = c(row$pct_overlap_jaccard, signif(row$logfc_spearman_rho_sig_both, 4),
              signif(row$logfc_spearman_rho_sig_memento, 4)),
    provenance = "figure_reproduction/scripts/basics_scran_crossmethod/13_crosscheck_memento_pairwise_dd.R",
    stringsAsFactors = FALSE)
}))
reportable <- rbind(reportable, data.frame(
  key = c("pooled_logfc_spearman_rho_sig_both", "pooled_logfc_spearman_rho_sig_memento",
          "pooled_mean_pct_overlap_jaccard"),
  value = c(signif(pooled_summary$logfc_spearman_rho_pooled_sig_both, 4),
            signif(pooled_summary$logfc_spearman_rho_pooled_sig_memento, 4),
            pooled_summary$mean_pct_overlap_jaccard),
  provenance = "figure_reproduction/scripts/basics_scran_crossmethod/13_crosscheck_memento_pairwise_dd.R",
  stringsAsFactors = FALSE))
write.csv(reportable, file.path(OUT_DIR, "reportable_values_13.csv"), row.names = FALSE)
