#!/usr/bin/env Rscript

stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(WGCNA)
  library(dynamicTreeCut)
})

FLOOR_STR <- Sys.getenv("CC_MODULE_FLOOR")
if (!nzchar(FLOOR_STR)) stop("CC_MODULE_FLOOR environment variable must be set -- aborting (no default).")
IS_NOFLOOR <- identical(FLOOR_STR, "nofloor")
FLOOR <- if (IS_NOFLOOR) NA_real_ else as.numeric(FLOOR_STR)
if (!IS_NOFLOOR) stopifnot("CC_MODULE_FLOOR must be numeric or the literal string 'nofloor'" = !is.na(FLOOR))

IN_MODULES_DIR <- "figure_reproduction/outputs/coexpression_plasticity_lowfidelity/cscore_modules"
DIAG_DIR <- "figure_reproduction/outputs/coexpression_plasticity_lowfidelity/diagnostics"
OUT_MODULES_DIR <- file.path(
  "figure_reproduction/outputs/coexpression_plasticity_lowfidelity",
  if (IS_NOFLOOR) "cscore_modules_nofloor_fcpower" else sprintf("cscore_modules_floor%s_fcpower", FLOOR_STR)
)
dir.create(DIAG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_MODULES_DIR, recursive = TRUE, showWarnings = FALSE)

CELL_TYPE <- "Epithelial_cells"
INDIVIDUALS <- c("NA19099", "NA19093", "NA18522")
POWER_VECTOR <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))  # matches source scripts 05-14
RSQUARED_CUT <- 0.85
DEEP_SPLIT <- 2
MIN_CLUSTER_SIZE <- 10
KNOWN_WEAK_MEAN_K_THRESHOLD <- 1  # matches source scripts' flagging convention

log_msg <- function(...) message("[04_floor_sweep_first_crossing_epithelial] ", ...)

summary_rows <- list()

for (ind in INDIVIDUALS) {
  sample_label <- paste0(CELL_TYPE, "_", ind)
  log_msg("=== ", sample_label, " (floor=", FLOOR_STR, ") ===")

  coexp_rds <- file.path(IN_MODULES_DIR, paste0(sample_label, "_coexp_thresholded.rds"))
  if (!file.exists(coexp_rds)) stop("Missing expected file: ", coexp_rds, " -- run script 03 first.")
  coexp <- readRDS(coexp_rds)
  genes_selected <- rownames(coexp)
  n_genes <- length(genes_selected)

  sim <- abs(coexp)
  if (!IS_NOFLOOR) {
    # ANALYSIS_OK[feature-filter]: applies the same effect-size floor
    # mechanism as source scripts 07/09/11/14 (Brendan's method); retained-
    # pair count logged below.
    sim[sim < FLOOR] <- 0
  }
  diag(sim) <- 1
  n_retained <- sum(sim[upper.tri(sim)] != 0)

  sft <- WGCNA::pickSoftThreshold.fromSimilarity(
    similarity = sim, RsquaredCut = RSQUARED_CUT, powerVector = POWER_VECTOR, verbose = 0
  )
  fit_tab <- sft$fitIndices
  crosses <- fit_tab$SFT.R.sq >= RSQUARED_CUT
  best_idx <- which.max(fit_tab$SFT.R.sq)
  if (any(crosses)) {
    first_idx <- which(crosses)[1]
    power <- fit_tab$Power[first_idx]
    power_source <- "first-crossing"
    chosen_r2 <- fit_tab$SFT.R.sq[first_idx]
    chosen_mean_k <- fit_tab$mean.k.[first_idx]
  } else {
    power <- fit_tab$Power[best_idx]
    power_source <- "global-best (never crosses 0.85)"
    chosen_r2 <- fit_tab$SFT.R.sq[best_idx]
    chosen_mean_k <- fit_tab$mean.k.[best_idx]
  }
  log_msg("  power selection: ", power, " (", power_source, "), R2=", round(chosen_r2, 3),
          ", mean.k=", round(chosen_mean_k, 3))

  # --- immediately rerun WGCNA module detection at this power (matches
  #     source scripts' identical recipe) ---
  adj <- WGCNA::adjacency.fromSimilarity(sim, power = power)
  mean_k <- mean(rowSums(adj) - diag(adj))
  if (mean_k < KNOWN_WEAK_MEAN_K_THRESHOLD) {
    log_msg("  FLAGGED: mean connectivity k=", round(mean_k, 3),
            " is low (<", KNOWN_WEAK_MEAN_K_THRESHOLD, ") even at the chosen power -- ",
            "treat this sample's module calls as marginal.")
  }

  tom <- WGCNA::TOMsimilarity(adj, verbose = 0)
  dissTOM <- 1 - tom
  rownames(dissTOM) <- colnames(dissTOM) <- genes_selected
  hclust_dist <- stats::hclust(stats::as.dist(dissTOM), method = "average")
  memb <- dynamicTreeCut::cutreeDynamic(
    dendro = hclust_dist, distM = dissTOM, deepSplit = DEEP_SPLIT,
    pamRespectsDendro = FALSE, minClusterSize = MIN_CLUSTER_SIZE
  )
  names(memb) <- genes_selected
  stopifnot("module membership vector length mismatch" = length(memb) == n_genes)

  n_modules <- length(unique(memb[memb != 0]))
  n_unassigned <- sum(memb == 0)
  log_msg("  ", n_modules, " module(s) (excl. 0); ", n_unassigned, " of ", n_genes,
          " genes (", round(100 * n_unassigned / n_genes, 2), "%) unassigned (module 0); ",
          "mean k=", round(mean_k, 2))

  module_df <- data.frame(
    gene = genes_selected, module = as.integer(memb), in_module = memb != 0,
    cell_type = CELL_TYPE, individual = ind, floor = FLOOR_STR, power = power,
    power_source = power_source, mean_connectivity = mean_k,
    row.names = NULL
  )
  out_csv <- file.path(OUT_MODULES_DIR, paste0(sample_label, "_module_membership.csv"))
  utils::write.csv(module_df, out_csv, row.names = FALSE)
  log_msg("  Wrote: ", out_csv)

  summary_rows[[sample_label]] <- data.frame(
    cell_type = CELL_TYPE, individual = ind, floor = FLOOR_STR,
    power = power, power_source = power_source,
    crosses_at_all = any(crosses), chosen_r2 = chosen_r2, chosen_mean_k = chosen_mean_k,
    n_genes = n_genes, n_pairs_retained = n_retained, n_modules = n_modules,
    n_unassigned = n_unassigned, pct_unassigned = 100 * n_unassigned / n_genes
  )

  rm(coexp, sim, sft, fit_tab, adj, tom, dissTOM, hclust_dist, memb, module_df)
  gc(FALSE)
}

summary_df <- do.call(rbind, summary_rows); rownames(summary_df) <- NULL
stopifnot("expected 3 (one per individual) rows" = nrow(summary_df) == length(INDIVIDUALS))
summary_csv <- file.path(DIAG_DIR, sprintf("first_crossing_power_summary_%s.csv", FLOOR_STR))
write.csv(summary_df, summary_csv, row.names = FALSE)
log_msg("Wrote: ", summary_csv)

log_msg("DONE for floor=", FLOOR_STR, ".")
