stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(BASiCS)
})

# ---------------------------------------------------------------------------
# NAMED CONSTANTS
# ---------------------------------------------------------------------------
CELL_TYPES <- c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                 "Epithelial_cells", "Activated_fibroblasts",
                 "Foregut_endodermal_cells")
N_CT <- length(CELL_TYPES)

MIN_SHARED_GENES <- 10000L  # canary: cell types share ~11,600-11,700 genes; a
                             # count far below that would mean gene-name
                             # alignment (e.g. a symbol format mismatch) broke

# BASiCS_TestDE's own package defaults, made explicit as named constants per
# repo convention (no unlabeled magic numbers) rather than left as implicit
# argument defaults.
EPSILON_M        <- log2(1.5)
EPSILON_D        <- log2(1.5)
EPSILON_R        <- log2(1.5) / log2(exp(1))
PROB_THRESHOLD_M <- 2 / 3
PROB_THRESHOLD_D <- 2 / 3
PROB_THRESHOLD_R <- 2 / 3
EFDR_M           <- 0.05
EFDR_D           <- 0.05
EFDR_R           <- 0.05
MIN_ESS_TESTDE   <- 100  # matches 02_run_basics_mcmc.R's own MCMC convergence-gate MIN_ESS

# Port change (figure_reproduction): Epithelial_cells' chain is the NEW
# sub-cluster-24-excluded dBASiCS chain (script 04 here); the other 5 cell
# types are unchanged and read read-only from the original analysis.
CHAIN_DIR_ORIGINAL   <- "analysis/basics_dispersion_reanalysis/outputs/chains"
CHAIN_DIR_EPITHELIAL <- "figure_reproduction/outputs/basics_scran_crossmethod/chains"
OUT_DIR   <- "figure_reproduction/outputs/basics_scran_crossmethod/test_de"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Stitch a dBASiCS gene-wise divide-and-conquer subset list into one synthetic
# BASiCS_Chain, combining only the gene-indexed parameters. See TECHNICAL NOTE
# 1 above.
# ---------------------------------------------------------------------------
combine_dbasics_subsets <- function(chain_list) {
  gene_names <- unlist(lapply(chain_list, function(ch) colnames(ch@parameters$mu)),
                        use.names = FALSE)
  stopifnot(
    "expected a list of BASiCS_Chain objects" =
      is.list(chain_list) && all(vapply(chain_list, inherits, logical(1), "BASiCS_Chain")),
    "gene names missing from a subset's mu matrix" = !is.null(gene_names) && length(gene_names) > 0,
    "duplicate genes across dBASiCS gene-wise subsets -- subsets should be disjoint" =
      !anyDuplicated(gene_names)
  )
  combined <- chain_list[[1]]
  combined@parameters$mu      <- do.call(cbind, lapply(chain_list, function(ch) ch@parameters$mu))
  combined@parameters$delta   <- do.call(cbind, lapply(chain_list, function(ch) ch@parameters$delta))
  combined@parameters$epsilon <- do.call(cbind, lapply(chain_list, function(ch) ch@parameters$epsilon))
  combined
}

# ---------------------------------------------------------------------------
# Restrict two combined chains to their shared gene set, in identical column
# order, for the gene-indexed parameters BASiCS_TestDE reads. See TECHNICAL
# NOTE 2 above.
# ---------------------------------------------------------------------------
restrict_to_shared_genes <- function(chain_a, chain_b) {
  shared <- intersect(colnames(chain_a@parameters$mu), colnames(chain_b@parameters$mu))
  stopifnot("fewer than MIN_SHARED_GENES shared genes between cell types -- gene-name alignment likely broken" =
              length(shared) >= MIN_SHARED_GENES)
  for (p in c("mu", "delta", "epsilon")) {
    chain_a@parameters[[p]] <- chain_a@parameters[[p]][, shared, drop = FALSE]
    chain_b@parameters[[p]] <- chain_b@parameters[[p]][, shared, drop = FALSE]
  }
  list(chain_a = chain_a, chain_b = chain_b, n_shared = length(shared))
}

# ---------------------------------------------------------------------------
# LOAD + STITCH all 6 cell types once (reused across their 5 pairings each)
# ---------------------------------------------------------------------------
message("=== Loading + stitching dBASiCS chains for all ", N_CT, " cell types ===")
combined_chains <- setNames(lapply(CELL_TYPES, function(ct) {
  chain_dir <- if (ct == "Epithelial_cells") CHAIN_DIR_EPITHELIAL else CHAIN_DIR_ORIGINAL
  f <- file.path(chain_dir, paste0(ct, "_1000subset_chain.rds"))
  stopifnot("chain file not found -- run the dBASiCS step (04) / original analysis first" = file.exists(f))
  chain_list <- readRDS(f)
  combined <- combine_dbasics_subsets(chain_list)
  message("  ", ct, ": ", length(chain_list), " subsets stitched -> ",
          ncol(combined@parameters$mu), " genes")
  combined
}), CELL_TYPES)

# ---------------------------------------------------------------------------
# PAIRWISE BASiCS_TestDE, Parameter = "ResDisp" (epsilon)
# ---------------------------------------------------------------------------
pairs <- combn(CELL_TYPES, 2, simplify = FALSE)
message("\n=== Running BASiCS_TestDE for all ", length(pairs), " pairwise cell-type comparisons ===")

summary_rows <- vector("list", length(pairs))
for (i in seq_along(pairs)) {
  ct1 <- pairs[[i]][1]; ct2 <- pairs[[i]][2]
  pair_label <- paste0(ct1, "_vs_", ct2)
  message("\n--- ", pair_label, " ---")

  restricted <- restrict_to_shared_genes(combined_chains[[ct1]], combined_chains[[ct2]])
  message("  shared genes: ", restricted$n_shared)

  test <- BASiCS_TestDE(
    Chain1 = restricted$chain_a, Chain2 = restricted$chain_b,
    GroupLabel1 = ct1, GroupLabel2 = ct2,
    EpsilonM = EPSILON_M, EpsilonD = EPSILON_D, EpsilonR = EPSILON_R,
    ProbThresholdM = PROB_THRESHOLD_M, ProbThresholdD = PROB_THRESHOLD_D,
    ProbThresholdR = PROB_THRESHOLD_R,
    EFDR_M = EFDR_M, EFDR_D = EFDR_D, EFDR_R = EFDR_R,
    MinESS = MIN_ESS_TESTDE, Offset = TRUE, Plot = FALSE
  )

  resdisp_df <- as.data.frame(test, Parameter = "ResDisp", Filter = FALSE)
  stopifnot(
    "ResDisp table missing expected result column" = "ResultDiffResDisp" %in% colnames(resdisp_df),
    "ResDisp table row count does not match shared-gene count" = nrow(resdisp_df) == restricted$n_shared
  )
  resdisp_df$cell_type_1 <- ct1
  resdisp_df$cell_type_2 <- ct2
  # ANALYSIS_OK[reorder]: column reorder only (move cell_type_1/2 to the
  # front); setdiff retains every other column, so no rows or columns are
  # dropped -- same pattern as 02_run_basics_mcmc.R's eps_df reorder.
  resdisp_df <- resdisp_df[, c("cell_type_1", "cell_type_2", setdiff(colnames(resdisp_df), c("cell_type_1", "cell_type_2")))]

  write.csv(resdisp_df, file.path(OUT_DIR, paste0(pair_label, "_resdisp.csv")), row.names = FALSE)
  saveRDS(test, file.path(OUT_DIR, paste0(pair_label, "_testde.rds")))

  # Tabulate dynamically rather than hardcoding category labels: possible
  # values include "NoDiff", "<GroupLabel>+" (direction of higher epsilon),
  # and "ExcludedLowESS" (gene failed the MinESS gate for this pair's shared
  # gene subset -- confirmed via a real-data smoke test, 2026-09-08).
  result_tab <- table(resdisp_df$ResultDiffResDisp, useNA = "ifany")
  n_higher_1 <- sum(resdisp_df$ResultDiffResDisp == paste0(ct1, "+"), na.rm = TRUE)
  n_higher_2 <- sum(resdisp_df$ResultDiffResDisp == paste0(ct2, "+"), na.rm = TRUE)
  n_nodiff   <- sum(resdisp_df$ResultDiffResDisp == "NoDiff", na.rm = TRUE)
  n_sig      <- n_higher_1 + n_higher_2
  n_excluded <- restricted$n_shared - n_sig - n_nodiff

  message("  ResultDiffResDisp: ", paste(names(result_tab), result_tab, sep = "=", collapse = ", "))
  message("  epsilon DE: ", n_sig, " / ", restricted$n_shared, " genes significant (",
          n_higher_1, " higher in ", ct1, ", ", n_higher_2, " higher in ", ct2, ")")

  summary_rows[[i]] <- data.frame(
    cell_type_1 = ct1, cell_type_2 = ct2,
    n_shared_genes = restricted$n_shared,
    offset_estimate = test@Offset,
    n_higher_in_1 = n_higher_1, n_higher_in_2 = n_higher_2,
    n_nodiff = n_nodiff, n_excluded = n_excluded,
    n_significant_resdisp = n_sig,
    pct_significant_resdisp = round(100 * n_sig / restricted$n_shared, 2),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
stopifnot("pairwise summary row count mismatch" = nrow(summary_df) == length(pairs))
write.csv(summary_df, file.path(OUT_DIR, "pairwise_resdisp_summary.csv"), row.names = FALSE)

message("\nDone. ", length(pairs), " pairwise epsilon differential tests written to ", OUT_DIR)
print(summary_df)

# ---- reportable values (R-side register_value shim; see dd_excl_topexpressed/
# scripts/06_loess_mean_dispersion_diagnostic.R's identical-pattern precedent --
# register_value.py infers namespace from the caller's Python file path, so an
# R script can't call it directly) ------------------------------------------
reportable <- do.call(rbind, lapply(seq_len(nrow(summary_df)), function(i) {
  row <- summary_df[i, ]
  pair_key <- paste0(row$cell_type_1, "_vs_", row$cell_type_2)
  data.frame(
    key = c(paste0("n_significant_resdisp_", pair_key),
            paste0("pct_significant_resdisp_", pair_key)),
    value = c(row$n_significant_resdisp, row$pct_significant_resdisp),
    provenance = "figure_reproduction/scripts/basics_scran_crossmethod/12_test_de_pairwise_celltypes.R",
    stringsAsFactors = FALSE)
}))
write.csv(reportable, file.path(dirname(OUT_DIR), "reportable_values_12.csv"), row.names = FALSE)
