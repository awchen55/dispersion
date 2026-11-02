#!/usr/bin/env Rscript

stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(WGCNA)
  library(dynamicTreeCut)
})

FR_DIR   <- "figure_reproduction/outputs/coexpression_plasticity_lowfidelity"
ORIG_DIR <- "analysis/coexpression_plasticity_lowfidelity/outputs"
TOM_DIR  <- file.path(FR_DIR, "consensus_tom")
QC_DIR   <- file.path(FR_DIR, "consensus_module_qc")
MEMB_DIR <- file.path(FR_DIR, "consensus_modules_prototype")
dir.create(TOM_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

CELL_TYPES <- c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells")
INDIVIDUALS <- c("NA19099", "NA19093", "NA18522")

FLOOR <- 0.25
DEEP_SPLIT <- 2
MIN_CLUSTER_SIZE <- 10
CALIBRATION <- "single quantile"
CALIBRATION_QUANTILE <- 0.95
CONSENSUS_QUANTILE <- 0.5
NPERM <- 200
PERM_SEED <- 20260929

log_msg <- function(...) message("[08_consensus_tom_and_cohesion] ", ...)
t_start <- Sys.time()

orig_power_df <- read.csv(file.path(ORIG_DIR, "diagnostics", "first_crossing_power_summary.csv"),
                          stringsAsFactors = FALSE)
epi_power_df <- read.csv(file.path(FR_DIR, "diagnostics", "first_crossing_power_summary_0.25.csv"),
                         stringsAsFactors = FALSE)

get_power <- function(ct, ind) {
  if (ct == "Epithelial_cells") {
    prow <- epi_power_df[epi_power_df$individual == ind, ]
    if (nrow(prow) != 1L) stop(ct, "/", ind, ": expected 1 row in Epithelial power CSV, got ", nrow(prow))
    return(prow$power)
  }
  prow <- orig_power_df[orig_power_df$cell_type == ct & orig_power_df$individual == ind &
                          abs(orig_power_df$floor - FLOOR) < 1e-9, ]
  if (nrow(prow) != 1L) stop(ct, "/", ind, ": expected 1 row in power CSV for floor=", FLOOR, ", got ", nrow(prow))
  if (prow$crosses_at_all) prow$first_crossing_power else prow$global_best_power
}

mean_offdiag <- function(sub) (sum(sub) - sum(diag(sub))) / (nrow(sub) * (nrow(sub) - 1))

module_rows <- list(); summary_rows <- list()
set.seed(PERM_SEED)

for (ct in CELL_TYPES) {
  log_msg("##### ", ct, " #####")
  in_dir <- if (ct == "Epithelial_cells") file.path(FR_DIR, "cscore_modules") else file.path(ORIG_DIR, "cscore_modules")
  powers <- setNames(vapply(INDIVIDUALS, function(ind) get_power(ct, ind), numeric(1)), INDIVIDUALS)
  log_msg(ct, ": powers ", paste(names(powers), powers, sep = "=", collapse = ", "))

  coexp_list <- list()
  for (ind in INDIVIDUALS) {
    rds <- file.path(in_dir, paste0(ct, "_", ind, "_coexp_thresholded.rds"))
    if (!file.exists(rds)) stop("Missing expected file: ", rds)
    coexp_list[[ind]] <- readRDS(rds)
  }
  shared_genes <- sort(Reduce(intersect, lapply(coexp_list, rownames)))
  n_shared <- length(shared_genes)
  stopifnot("shared gene universe implausibly small (<1000)" = n_shared >= 1000)

  tom_list <- list()
  for (ind in INDIVIDUALS) {
    coexp <- coexp_list[[ind]][shared_genes, shared_genes]
    sim <- abs(coexp)
    sim[sim < FLOOR] <- 0
    diag(sim) <- 1
    adj <- WGCNA::adjacency.fromSimilarity(sim, power = unname(powers[ind]))
    tom_list[[ind]] <- WGCNA::TOMsimilarity(adj, verbose = 0)
    rm(coexp, sim, adj); gc(FALSE)
  }
  rm(coexp_list)

  opts <- WGCNA::newConsensusOptions(calibration = CALIBRATION, calibrationQuantile = CALIBRATION_QUANTILE,
                                     sampleForCalibration = FALSE, consensusQuantile = CONSENSUS_QUANTILE)
  res <- WGCNA::consensusCalculation(individualData = unname(tom_list), consensusOptions = opts, verbose = 0)
  consensus_tom <- res$consensusData$data[[1]]
  rownames(consensus_tom) <- colnames(consensus_tom) <- shared_genes
  stopifnot("consensus TOM dimensions don't match shared gene universe" =
              identical(dim(consensus_tom), c(n_shared, n_shared)))
  rm(tom_list, res); gc(FALSE)

  # --- verify recomputed modules == saved consensus membership --------------
  diss <- 1 - consensus_tom
  memb <- dynamicTreeCut::cutreeDynamic(
    dendro = stats::hclust(stats::as.dist(diss), method = "average"), distM = diss,
    deepSplit = DEEP_SPLIT, pamRespectsDendro = FALSE, minClusterSize = MIN_CLUSTER_SIZE
  )
  names(memb) <- shared_genes
  rm(diss)
  saved_csv <- file.path(MEMB_DIR, paste0(ct, "_consensus_q0.5_module_membership.csv"))
  if (!file.exists(saved_csv)) stop("Missing ", saved_csv)
  saved <- read.csv(saved_csv, stringsAsFactors = FALSE)
  if (!setequal(saved$gene, shared_genes)) stop(ct, ": saved membership gene set differs from recomputed shared universe")
  saved_memb <- setNames(saved$module, saved$gene)[shared_genes]
  n_diff <- sum(saved_memb != memb)
  if (n_diff > 0) {
    stop(ct, ": recomputed consensus modules differ from the saved membership CSV for ", n_diff, " of ",
         n_shared, " genes -- cohesion would describe a different module set; investigate before proceeding.")
  }
  log_msg(ct, ": recomputed modules match saved membership exactly (", n_shared, " genes, ",
          length(unique(memb[memb != 0])), " modules).")

  tom_path <- file.path(TOM_DIR, paste0(ct, "_consensus_tom.rds"))
  saveRDS(consensus_tom, tom_path)
  log_msg(ct, ": wrote ", tom_path, " (", round(file.size(tom_path) / 1e6), " MB)")

  # --- cohesion ---------------------------------------------------------------
  mods <- sort(unique(memb[memb != 0]))
  assigned <- which(memb != 0)
  for (m in mods) {
    idx <- which(memb == m)
    other <- setdiff(assigned, idx)
    n <- length(idx)
    within <- mean_offdiag(consensus_tom[idx, idx, drop = FALSE])
    between <- if (length(other) > 0) mean(consensus_tom[idx, other, drop = FALSE]) else NA_real_
    null <- vapply(seq_len(NPERM), function(i) {
      r <- sample.int(n_shared, n)
      mean_offdiag(consensus_tom[r, r, drop = FALSE])
    }, numeric(1))
    module_rows[[paste(ct, m)]] <- data.frame(
      cell_type = ct, module = m, n_genes = n,
      mean_within = within, mean_between = between,
      within_minus_between = within - between,
      null_mean = mean(null), null_p95 = unname(quantile(null, 0.95)),
      perm_p = (1 + sum(null >= within)) / (1 + NPERM)
    )
  }
  ct_df <- do.call(rbind, module_rows[startsWith(names(module_rows), paste0(ct, " "))])
  summary_rows[[ct]] <- data.frame(
    cell_type = ct, n_modules = nrow(ct_df),
    median_within = median(ct_df$mean_within), median_between = median(ct_df$mean_between, na.rm = TRUE),
    n_within_gt_between = sum(ct_df$mean_within > ct_df$mean_between, na.rm = TRUE),
    n_perm_p_lt_0.05 = sum(ct_df$perm_p < 0.05)
  )
  log_msg(ct, ": ", nrow(ct_df), " modules; within > between in ", summary_rows[[ct]]$n_within_gt_between,
          "; perm p<0.05 in ", summary_rows[[ct]]$n_perm_p_lt_0.05, ".")
  rm(consensus_tom); gc(FALSE)
}

module_df <- do.call(rbind, module_rows); rownames(module_df) <- NULL
summary_df <- do.call(rbind, summary_rows); rownames(summary_df) <- NULL
stopifnot("NA within-module cohesion" = !anyNA(module_df$mean_within))
write.csv(module_df, file.path(QC_DIR, "consensus_module_cohesion.csv"), row.names = FALSE)
write.csv(summary_df, file.path(QC_DIR, "consensus_module_cohesion_summary.csv"), row.names = FALSE)
log_msg("DONE in ", round(difftime(Sys.time(), t_start, units = "mins"), 2), " min.")
