#!/usr/bin/env Rscript

stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(WGCNA)
  library(dynamicTreeCut)
})

BASE_DIR       <- "figure_reproduction/outputs/coexpression_plasticity_lowfidelity"
MODULES_DIR    <- file.path(BASE_DIR, "cscore_modules")
INDIV_MOD_DIR  <- file.path(BASE_DIR, "cscore_modules_floor0.25_fcpower")
DIAG_DIR       <- file.path(BASE_DIR, "diagnostics")
POWER_CSV      <- file.path(DIAG_DIR, "first_crossing_power_summary_0.25.csv")
OUT_DIR        <- file.path(BASE_DIR, "consensus_modules_prototype")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CELL_TYPE <- "Epithelial_cells"
INDIVIDUALS <- c("NA19099", "NA19093", "NA18522")

FLOOR <- 0.25  # finalized go-forward parameter, .living/decisions.md 2026-09-17
DEEP_SPLIT <- 2
MIN_CLUSTER_SIZE <- 10

CALIBRATION <- "single quantile"
CALIBRATION_QUANTILE <- 0.95
CONSENSUS_QUANTILE <- 0.5

log_msg <- function(...) message("[06_consensus_wgcna_epithelial] ", ...)
t_start <- Sys.time()

power_df <- read.csv(POWER_CSV, stringsAsFactors = FALSE)
stopifnot(
  "power CSV missing expected columns" =
    all(c("cell_type", "individual", "power") %in% colnames(power_df)),
  "expected 3 rows (one per individual)" = nrow(power_df) == length(INDIVIDUALS)
)

get_power <- function(ind) {
  prow <- power_df[power_df$individual == ind, ]
  if (nrow(prow) != 1L) stop(ind, ": expected exactly 1 row in ", POWER_CSV, ", got ", nrow(prow))
  prow$power
}

powers <- setNames(vapply(INDIVIDUALS, get_power, numeric(1)), INDIVIDUALS)
log_msg(CELL_TYPE, ": per-individual power (floor=", FLOOR, "): ",
        paste(names(powers), powers, sep = "=", collapse = ", "))

coexp_list <- list()
for (ind in INDIVIDUALS) {
  rds <- file.path(MODULES_DIR, paste0(CELL_TYPE, "_", ind, "_coexp_thresholded.rds"))
  if (!file.exists(rds)) stop("Missing expected file: ", rds)
  coexp_list[[ind]] <- readRDS(rds)
}

n_each <- vapply(coexp_list, nrow, integer(1))
stopifnot(
  "individuals have different total gene-universe sizes -- investigate before computing % retained" =
    length(unique(n_each)) == 1
)
shared_genes <- sort(Reduce(intersect, lapply(coexp_list, rownames)))
n_shared <- length(shared_genes)
log_msg(CELL_TYPE, ": shared gene universe (intersection of all 3): ", n_shared, "/",
        n_each[1], " (", round(100 * n_shared / n_each[1], 1), "%).")
stopifnot(
  "shared gene universe implausibly small (<1000) -- investigate before proceeding" =
    n_shared >= 1000
)

tom_list <- list()
mean_k_list <- c()
for (ind in INDIVIDUALS) {
  coexp <- coexp_list[[ind]][shared_genes, shared_genes]
  stopifnot(
    "coexp subset not square/aligned to shared_genes" =
      identical(rownames(coexp), shared_genes) && identical(colnames(coexp), shared_genes)
  )
  sim <- abs(coexp)
  sim[sim < FLOOR] <- 0
  diag(sim) <- 1

  power <- unname(powers[ind])
  adj <- WGCNA::adjacency.fromSimilarity(sim, power = power)
  mean_k <- mean(rowSums(adj) - diag(adj))
  mean_k_list[ind] <- mean_k
  log_msg(CELL_TYPE, "/", ind, ": power=", power, ", mean connectivity k=", round(mean_k, 2),
          " (adjacency on ", n_shared, " shared genes).")

  tom_list[[ind]] <- WGCNA::TOMsimilarity(adj, verbose = 0)
  rm(coexp, sim, adj)
  gc(FALSE)
}
log_msg(CELL_TYPE, ": per-individual mean connectivity spread: ",
        paste(names(mean_k_list), round(mean_k_list, 2), sep = "=", collapse = ", "))

opts <- WGCNA::newConsensusOptions(
  calibration = CALIBRATION,
  calibrationQuantile = CALIBRATION_QUANTILE,
  sampleForCalibration = FALSE,
  consensusQuantile = CONSENSUS_QUANTILE
)
res <- WGCNA::consensusCalculation(individualData = unname(tom_list), consensusOptions = opts, verbose = 0)
consensus_tom <- res$consensusData$data[[1]]
rownames(consensus_tom) <- colnames(consensus_tom) <- shared_genes
stopifnot(
  "consensus TOM dimensions don't match shared gene universe" =
    identical(dim(consensus_tom), c(n_shared, n_shared))
)

diss_consensus <- 1 - consensus_tom
hclust_dist <- stats::hclust(stats::as.dist(diss_consensus), method = "average")
memb <- dynamicTreeCut::cutreeDynamic(
  dendro = hclust_dist, distM = diss_consensus, deepSplit = DEEP_SPLIT,
  pamRespectsDendro = FALSE, minClusterSize = MIN_CLUSTER_SIZE
)
names(memb) <- shared_genes
stopifnot("consensus module membership vector length mismatch" = length(memb) == n_shared)

n_modules <- length(unique(memb[memb != 0]))
n_unassigned <- sum(memb == 0)
tab <- sort(table(memb), decreasing = TRUE)
nonzero_tab <- tab[names(tab) != "0"]
largest_module_size <- if (length(nonzero_tab) > 0) unname(nonzero_tab[1]) else 0L
log_msg(CELL_TYPE, ": ", n_modules, " consensus module(s) (excl. 0); ", n_unassigned, " of ",
        n_shared, " genes (", round(100 * n_unassigned / n_shared, 2), "%) unassigned; largest module = ",
        largest_module_size, " genes (", round(100 * largest_module_size / n_shared, 1), "%).")

module_df <- data.frame(
  gene = shared_genes, module = as.integer(memb), in_module = memb != 0,
  cell_type = CELL_TYPE, consensus_quantile = CONSENSUS_QUANTILE, consensus_quantile_name = "majority_2of3",
  calibration = CALIBRATION, calibration_quantile = CALIBRATION_QUANTILE, floor = FLOOR,
  row.names = NULL
)
out_csv <- file.path(OUT_DIR, paste0(CELL_TYPE, "_consensus_q0.5_module_membership.csv"))
utils::write.csv(module_df, out_csv, row.names = FALSE)
log_msg(CELL_TYPE, ": wrote ", out_csv)

# Sanity check: agreement with each individual's own independent per-sample call.
indiv_in_module <- list()
for (ind in INDIVIDUALS) {
  csv <- file.path(INDIV_MOD_DIR, paste0(CELL_TYPE, "_", ind, "_module_membership.csv"))
  if (!file.exists(csv)) stop("Missing expected file: ", csv)
  df <- read.csv(csv, stringsAsFactors = FALSE)
  m <- setNames(df$in_module, df$gene)
  indiv_in_module[[ind]] <- unname(m[shared_genes])
  if (any(is.na(indiv_in_module[[ind]]))) {
    stop(CELL_TYPE, "/", ind, ": missing in_module calls for some shared genes")
  }
}
n_all3_agree_indiv <- sum(Reduce(`&`, indiv_in_module))
log_msg(CELL_TYPE, ": reference baseline -- ", n_all3_agree_indiv, " of ", n_shared,
        " shared genes are in-module in ALL 3 individuals' OWN independent per-sample calls.")

log_msg("DONE in ", round(difftime(Sys.time(), t_start, units = "mins"), 2), " min.")
