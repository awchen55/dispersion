stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(BASiCS)
  library(SingleCellExperiment)
  library(BiocParallel)
  library(coda)
})

MCMC_N     <- 20000L
MCMC_THIN  <- 10L
MCMC_BURN  <- 10000L
MCMC_SEED  <- 42L

NSUBSETS <- as.integer(Sys.getenv("BASICS_NSUBSETS", unset = "16"))
MIN_GENES_PER_SUBSET <- 200L

MIN_ESS              <- 100
GEWEKE_ABS_CUTOFF    <- 1.96
MIN_FRAC_ESS_PASS    <- 0.90
MIN_FRAC_GEWEKE_PASS <- 0.90
MIN_BATCHES          <- 2L

CELL_TYPE <- "Epithelial_cells_1000subset"  # matches original's "{ct}_1000subset" convention

REPO_DIR <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
BASE_DIR <- file.path(REPO_DIR, "figure_reproduction/outputs/basics_scran_crossmethod")
SCE_RDS   <- file.path(BASE_DIR, "sce", "Epithelial_cells_1000subset_sce_raw.rds")
OUT_CHAIN <- file.path(BASE_DIR, "chains")
OUT_EPS   <- file.path(BASE_DIR, "epsilon")
OUT_DIAG  <- file.path(BASE_DIR, "diagnostics")
for (d in c(OUT_CHAIN, OUT_EPS, OUT_DIAG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

stopifnot("SCE not found -- run 03_build_epithelial_1000subset_sce.R first" = file.exists(SCE_RDS))

sce <- readRDS(SCE_RDS)
stopifnot(
  inherits(sce, "SingleCellExperiment"),
  "counts" %in% assayNames(sce),
  "BatchInfo" %in% colnames(colData(sce))
)
n_batches <- length(unique(colData(sce)$BatchInfo))
stopifnot("no-spike-in mode requires >= MIN_BATCHES batches" = n_batches >= MIN_BATCHES)
message(CELL_TYPE, ": ", nrow(sce), " genes x ", ncol(sce), " cells, ",
        n_batches, " batches (", paste(sort(unique(colData(sce)$BatchInfo)), collapse = ", "), ")")

message("Running BASiCS_DivideAndConquer (heavy): N=", MCMC_N, " Thin=", MCMC_THIN,
        " Burn=", MCMC_BURN, " NSubsets=", NSUBSETS)
set.seed(MCMC_SEED)
bp <- BiocParallel::SnowParam(workers = NSUBSETS, RNGseed = MCMC_SEED)
chain_list <- BASiCS_DivideAndConquer(
  Data          = sce,
  NSubsets      = NSUBSETS,
  SubsetBy      = "gene",
  WithSpikes    = FALSE,
  Regression    = TRUE,
  BPPARAM       = bp,
  N             = MCMC_N,
  Thin          = MCMC_THIN,
  Burn          = MCMC_BURN,
  PrintProgress = FALSE
)
stopifnot(is.list(chain_list) && length(chain_list) == NSUBSETS)
genes_per_subset <- vapply(chain_list, function(ch) ncol(ch@parameters$epsilon), integer(1))
message("Genes per subset: min=", min(genes_per_subset), " median=", median(genes_per_subset),
        " max=", max(genes_per_subset), " (total=", sum(genes_per_subset), ")")
stopifnot(min(genes_per_subset) >= MIN_GENES_PER_SUBSET)
saveRDS(chain_list, file.path(OUT_CHAIN, paste0(CELL_TYPE, "_chain.rds")))

ess_list <- vector("list", NSUBSETS)
geweke_z_list <- vector("list", NSUBSETS)
eps_df_list <- vector("list", NSUBSETS)
n_finite_total <- 0L
n_genes_total <- 0L

for (i in seq_len(NSUBSETS)) {
  chain_i <- chain_list[[i]]
  stopifnot("epsilon" %in% names(chain_i@parameters))
  eps_chain_i <- chain_i@parameters$epsilon
  stopifnot(nrow(eps_chain_i) > 0 && ncol(eps_chain_i) > 0)
  finite_gene_i <- apply(eps_chain_i, 2, function(col) all(is.finite(col)))
  n_finite_total <- n_finite_total + sum(finite_gene_i)
  n_genes_total <- n_genes_total + ncol(eps_chain_i)
  eps_finite_i <- eps_chain_i[, finite_gene_i, drop = FALSE]
  stopifnot(ncol(eps_finite_i) > 0)
  ess_i <- BASiCS_EffectiveSize(chain_i, Parameter = "epsilon")
  ess_list[[i]] <- ess_i[is.finite(ess_i)]
  geweke_z_i <- coda::geweke.diag(coda::mcmc(eps_finite_i))$z
  geweke_z_list[[i]] <- geweke_z_i[is.finite(geweke_z_i)]
  summ_i <- Summary(chain_i)
  eps_summary_i <- displaySummaryBASiCS(summ_i, Parameter = "epsilon")
  eps_df_i <- as.data.frame(eps_summary_i)
  eps_df_i$gene <- rownames(eps_summary_i)
  eps_df_list[[i]] <- eps_df_i
}

message("epsilon genes with finite chain (pooled across ", NSUBSETS, " subsets): ",
        n_finite_total, " / ", n_genes_total)

ess <- unlist(ess_list, use.names = FALSE)
geweke_z <- unlist(geweke_z_list, use.names = FALSE)
frac_ess_pass <- mean(ess >= MIN_ESS)
frac_geweke_pass <- mean(abs(geweke_z) < GEWEKE_ABS_CUTOFF)

diag_lines <- c(
  paste0("cell_type: ", CELL_TYPE),
  paste0("method: BASiCS_DivideAndConquer (gene-wise), NSubsets=", NSUBSETS),
  paste0("mcmc: N=", MCMC_N, " Thin=", MCMC_THIN, " Burn=", MCMC_BURN, " seed=", MCMC_SEED),
  paste0("ESS: min=", round(min(ess), 1), " median=", round(median(ess), 1),
         " frac>=", MIN_ESS, ": ", round(frac_ess_pass, 3)),
  paste0("Geweke |z|: median=", round(median(abs(geweke_z)), 2),
         " frac<", GEWEKE_ABS_CUTOFF, ": ", round(frac_geweke_pass, 3))
)
writeLines(diag_lines, file.path(OUT_DIAG, paste0(CELL_TYPE, "_convergence.txt")))
message(paste(diag_lines, collapse = "\n"))

if (frac_ess_pass < MIN_FRAC_ESS_PASS || frac_geweke_pass < MIN_FRAC_GEWEKE_PASS) {
  stop(CELL_TYPE, ": MCMC NOT CONVERGED (ESS pass=", round(frac_ess_pass, 3),
       ", Geweke pass=", round(frac_geweke_pass, 3), "). Chains saved but epsilon NOT written.")
}
message("Convergence gate PASSED for ", CELL_TYPE)

eps_df <- do.call(rbind, eps_df_list)
rownames(eps_df) <- NULL
eps_df <- eps_df[, c("gene", setdiff(colnames(eps_df), "gene"))]
stopifnot(nrow(eps_df) > 0, !any(duplicated(eps_df$gene)))

out_csv <- file.path(OUT_EPS, paste0(CELL_TYPE, "_epsilon.csv"))
write.csv(eps_df, out_csv, row.names = FALSE)
message("Saved epsilon table: ", out_csv)
