# SETUP ########################################################################

library(tidyverse)
library(rhdf5)
library(ashr)
library(mashr)
library(flashier)
library(udr)

options(future.globals.maxSize = 500e+09)

set.seed(1)
################################################################################




# READ IN PREMASH FORMATTED OBJECT #############################################

# Fully run-versioned as of 2026-08-11 (Brendan: "fully version both scripts
# so they can be run in parallel"): reads from RUN_ID's own isolated
# reformatted_dqtls/{RUN_ID}/ready_for_mash/ (written by the now also
# RUN_ID-isolated 06a_run_fastqtl2mash.sbatch) instead of the old shared
# top-level reformatted_dqtls/ready_for_mash/ -- required, no fallback,
# matching this pipeline's established DQTL_RUN_ID convention (e.g. 04/05).
run_id <- Sys.getenv("DQTL_RUN_ID", unset = NA_character_)
if (is.na(run_id) || run_id == "") {
  stop("Set DQTL_RUN_ID explicitly to the RUN_ID 06a_run_fastqtl2mash.sbatch just processed.")
}

base_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data/analysis/dqtl_mapping"
input_dir <- file.path(base_dir, "outputs/reformatted_dqtls", run_id, "ready_for_mash")
output_dir <- file.path(base_dir, "outputs/mash")

input_mash_path <- file.path(input_dir, "all_celltypes.mash.rds")
if (!file.exists(input_mash_path)) {
  stop(sprintf("Missing %s. Run 06a_run_fastqtl2mash.sbatch (DQTL_RUN_ID=%s) first.",
               input_mash_path, run_id))
}
input_mash <- readRDS(input_mash_path)
################################################################################




# ACCOUNT FOR CORRELATIONS AMONG MEASUREMENTS ##################################

data.temp = mash_set_data(input_mash$random.b,
                          input_mash$random.s,
                          alpha = 0)

U.c = cov_canonical(data.temp)

V.em = mash_estimate_corr_em(data.temp, U.c, details = TRUE)

Vhat = V.em$V

rm(data.temp)
################################################################################




# SETUP STRONG AND RANDOM SETS #################################################

data.random = mash_set_data(input_mash$random.b,
                            input_mash$random.s,
                            V = Vhat,
                            alpha = 0)

data.strong = mash_set_data(input_mash$strong.b,
                            input_mash$strong.s,
                            V = Vhat,
                            alpha = 0)
################################################################################




# DATA-DRIVEN COVARIANCES FROM STRONG SET ######################################

n_conditions <- ncol(data.strong$Bhat)
cat("Number of conditions:", n_conditions, "\n")

n_pcs <- min(5, n_conditions)

U.pca = cov_pca(data.strong, n_pcs)

U.f = cov_flash(data.strong, factors = "nonneg", tag = "non_neg")
################################################################################




# RUN EXTREME DECONVOLUTION ALGORITHM ##########################################

X.center = apply(data.strong$Bhat, 2, function(x) x - mean(x))
U.emp <- list("XX" = t(X.center) %*% X.center / nrow(X.center))

U.ed <- cov_ed(data.strong, c(U.f, U.pca, U.emp), maxiter = 1000)
################################################################################




# FIT MASH MODEL ON LARGE RANDOM SUBSET WITH COVARIANCES #######################

m = mash(data.random, Ulist = c(U.ed, U.c), outputlevel = 2)

m2 = mash(data.strong, g = get_fitted_g(m), fixg = TRUE)
################################################################################




# SAVE M2 POSTERIOR SUMMARIES ##################################################

# Default filename is RUN_ID-tagged (not date-tagged) as of 2026-08-11 --
# a date-based default silently collides/overwrites when two different
# RUN_IDs are mashed on the same day (confirmed near-miss: 20260810_133901
# and 20260810_142602 both mashed 2026-08-11 would have both defaulted to
# m2_dQTL_250kb_11Aug26.rds). Still overridable via an explicit CLI arg.
args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args) >= 1) args[1] else
  file.path(output_dir, paste0("m2_dQTL_250kb_", run_id, ".rds"))
saveRDS(m2, file = out_file)
cat("\nSaved:", out_file, "\n")

cat("\nSCRIPT COMPLETE\n")
################################################################################




# POST-RUN DIAGNOSTIC WARNINGS #################################################

diagnostics <- local({
  warn <- character()
  add_warn <- function(msg) { warn <<- c(warn, msg); invisible(NULL) }

  .eigvals <- function(M) {
    M <- (M + t(M)) / 2
    tryCatch(eigen(M, symmetric = TRUE, only.values = TRUE)$values,
             error = function(e) rep(NA_real_, ncol(M)))
  }

  .rank_from_eigs <- function(ev, rel_tol = 1e-8) {
    if (all(is.na(ev))) return(NA_integer_)
    lam_max <- max(ev, na.rm = TRUE)
    tol <- rel_tol * max(1, lam_max)
    sum(ev > tol)
  }

  .psd_issues <- function(ev) {
    if (all(is.na(ev))) return(list(bad = TRUE, near = FALSE, min_ev = NA))
    min_ev <- min(ev, na.rm = TRUE)
    list(bad = is.finite(min_ev) && min_ev < -1e-8,
         near = is.finite(min_ev) && min_ev >= -1e-8 && min_ev < 1e-8,
         min_ev = min_ev)
  }

  K <- tryCatch(ncol(data.strong$Bhat), error = function(e) NA_integer_)

  if (exists("U.emp")) {
    Ue <- U.emp[[1]]
    if (!all(is.finite(Ue))) {
      add_warn("U.emp contains non-finite entries; ED may be unstable.")
    }
    evUe <- .eigvals(Ue)
    rUe <- .rank_from_eigs(evUe)
    if (!is.na(K) && !is.na(rUe) && rUe < K) {
      add_warn(sprintf("U.emp rank %d < K = %d; ED may be confined to a low-dimensional subspace.", rUe, K))
    }
  } else {
    add_warn("U.emp not found; ED initializers may be low-rank only.")
  }

  Ui <- c(if (exists("U.f")) U.f else list(),
          if (exists("U.pca")) U.pca else list(),
          if (exists("U.emp")) U.emp else list())
  if (length(Ui) > 0L && !is.na(K)) {
    sumU <- Reduce(`+`, lapply(Ui, function(M) (M + t(M)) / 2))
    evSum <- .eigvals(sumU)
    rSum <- .rank_from_eigs(evSum)
    if (!is.na(rSum) && rSum < K) {
      add_warn(sprintf("Initializer span rank %d < K = %d; ED may be trapped in a subspace.", rSum, K))
    }
  }

  if (exists("Vhat")) {
    evV <- .eigvals(Vhat)
    psdV <- .psd_issues(evV)
    if (psdV$bad) {
      add_warn(sprintf("Vhat has negative eigenvalue (min = %.3e); residual correlation may be non-PSD.", psdV$min_ev))
    } else if (psdV$near) {
      add_warn("Vhat is nearly singular (min eigenvalue < 1e-8).")
    }
    Vsym <- (Vhat + t(Vhat)) / 2
    od <- Vsym - diag(diag(Vsym))
    max_off <- suppressWarnings(max(abs(od), na.rm = TRUE))
    if (is.finite(max_off) && max_off > 0.5) {
      add_warn(sprintf("Large residual correlations (max |offdiag| = %.3f); random subset may not be mostly null.", max_off))
    }
  } else {
    add_warn("Vhat not found; cannot assess residual correlation.")
  }

  if (exists("data.random")) {
    Br <- data.random$Bhat
    Sr <- data.random$Shat
    nr <- nrow(Br)
    idx <- if (!is.null(nr) && nr > 2e5) sample.int(nr, 2e5) else seq_len(nr)
    Zr <- Br[idx, , drop = FALSE] / Sr[idx, , drop = FALSE]
    frac_big <- mean(abs(Zr) > 4, na.rm = TRUE)
    if (is.finite(frac_big) && frac_big > 0.05) {
      add_warn(sprintf("Random set may not be mostly null: Pr(|Z|>4) ≈ %.3f.", frac_big))
    }
  }

  if (exists("data.strong")) {
    Bs <- data.strong$Bhat
    Cs <- suppressWarnings(stats::cor(Bs, use = "pairwise.complete.obs"))
    if (is.matrix(Cs) && all(dim(Cs) == c(K, K))) {
      up <- Cs[upper.tri(Cs)]
      neg_frac <- mean(up < -0.10, na.rm = TRUE)
      min_c <- suppressWarnings(min(up, na.rm = TRUE))
      if (is.finite(neg_frac) && is.finite(min_c) &&
          (neg_frac > 0.25 || min_c < -0.30)) {
        add_warn(sprintf("Substantial antagonism among conditions (neg pair fraction ≈ %.2f, min corr ≈ %.2f); non-negative factors cannot capture such patterns.", neg_frac, min_c))
      }
    }
  }

  a_r <- tryCatch(data.random$alpha, error = function(e) NA_real_)
  a_s <- tryCatch(data.strong$alpha, error = function(e) NA_real_)
  if (is.finite(a_r) && is.finite(a_s) && (a_r != 0 || a_s != 0)) {
    add_warn(sprintf("alpha not EE: data.random=%.3g, data.strong=%.3g.", a_r, a_s))
  }

  if (exists("U.ed")) {
    bad <- 0L; near <- 0L
    for (j in seq_along(U.ed)) {
      ev <- .eigvals(U.ed[[j]])
      chk <- .psd_issues(ev)
      if (chk$bad) bad <- bad + 1L else if (chk$near) near <- near + 1L
    }
    if (bad > 0L)
      add_warn(sprintf("%d ED covariance(s) non-PSD.", bad))
    if (near > 0L)
      add_warn(sprintf("%d ED covariance(s) nearly singular (min eig < 1e-8).", near))
  }

  if (exists("Vhat") && exists("data.strong") && is.finite(K)) {
    if (!all(dim(Vhat) == c(K, K))) {
      add_warn("Dimension mismatch: dim(Vhat) != K x K.")
    }
  }

  list(warn = warn)
})

if (length(diagnostics$warn)) {
  cat("\n=== POST-RUN DIAGNOSTIC WARNINGS ===\n")
  for (w in diagnostics$warn) cat("* ", w, "\n", sep = "")
  cat("=== END WARNINGS ===\n")
} else {
  cat("\nNo diagnostic warnings triggered.\n")
}
################################################################################
