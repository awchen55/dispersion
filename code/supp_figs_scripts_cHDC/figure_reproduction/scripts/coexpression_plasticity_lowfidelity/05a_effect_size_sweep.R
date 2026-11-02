#!/usr/bin/env Rscript

stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages(library(WGCNA))

BASE_DIR <- "figure_reproduction/outputs/coexpression_plasticity_lowfidelity"
ORIG_MODULES_DIR <- "analysis/coexpression_plasticity_lowfidelity/outputs/cscore_modules"
OUT_CSV  <- file.path(BASE_DIR, "diagnostics", "effect_size_sweep_all_celltypes.csv")

CELL_TYPES <- c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells")
INDIVIDUALS <- c("NA19099", "NA19093", "NA18522")
EFFECT_SIZE_FLOORS <- c(0.1, 0.15, 0.2, 0.25, 0.3)  # matches source script 07
POWER_VECTOR <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))  # matches source scripts 05-14
RSQUARED_CUT <- 0.85

log_msg <- function(...) message("[05a_effect_size_sweep] ", ...)

rows <- list()
for (ct in CELL_TYPES) {
  in_dir <- if (ct == "Epithelial_cells") file.path(BASE_DIR, "cscore_modules") else ORIG_MODULES_DIR
  for (ind in INDIVIDUALS) {
    rds <- file.path(in_dir, paste0(ct, "_", ind, "_coexp_thresholded.rds"))
    if (!file.exists(rds)) stop("Missing expected file: ", rds)
    sim_base <- abs(readRDS(rds))
    for (floor in EFFECT_SIZE_FLOORS) {
      sim <- sim_base
      sim[sim < floor] <- 0
      diag(sim) <- 1
      sft <- WGCNA::pickSoftThreshold.fromSimilarity(
        similarity = sim, RsquaredCut = RSQUARED_CUT, powerVector = POWER_VECTOR, verbose = 0
      )
      fit_tab <- sft$fitIndices
      crosses <- fit_tab$SFT.R.sq >= RSQUARED_CUT
      best_idx <- which.max(fit_tab$SFT.R.sq)
      chosen_idx <- if (any(crosses)) which(crosses)[1] else best_idx
      rows[[paste(ct, ind, floor)]] <- data.frame(
        cell_type = ct, individual = ind, floor = floor,
        chosen_power = fit_tab$Power[chosen_idx], chosen_r2 = fit_tab$SFT.R.sq[chosen_idx],
        crosses_at_all = any(crosses),
        best_power = fit_tab$Power[best_idx], best_r2 = fit_tab$SFT.R.sq[best_idx]
      )
      log_msg(ct, " ", ind, " floor=", floor, ": chosen power=", fit_tab$Power[chosen_idx],
              " R2=", round(fit_tab$SFT.R.sq[chosen_idx], 3), " (best R2=", round(fit_tab$SFT.R.sq[best_idx], 3), ")")
      rm(sim, sft, fit_tab)
    }
    rm(sim_base); gc(FALSE)
  }
}
out <- do.call(rbind, rows); rownames(out) <- NULL
stopifnot("expected 18 samples x 5 floors" = nrow(out) == length(CELL_TYPES) * length(INDIVIDUALS) * length(EFFECT_SIZE_FLOORS),
          "NA R^2 in sweep" = !anyNA(out$chosen_r2))
write.csv(out, OUT_CSV, row.names = FALSE)
log_msg("Wrote: ", OUT_CSV, "\n[05a_effect_size_sweep] DONE.")
