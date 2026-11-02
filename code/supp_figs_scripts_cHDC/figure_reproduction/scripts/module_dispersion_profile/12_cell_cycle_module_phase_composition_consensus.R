#!/usr/bin/env Rscript

stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(Seurat)
})

set.seed(42)

CONSENSUS_DIR <- "figure_reproduction/outputs/coexpression_plasticity_lowfidelity/consensus_modules_prototype"
OUT_DIR       <- "figure_reproduction/outputs/module_dispersion_profile/phase2_consensus"
CONSENSUS_QUANTILE_SUFFIX <- "q0.5"

CELL_TYPES <- c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                "Epithelial_cells", "Activated_fibroblasts",
                "Foregut_endodermal_cells")
N_PERM <- 10000L
LARGE_MODULE_THRESHOLD <- 10L

log_msg <- function(...) message("[12_cell_cycle_module_phase_composition_consensus] ", ...)

cc_gene_list <- Seurat::cc.genes.updated.2019
stopifnot("S and G2M gene lists unexpectedly overlap" =
            length(intersect(cc_gene_list$s.genes, cc_gene_list$g2m.genes)) == 0)
phase_lookup_all <- setNames(
  c(rep("S", length(cc_gene_list$s.genes)), rep("G2M", length(cc_gene_list$g2m.genes))),
  c(cc_gene_list$s.genes, cc_gene_list$g2m.genes)
)
log_msg(length(cc_gene_list$s.genes), " S-phase + ", length(cc_gene_list$g2m.genes),
        " G2M-phase = ", length(phase_lookup_all), " unique genes (no G1 list in this source)")

same_phase_pair_counts <- function(module_ids, phases) {
  stopifnot(length(module_ids) == length(phases))
  same <- 0L
  total <- 0L
  for (m in unique(module_ids)) {
    p <- phases[module_ids == m]
    n_s <- sum(p == "S")
    n_g <- sum(p == "G2M")
    n_tot <- n_s + n_g
    same <- same + choose(n_s, 2) + choose(n_g, 2)
    total <- total + choose(n_tot, 2)
  }
  list(same = same, total = total)
}

composition_rows <- list()
test_rows <- list()

for (ct in CELL_TYPES) {
  log_msg("=== ", ct, " ===")
  module_csv <- file.path(CONSENSUS_DIR, paste0(ct, "_consensus_", CONSENSUS_QUANTILE_SUFFIX, "_module_membership.csv"))
  if (!file.exists(module_csv)) stop("Missing consensus module membership file: ", module_csv)
  module_df <- read.csv(module_csv, stringsAsFactors = FALSE)
  stopifnot(
    "consensus module_membership.csv missing expected columns" =
      all(c("gene", "module", "in_module") %in% colnames(module_df)),
    "consensus module_membership.csv has duplicate genes" = !anyDuplicated(module_df$gene)
  )

  cc_found <- intersect(names(phase_lookup_all), module_df$gene)
  if (length(cc_found) == 0) {
    stop("[", ct, "] zero cell-cycle genes found in this cell type's consensus universe -- unexpected, investigate")
  }

  module_lookup <- setNames(module_df$module, module_df$gene)
  in_module_lookup <- setNames(module_df$in_module, module_df$gene)
  is_clustered <- in_module_lookup[cc_found] & module_lookup[cc_found] != 0
  clustered_genes <- cc_found[is_clustered]
  n_unassigned <- length(cc_found) - length(clustered_genes)
  log_msg("  cc genes found: ", length(cc_found), "; clustered (in a real module): ",
          length(clustered_genes), "; unassigned/module-0: ", n_unassigned)

  if (length(clustered_genes) == 0) {
    log_msg("  WARNING: no clustered cell-cycle genes -- phase composition/test skipped for ", ct)
    next
  }

  clustered_modules <- module_lookup[clustered_genes]
  clustered_phases <- phase_lookup_all[clustered_genes]

  comp_df <- do.call(rbind, lapply(sort(unique(clustered_modules)), function(m) {
    p <- clustered_phases[clustered_modules == m]
    n_s <- sum(p == "S"); n_g <- sum(p == "G2M"); n_tot <- n_s + n_g
    dominant <- if (n_s > n_g) "S" else if (n_g > n_s) "G2M" else "tied"
    data.frame(cell_type = ct, module = m,
               n_S = n_s, n_G2M = n_g, n_total_cc = n_tot,
               dominant_phase = dominant, purity = max(n_s, n_g) / n_tot,
               is_large_module = n_tot > LARGE_MODULE_THRESHOLD,
               stringsAsFactors = FALSE)
  }))
  composition_rows[[ct]] <- comp_df

  n_large <- sum(comp_df$is_large_module)
  log_msg("  ", nrow(comp_df), " cell-cycle-containing module(s); ", n_large,
          " with >", LARGE_MODULE_THRESHOLD, " cc genes: ",
          paste(sprintf("module %d (%s, %d/%d %s, purity=%.2f)",
                         comp_df$module[comp_df$is_large_module],
                         comp_df$dominant_phase[comp_df$is_large_module],
                         pmax(comp_df$n_S[comp_df$is_large_module], comp_df$n_G2M[comp_df$is_large_module]),
                         comp_df$n_total_cc[comp_df$is_large_module],
                         comp_df$dominant_phase[comp_df$is_large_module],
                         comp_df$purity[comp_df$is_large_module]),
                collapse = "; "))

  obs <- same_phase_pair_counts(clustered_modules, clustered_phases)
  if (obs$total == 0) {
    log_msg("  WARNING: every cell-cycle-containing module has <2 clustered cc genes -- ",
            "no same-module pairs exist; phase-segregation test skipped for ", ct)
    test_rows[[ct]] <- data.frame(
      cell_type = ct, n_clustered_cc_genes = length(clustered_genes),
      n_cc_containing_modules = length(unique(clustered_modules)),
      n_same_module_pairs = 0L, observed_frac_same_phase = NA_real_,
      perm_null_mean = NA_real_, perm_p_segregation = NA_real_, n_permutations = N_PERM,
      stringsAsFactors = FALSE
    )
    next
  }
  observed_frac <- obs$same / obs$total

  null_fracs <- vapply(seq_len(N_PERM), function(i) {
    perm_phases <- sample(clustered_phases)
    same_phase_pair_counts(clustered_modules, perm_phases)$same / obs$total
  }, numeric(1))
  perm_p <- (1 + sum(null_fracs >= observed_frac)) / (N_PERM + 1)

  log_msg(sprintf("  same-module pairs: %d, observed same-phase frac: %.3f, perm null mean: %.3f, perm p: %.4f",
                  obs$total, observed_frac, mean(null_fracs), perm_p))

  test_rows[[ct]] <- data.frame(
    cell_type = ct, n_clustered_cc_genes = length(clustered_genes),
    n_cc_containing_modules = length(unique(clustered_modules)),
    n_same_module_pairs = obs$total, observed_frac_same_phase = observed_frac,
    perm_null_mean = mean(null_fracs), perm_p_segregation = perm_p, n_permutations = N_PERM,
    stringsAsFactors = FALSE
  )
}

composition_df <- do.call(rbind, composition_rows); rownames(composition_df) <- NULL
composition_csv <- file.path(OUT_DIR, "cc_module_phase_composition_consensus.csv")
write.csv(composition_df, composition_csv, row.names = FALSE)
log_msg("Wrote: ", composition_csv, " (", nrow(composition_df), " module rows across ",
        length(composition_rows), " cell types)")

test_df <- do.call(rbind, test_rows); rownames(test_df) <- NULL
stopifnot("expected one test row per cell type" = nrow(test_df) == length(CELL_TYPES))
test_csv <- file.path(OUT_DIR, "cc_module_phase_segregation_test_consensus.csv")
write.csv(test_df, test_csv, row.names = FALSE)
log_msg("Wrote: ", test_csv)
print(test_df)

large_df <- composition_df[composition_df$is_large_module, ]
log_msg("=== Large (>", LARGE_MODULE_THRESHOLD, "-gene) cell-cycle-containing modules only ===")
print(large_df[, c("cell_type", "module", "n_S", "n_G2M", "n_total_cc", "dominant_phase", "purity")])

log_msg("DONE.")
