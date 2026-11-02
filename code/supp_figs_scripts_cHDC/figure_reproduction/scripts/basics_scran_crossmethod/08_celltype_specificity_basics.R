stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(ggplot2)
})

# ---------------------------------------------------------------------------
# NAMED CONSTANTS
# ---------------------------------------------------------------------------
MIN_SHARED_GENES <- 10L  # too few shared genes => alignment likely broken (symbol mismatch)

CELL_TYPES <- c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                 "Epithelial_cells", "Activated_fibroblasts",
                 "Foregut_endodermal_cells")
N_CT <- length(CELL_TYPES)

# ANALYSIS_OK[fixed-design]: the 3 pooled iPSC lines are the fixed experimental
# design of the cHDC dataset (not a tunable exclusion list); matches
# compare_1000subset_to_memento.R's INDIVIDUALS constant.
INDIVIDUALS <- c("NA19099", "NA19093", "NA18522")

# PATCHED (figure_reproduction/PORT_CHANGELOG_crossmethod.md): EPS_DIR and
# MEMENTO_MATCHED_DIR point at this fork's mirror directories (5 cell types
# symlinked from the original analysis, Epithelial_cells rebuilt
# excl-cluster24 -- see setup_mirror_dirs.sh / 02_new_epithelial_1000subset_
# memento.py); MEMENTO_FULL_CSV points at the canonical excl-cluster24 matrix
# (already has all 6 cell types in one file -- only Epithelial's 3 columns
# differ from published). Model/statistics code below is otherwise
# byte-identical to analysis/basics_dispersion_reanalysis/scripts/
# celltype_specificity_matched_vs_unmatched.R.
EPS_DIR <- "figure_reproduction/outputs/basics_scran_crossmethod/epsilon"
MEMENTO_MATCHED_DIR <- "figure_reproduction/outputs/basics_scran_crossmethod/memento_1000subset"
MEMENTO_FULL_CSV <- "figure_reproduction/data/res_var_excl_cluster24.csv"
OUT_DIR <- "figure_reproduction/outputs/basics_scran_crossmethod/comparison/celltype_specificity"
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# LOAD: BASiCS epsilon (1000-cell subset), one named vector per cell type
# ---------------------------------------------------------------------------
load_basics_vec <- function(cell_type) {
  f <- file.path(EPS_DIR, paste0(cell_type, "_1000subset_epsilon.csv"))
  stopifnot("BASiCS 1000-subset epsilon CSV not found -- run 02 for all 6 cell types first" =
              file.exists(f))
  eps <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot("epsilon CSV missing 'gene' column" = "gene" %in% colnames(eps))
  num_cols <- setdiff(colnames(eps)[vapply(eps, is.numeric, logical(1))], "gene")
  stopifnot("no numeric epsilon columns" = length(num_cols) > 0)
  eps_col <- if ("epsilon" %in% num_cols) "epsilon" else num_cols[1]
  vec <- setNames(eps[[eps_col]], eps$gene)
  n_dropped <- sum(!is.finite(vec))
  message("  BASiCS[", cell_type, "] (col='", eps_col, "'): ", length(vec),
          " genes, dropping ", n_dropped, " non-finite")
  vec[is.finite(vec)]
}

# ---------------------------------------------------------------------------
# LOAD: Memento residual variance, one named vector per cell type, per
# comparator ("matched_1000cell" = same 3000-cell subset per cell type;
# "fullcell" = full cell-type cell set).
# ---------------------------------------------------------------------------
load_memento_matched_vec <- function(cell_type) {
  f <- file.path(MEMENTO_MATCHED_DIR, paste0("dispersion_by_ind_", cell_type, ".csv"))
  stopifnot("matched-cell-count Memento CSV not found" = file.exists(f))
  df <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot("matched Memento CSV missing 'gene' column" = "gene" %in% colnames(df))
  cols <- paste0(INDIVIDUALS, "_", cell_type)
  present <- cols %in% colnames(df)
  stopifnot("no matched-Memento columns found for this cell type" = any(present))
  vec <- setNames(rowMeans(df[, cols[present], drop = FALSE], na.rm = TRUE), df$gene)
  n_dropped <- sum(!is.finite(vec))
  message("  Memento-matched[", cell_type, "]: ", length(vec),
          " genes (", sum(present), "/", length(cols), " individuals present), dropping ",
          n_dropped, " non-finite")
  vec[is.finite(vec)]
}

message("Loading full-cell Memento matrix once (shared across all fullcell vectors)...")
mem_full_mat <- read.csv(MEMENTO_FULL_CSV, row.names = 1, check.names = FALSE)
stopifnot("full-cell Memento matrix empty" = nrow(mem_full_mat) > 0 && ncol(mem_full_mat) > 0)

load_memento_full_vec <- function(cell_type) {
  cols <- paste0(INDIVIDUALS, "_", cell_type)
  present <- cols %in% colnames(mem_full_mat)
  stopifnot("no full-cell-count Memento columns found for this cell type" = any(present))
  vec <- rowMeans(mem_full_mat[, cols[present], drop = FALSE], na.rm = TRUE)
  n_dropped <- sum(!is.finite(vec))
  message("  Memento-fullcell[", cell_type, "]: ", length(vec),
          " genes (", sum(present), "/", length(cols), " individuals present), dropping ",
          n_dropped, " non-finite")
  vec[is.finite(vec)]
}

message("=== Loading BASiCS epsilon (all ", N_CT, " cell types) ===")
basics_list <- setNames(lapply(CELL_TYPES, load_basics_vec), CELL_TYPES)

message("=== Loading Memento (matched_1000cell comparator) ===")
memento_matched_list <- setNames(lapply(CELL_TYPES, load_memento_matched_vec), CELL_TYPES)

message("=== Loading Memento (fullcell comparator) ===")
memento_fullcell_list <- setNames(lapply(CELL_TYPES, load_memento_full_vec), CELL_TYPES)

# ---------------------------------------------------------------------------
# Exact permutation generator: all n! permutations of 1:n as an (n!) x n
# integer matrix. n = 6 -> 720 rows, trivial to enumerate and hold in memory.
# ---------------------------------------------------------------------------
all_perms <- function(n) {
  stopifnot(n >= 1L)
  if (n == 1L) return(matrix(1L, nrow = 1L, ncol = 1L))
  sub <- all_perms(n - 1L)
  out <- matrix(NA_integer_, nrow = nrow(sub) * n, ncol = n)
  row_i <- 1L
  for (i in seq_len(nrow(sub))) {
    for (pos in seq_len(n)) {
      out[row_i, ] <- append(sub[i, ], n, after = pos - 1L)
      row_i <- row_i + 1L
    }
  }
  out
}
PERM_MATRIX <- all_perms(N_CT)
stopifnot("permutation enumeration count mismatch" = nrow(PERM_MATRIX) == factorial(N_CT))
message("Exact permutation null: ", nrow(PERM_MATRIX), " label permutations of ", N_CT, " cell types")

# ---------------------------------------------------------------------------
# Build the full BASiCS-cell-type x Memento-cell-type Spearman rho matrix for
# one comparator, run the exact permutation test, and write matrix/heatmap/
# boxplot outputs. Returns the per-cell-type detail rows and one summary row.
# ---------------------------------------------------------------------------
run_specificity_check <- function(memento_list, comparator_label) {
  message("\n=== Comparator: ", comparator_label, " ===")
  rho_mat <- matrix(NA_real_, nrow = N_CT, ncol = N_CT,
                     dimnames = list(CELL_TYPES, CELL_TYPES))
  n_shared_mat <- matrix(NA_integer_, nrow = N_CT, ncol = N_CT,
                          dimnames = list(CELL_TYPES, CELL_TYPES))

  for (bi in seq_len(N_CT)) {
    for (mi in seq_len(N_CT)) {
      bt <- CELL_TYPES[bi]; mt <- CELL_TYPES[mi]
      x <- basics_list[[bt]]; y <- memento_list[[mt]]
      shared <- intersect(names(x), names(y))
      stopifnot("fewer than MIN_SHARED_GENES shared genes -- alignment likely broken" =
                  length(shared) >= MIN_SHARED_GENES)
      # ANALYSIS_OK[warning-suppression]: cor() emits a zero-variance warning
      # only in a degenerate all-tied-ranks edge case, which does not occur
      # for real epsilon/residual-variance data; guards against that edge
      # case without masking a real computation error (rho would be NA, and
      # the downstream stopifnot on rho_mat range below would catch it).
      rho_mat[bi, mi] <- suppressWarnings(cor(x[shared], y[shared], method = "spearman"))
      n_shared_mat[bi, mi] <- length(shared)
    }
  }
  stopifnot("Spearman rho matrix has non-finite entries" = all(is.finite(rho_mat)))
  stopifnot("Spearman rho matrix has out-of-range entries" = all(rho_mat >= -1 & rho_mat <= 1))
  message("  rho matrix: ", N_CT, "x", N_CT, ", shared genes per pair range [",
          min(n_shared_mat), ", ", max(n_shared_mat), "]")

  write.csv(rho_mat, file.path(OUT_DIR, paste0("rho_matrix_", comparator_label, ".csv")))

  # ---- Exact permutation test on mean(diagonal under permuted column labels) ----
  # ANALYSIS_OK[positional-index]: rho_mat's rows/cols are already identifier-
  # aligned to CELL_TYPES (set as dimnames above; each cell was itself computed
  # from a name-based gene intersection, asserted via MIN_SHARED_GENES). This
  # permutation test deliberately enumerates all bijections of that fixed
  # index set onto itself -- integer position IS the identifier here, and a
  # "shared identifier join" would be a no-op since the index set never
  # changes across permutations.
  null_means <- apply(PERM_MATRIX, 1, function(p) mean(rho_mat[cbind(seq_len(N_CT), p)]))
  observed_mean <- mean(diag(rho_mat))
  identity_row <- which(apply(PERM_MATRIX, 1, function(p) all(p == seq_len(N_CT))))
  stopifnot("identity permutation not found in enumerated set" = length(identity_row) == 1)
  stopifnot("identity permutation's null mean does not match observed diagonal mean" =
              isTRUE(all.equal(null_means[identity_row], observed_mean)))
  p_value <- mean(null_means >= observed_mean)
  message(sprintf("  Observed mean matched (diagonal) rho = %.4f", observed_mean))
  message(sprintf("  Null (any pairing) mean rho: mean=%.4f, range=[%.4f, %.4f]",
                  mean(null_means), min(null_means), max(null_means)))
  message(sprintf("  Exact permutation p-value (matched >= chance pairing) = %.4f (n=%d perms)",
                  p_value, nrow(PERM_MATRIX)))

  # ---- Per-cell-type spot-check detail: rank of the matched entry within its row ----
  detail_rows <- lapply(seq_len(N_CT), function(bi) {
    row <- rho_mat[bi, ]
    off_diag <- row[-bi]
    rank_of_matched <- N_CT - rank(row, ties.method = "average")[bi] + 1L  # 1 = highest rho in row
    data.frame(
      comparator = comparator_label,
      basics_cell_type = CELL_TYPES[bi],
      matched_rho = row[bi],
      mean_unmatched_rho = mean(off_diag),
      min_unmatched_rho = min(off_diag),
      max_unmatched_rho = max(off_diag),
      rank_of_matched_of_n = rank_of_matched,
      n_cell_types = N_CT,
      stringsAsFactors = FALSE
    )
  })
  detail_df <- do.call(rbind, detail_rows)
  n_rank1 <- sum(detail_df$rank_of_matched_of_n == 1)
  message("  Matched cell type is the single best Memento correlate for ", n_rank1,
          "/", N_CT, " BASiCS cell types")

  # ---- Heatmap ----
  hm_df <- as.data.frame(as.table(rho_mat))
  colnames(hm_df) <- c("basics_cell_type", "memento_cell_type", "rho")
  hm_df$is_matched <- hm_df$basics_cell_type == hm_df$memento_cell_type
  p_hm <- ggplot(hm_df, aes(memento_cell_type, basics_cell_type, fill = rho)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", rho), fontface = ifelse(is_matched, "bold", "plain")),
              size = 3) +
    # ANALYSIS_OK[magic-threshold]: rho_mat is all-positive Spearman rho
    # (observed range ~0.82-0.96 across both comparators -- see
    # rho_matrix_{comparator}.csv), so the diverging blue-white-red scale
    # centered at 0 wasted its blue half and compressed all contrast into a
    # narrow red band. Fixed sequential limits at [0.8, 1] match that
    # observed range; oob = scales::squish clips any future out-of-range
    # cell to the nearest end color rather than rendering as NA/grey.
    scale_fill_gradient(low = "white", high = "red", limits = c(0.8, 1),
                         oob = scales::squish) +
    labs(title = paste0("BASiCS epsilon vs Memento residual variance (", comparator_label, ")"),
         subtitle = "Bold = matched (same) cell type; rows = BASiCS, cols = Memento",
         x = "Memento cell type", y = "BASiCS cell type", fill = "Spearman\nrho") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # PDF (not PNG) per Brendan's request, 2026-09-02 -- vector output for
  # editing in Affinity Designer; cairo_pdf matches the repo's device
  # convention for vector figures elsewhere.
  ggsave(file.path(FIG_DIR, paste0("heatmap_", comparator_label, ".pdf")),
         p_hm, width = 7.5, height = 6.5, device = cairo_pdf)

  # ---- Matched vs unmatched boxplot ----
  box_df <- rbind(
    data.frame(group = "matched", rho = diag(rho_mat)),
    data.frame(group = "unmatched", rho = rho_mat[row(rho_mat) != col(rho_mat)])
  )
  p_box <- ggplot(box_df, aes(group, rho)) +
    geom_boxplot(outlier.shape = NA, width = 0.4) +
    geom_jitter(width = 0.08, alpha = 0.5, size = 1.5) +
    labs(title = paste0("Matched vs. unmatched cell-type rho (", comparator_label, ")"),
         subtitle = sprintf("Exact permutation p = %.4f (n=%d cell types, %d perms)",
                            p_value, N_CT, nrow(PERM_MATRIX)),
         x = NULL, y = "Spearman rho (BASiCS epsilon vs Memento)") +
    theme_bw()
  # PDF (not PNG) per Brendan's request, 2026-09-02 -- see heatmap ggsave
  # above for rationale. Width bumped from the PNG version's 5in to 6.5in --
  # at 5in the "matched_1000cell" comparator's longer title clipped off the
  # right edge of the canvas (confirmed by rendering the saved PDF).
  ggsave(file.path(FIG_DIR, paste0("matched_vs_unmatched_", comparator_label, ".pdf")),
         p_box, width = 6.5, height = 5, device = cairo_pdf)

  list(
    detail = detail_df,
    summary = data.frame(
      comparator = comparator_label,
      n_cell_types = N_CT,
      mean_matched_rho = observed_mean,
      median_matched_rho = median(diag(rho_mat)),
      mean_unmatched_rho = mean(rho_mat[row(rho_mat) != col(rho_mat)]),
      median_unmatched_rho = median(rho_mat[row(rho_mat) != col(rho_mat)]),
      diff_matched_minus_unmatched = observed_mean - mean(rho_mat[row(rho_mat) != col(rho_mat)]),
      n_matched_rank1_of_n = n_rank1,
      permutation_p_value = p_value,
      n_permutations = nrow(PERM_MATRIX),
      stringsAsFactors = FALSE
    )
  )
}

result_matched <- run_specificity_check(memento_matched_list, "matched_1000cell")
result_fullcell <- run_specificity_check(memento_fullcell_list, "fullcell")

detail_all <- rbind(result_matched$detail, result_fullcell$detail)
summary_all <- rbind(result_matched$summary, result_fullcell$summary)

write.csv(detail_all, file.path(OUT_DIR, "per_celltype_detail.csv"), row.names = FALSE)
write.csv(summary_all, file.path(OUT_DIR, "specificity_summary.csv"), row.names = FALSE)

message("\nDone. Outputs written to ", OUT_DIR)
print(summary_all)
