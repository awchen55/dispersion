stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(ggplot2)
})

# ---------------------------------------------------------------------------
# NAMED CONSTANTS
# ---------------------------------------------------------------------------
CELL_TYPES <- c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                "Epithelial_cells", "Activated_fibroblasts",
                "Foregut_endodermal_cells")
N_CT <- length(CELL_TYPES)
N_QUANTILES <- 5L               # source Rmd line 96
N_CT_REQUIRED <- N_CT           # source Rmd: rowSums(. == 1) >= 6 (all 6, not a majority)
RIBO_PATTERN <- "^RP"           # source Rmd line 51 -- the ONLY filter in the canonical definition

# ANALYSIS_OK[fixed-design]: the 3 pooled iPSC lines are the fixed experimental
# design of the cHDC dataset; matches the source Rmd's column-name convention.
INDIVIDUALS <- c("NA19099", "NA19093", "NA18522")

# PATCHED (figure_reproduction/PORT_CHANGELOG.md, revised per Brendan's
# request): the reference Memento gene sets are read DIRECTLY from
# figure_reproduction's own already-staged 02_rerun_ranked_dispersion.R
# output -- ribosomal genes dropped AND Alexander's -6 top-expressed-gene
# exclusion applied, on the Epithelial-excl-cluster-24 matrix. This is the
# new canonical definition (the "main data for the manuscript" going
# forward), not a from-scratch re-derivation to be validated against
# anything older.
MEMENTO_LOW_TXT  <- "figure_reproduction/outputs/ranked/consistently_low_retained.txt"
MEMENTO_HIGH_TXT <- "figure_reproduction/outputs/ranked/consistently_high_retained.txt"

# PATCHED: mirror directories (5 cell types symlinked, Epithelial_cells new --
# see setup_mirror_dirs.sh) and this fork's freshly-regenerated scran outputs.
SCRAN_FULL_DIR   <- "figure_reproduction/outputs/basics_scran_crossmethod/scran/bio_variance"
BASICS_SUBSET_DIR <- "figure_reproduction/outputs/basics_scran_crossmethod/epsilon"
MEMENTO_SUBSET_DIR <- "figure_reproduction/outputs/basics_scran_crossmethod/memento_1000subset"

OUT_DIR <- "figure_reproduction/outputs/basics_scran_crossmethod/consistent_geneset_retention"
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Core algorithm, replicating main_figs_1-4_related_supp.Rmd lines 65-107 and
# 804-816 exactly.
# ---------------------------------------------------------------------------
# ANALYSIS_OK[within-file-helper]: called below; cross-file "dead code" is a
# false positive (no callers in OTHER files by design).
ranks_to_quintiles <- function(ranks) {
  cut(ranks,
      breaks = quantile(ranks, probs = seq(0, 1, length.out = N_QUANTILES + 1), na.rm = TRUE),
      labels = FALSE, include.lowest = TRUE)
}

# Given a genes x N_CT numeric matrix (one dispersion value per gene per cell
# type; every row must be fully present -- callers restrict to their own
# all-6-cell-types-finite gene universe before calling this), return the
# genes x N_CT integer quintile-bin matrix.
# ANALYSIS_OK[within-file-helper]: called below; cross-file "dead code" is a
# false positive (no callers in OTHER files by design).
build_quintile_matrix <- function(mat) {
  stopifnot("input matrix has non-finite values -- restrict to the all-6-cell-types-finite gene universe first" =
              all(is.finite(mat)))
  qmat <- apply(mat, 2, function(x) ranks_to_quintiles(rank(x, ties.method = "average")))
  rownames(qmat) <- rownames(mat)
  qmat
}

# ANALYSIS_OK[within-file-helper]: called below; cross-file "dead code" is a
# false positive (no callers in OTHER files by design).
consistent_genes <- function(quintile_mat, target_quintile) {
  hits <- rowSums(quintile_mat == target_quintile, na.rm = TRUE)
  rownames(quintile_mat)[hits >= N_CT_REQUIRED]
}

# ---------------------------------------------------------------------------
# Load Memento's reference canonical gene sets -- the -6-excluded, Epithelial-
# excl-cluster-24 retained sets already computed by 02_rerun_ranked_
# dispersion.R (see PATCHED note above). No re-derivation, no validation
# against an older baseline: this IS the reference definition.
# ---------------------------------------------------------------------------
message("=== Loading Memento's canonical low/high gene sets (excl-cluster-24, -6-threshold retained) ===")
stopifnot(
  "Memento consistently-low retained gene list not found -- run 02_rerun_ranked_dispersion.R first" =
    file.exists(MEMENTO_LOW_TXT),
  "Memento consistently-high retained gene list not found -- run 02_rerun_ranked_dispersion.R first" =
    file.exists(MEMENTO_HIGH_TXT)
)
memento_low  <- readLines(MEMENTO_LOW_TXT)
memento_high <- readLines(MEMENTO_HIGH_TXT)
stopifnot(
  "Memento consistently-low gene list is empty" = length(memento_low) > 0,
  "Memento consistently-high gene list is empty" = length(memento_high) > 0
)
message("Loaded: ", length(memento_low), " consistently-low, ",
        length(memento_high), " consistently-high genes")

# ---------------------------------------------------------------------------
# Load one candidate method's per-cell-type dispersion measure, restrict to
# the gene universe finite in ALL 6 cell types, and return that genes x
# N_CT matrix.
# ---------------------------------------------------------------------------
# ANALYSIS_OK[within-file-helper]: called below; cross-file "dead code" is a
# false positive (no callers in OTHER files by design).
load_method_matrix <- function(dir, file_suffix, value_col, method_label) {
  vecs <- setNames(lapply(CELL_TYPES, function(ct) {
    f <- file.path(dir, paste0(ct, file_suffix))
    stopifnot("input CSV not found" = file.exists(f))
    df <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
    stopifnot("input CSV missing 'gene' column" = "gene" %in% colnames(df))
    stopifnot("input CSV missing expected value column" = value_col %in% colnames(df))
    v <- setNames(df[[value_col]], df$gene)
    v[is.finite(v)]
  }), CELL_TYPES)
  universe <- Reduce(intersect, lapply(vecs, names))
  message("  [", method_label, "] genes finite in all ", N_CT, " cell types: ", length(universe))
  mat <- vapply(CELL_TYPES, function(ct) vecs[[ct]][universe], numeric(length(universe)))
  rownames(mat) <- universe
  mat
}

# ---------------------------------------------------------------------------
# Load Memento's OWN dispersion re-estimated on the 1000-cell/individual
# subset (the baseline-control candidate). Unlike BASiCS/scran (one value per
# gene per cell type already), this source has 3 per-individual columns per
# cell type, same as the full-scale Memento data -- median across individuals
# first (source Rmd line 88), THEN apply the same ribosomal-gene filter the
# canonical algorithm applies, before restricting to the all-6-cell-types
# gene universe.
# ---------------------------------------------------------------------------
# ANALYSIS_OK[within-file-helper]: called below; cross-file "dead code" is a
# false positive (no callers in OTHER files by design).
load_memento_subset_matrix <- function() {
  vecs <- setNames(lapply(CELL_TYPES, function(ct) {
    f <- file.path(MEMENTO_SUBSET_DIR, paste0("dispersion_by_ind_", ct, ".csv"))
    stopifnot("Memento 1000-cell subset dispersion CSV not found" = file.exists(f))
    df <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
    stopifnot("Memento subset CSV missing 'gene' column" = "gene" %in% colnames(df))
    cols <- paste0(INDIVIDUALS, "_", ct)
    stopifnot("expected individual columns missing for this cell type" = all(cols %in% colnames(df)))
    med <- apply(df[, cols, drop = FALSE], 1, median)
    v <- setNames(med, df$gene)
    v <- v[is.finite(v)]
    v[!grepl(RIBO_PATTERN, names(v))]
  }), CELL_TYPES)
  universe <- Reduce(intersect, lapply(vecs, names))
  message("  [memento_1000subset] genes finite (ribo-filtered) in all ", N_CT, " cell types: ", length(universe))
  mat <- vapply(CELL_TYPES, function(ct) vecs[[ct]][universe], numeric(length(universe)))
  rownames(mat) <- universe
  mat
}

# ---------------------------------------------------------------------------
# Retention: of a Memento canonical gene set, what fraction (a) survive the
# method's own gene filter and are present in its all-6-cell-types universe,
# and (b) are ALSO classified the same way (consistently low/high) by that
# method's own quintile ranking.
# ---------------------------------------------------------------------------
# ANALYSIS_OK[within-file-helper]: called below; cross-file "dead code" is a
# false positive (no callers in OTHER files by design).
retention_stats <- function(memento_set, universe, method_set, method_label, geneset_label) {
  n_memento <- length(memento_set)
  n_present <- sum(memento_set %in% universe)
  n_retained <- sum(memento_set %in% method_set)
  data.frame(
    method = method_label, geneset = geneset_label,
    n_memento_geneset = n_memento,
    n_present_in_method_universe = n_present,
    n_retained_consistent = n_retained,
    prop_retained_of_memento_geneset = n_retained / n_memento,
    prop_retained_of_present_genes = if (n_present > 0) n_retained / n_present else NA_real_,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Run all 3 candidates: scran (full cell count), BASiCS (1000-cell subset,
# the only scale available), and Memento itself re-ranked on the SAME
# 1000-cell subset (baseline control, isolates cell-count/subsampling noise
# from genuine method disagreement).
# ---------------------------------------------------------------------------
message("\n=== Loading candidate method matrices ===")
scran_full_mat     <- load_method_matrix(SCRAN_FULL_DIR,    "_bio_variance.csv",       "bio",    "scran_fullcell")
basics_subset_mat  <- load_method_matrix(BASICS_SUBSET_DIR, "_1000subset_epsilon.csv", "median", "basics_1000subset")
memento_subset_mat <- load_memento_subset_matrix()

candidates <- list(
  scran_fullcell        = scran_full_mat,
  basics_1000subset     = basics_subset_mat,
  memento_1000subset    = memento_subset_mat
)

message("\n=== Classifying each method's own consistently-low/high sets ===")
result_rows <- list()
for (method_label in names(candidates)) {
  mat <- candidates[[method_label]]
  qmat <- build_quintile_matrix(mat)
  low_set  <- consistent_genes(qmat, 1L)
  high_set <- consistent_genes(qmat, 5L)
  message("  [", method_label, "] own consistently-low: ", length(low_set),
          " | own consistently-high: ", length(high_set))

  result_rows[[paste0(method_label, "_low")]]  <- retention_stats(
    memento_low,  rownames(mat), low_set,  method_label, "consistently_low")
  result_rows[[paste0(method_label, "_high")]] <- retention_stats(
    memento_high, rownames(mat), high_set, method_label, "consistently_high")
}

summary_df <- do.call(rbind, result_rows)
rownames(summary_df) <- NULL
write.csv(summary_df, file.path(OUT_DIR, "retention_summary.csv"), row.names = FALSE)
message("\nDone. Retention summary written to ", file.path(OUT_DIR, "retention_summary.csv"))
print(summary_df)

# ---- Bar chart: proportion of Memento's canonical gene set retained ----
plot_df <- summary_df
plot_df$geneset <- factor(plot_df$geneset, levels = c("consistently_low", "consistently_high"),
                           labels = c("Consistently low", "Consistently high"))
p <- ggplot(plot_df, aes(method, prop_retained_of_memento_geneset, fill = geneset)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", prop_retained_of_memento_geneset)),
            position = position_dodge(width = 0.7), vjust = -0.4, size = 3) +
  labs(
    title = "Retention of Memento's consistently low/high dispersion genes",
    subtitle = "Proportion of the manuscript's canonical gene set (quintile 1/5 in all 6 cell types)\nalso classified the same way under an alternative dispersion-estimation method",
    x = NULL, y = "Proportion of Memento's gene set retained", fill = NULL
  ) +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
# PDF (not PNG) per Brendan's request, 2026-09-03 -- vector output for
# editing in Affinity Designer; cairo_pdf matches the repo's device
# convention for vector figures elsewhere.
ggsave(file.path(FIG_DIR, "retention_by_method.pdf"), p, width = 7, height = 5.5, device = cairo_pdf)

# ---- reportable values (R-side register_value shim, per dqtl_batch_structure
# precedent -- register_value.py infers namespace from the caller's Python
# file path, so an R script can't call it directly) --------
reportable <- do.call(rbind, lapply(seq_len(nrow(summary_df)), function(i) {
  row <- summary_df[i, ]
  key_prefix <- paste0(row$method, "_", row$geneset)
  data.frame(
    key = c(
      paste0("prop_retained_of_memento_geneset_", key_prefix),
      paste0("n_retained_consistent_", key_prefix)
    ),
    value = c(signif(row$prop_retained_of_memento_geneset, 4), row$n_retained_consistent),
    provenance = "figure_reproduction/scripts/basics_scran_crossmethod/10_consistent_geneset_retention.R",
    stringsAsFactors = FALSE
  )
}))
write.csv(reportable, file.path(OUT_DIR, "reportable_values_06.csv"), row.names = FALSE)
message("Reportable values written to ", file.path(OUT_DIR, "reportable_values_06.csv"))
