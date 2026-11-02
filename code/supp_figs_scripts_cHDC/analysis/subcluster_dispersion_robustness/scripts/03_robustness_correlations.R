#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(ComplexHeatmap)
  library(circlize)
})

# ---- Named constants (lowercase per R convention) ---------------------------
repo_dir     <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
trial_dir    <- file.path(repo_dir, "cellranger_cluster-mode_trial")
analysis_dir <- file.path(repo_dir, "analysis/subcluster_dispersion_robustness")

# Which step-02 config to read from / write results into. Must match a
# directory step 02 actually wrote (e.g. "mincells1000_minind2" or
# "mincells700_minind2" -- see 02_memento_per_subcluster.py's CONFIG_LABEL).
config_label <- Sys.getenv("SUBCLUSTER_CONFIG_LABEL", unset = "mincells1000_minind2")
out_dir      <- file.path(analysis_dir, "outputs", config_label)
if (!dir.exists(out_dir)) {
  stop(sprintf(
    "Config directory does not exist: %s. Run 02_memento_per_subcluster.py with SUBCLUSTER_MIN_CELLS_PER_GROUP/SUBCLUSTER_MIN_INDIVIDUALS set to produce it, or set SUBCLUSTER_CONFIG_LABEL to an existing config.",
    out_dir
  ))
}

# Original cell-type-level dispersion (residual variance), genes x 18 columns
# named {individual}_{cell_type}.
original_disp_csv <- file.path(
  trial_dir,
  "analysis/datasets/memento/new_res_var_by_cell_type_ind_1000_cts_new_predicted_id.csv"
)
# Parallel mean-expression matrix, same shape/columns -- source of the -6
# mean-threshold filter below.
original_mean_csv <- file.path(
  trial_dir,
  "analysis/datasets/memento/new_mean_expression_by_cell_type_ind_1000_cts_new_predicted_id.csv"
)

# Resolution -> per-sub-cluster residual-variance / mean-expression CSVs from
# step 02 (same directory, same column set per resolution).
resolution_files <- list(
  "0.5" = file.path(out_dir, "res_var_by_subcluster_res_0.5.csv"),
  "1.0" = file.path(out_dir, "res_var_by_subcluster_res_1.0.csv"),
  "2.0" = file.path(out_dir, "res_var_by_subcluster_res_2.0.csv")
)
mean_resolution_files <- list(
  "0.5" = file.path(out_dir, "mean_expression_by_subcluster_res_0.5.csv"),
  "1.0" = file.path(out_dir, "mean_expression_by_subcluster_res_1.0.csv"),
  "2.0" = file.path(out_dir, "mean_expression_by_subcluster_res_2.0.csv")
)

# Ribosomal/mitochondrial gene-symbol prefixes and Alexander's Memento
# mean-expression exclusion threshold -- repo convention established in
# dd_excl_topexpressed/scripts/01_exclude_topexpressed_genes.R (ribo, -6
# any-individual rule) and basics_dispersion_reanalysis/scripts/
# 01_prepare_basics_inputs.R (adds mito). Wired into this analysis 2026-08-11
# per todo/wire-neg6-exclusion-remaining-revisions.md.
ribo_pattern        <- "^RP"
mito_pattern        <- "^MT-"
mean_excl_threshold <- -6

cell_types <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)
# ANALYSIS_OK[fixed-design]: the 3 pooled cHDC iPSC lines (the entire study
# cohort), not an exclusion filter; expected columns asserted in the orig_ct
# aggregation below.
individuals <- c("NA19099", "NA19093", "NA18522")

n_quintiles      <- 5L
cor_method       <- "spearman"
# Minimum shared-gene overlap fraction between original and sub-cluster
# estimates below which we assume a gene-ID mismatch and abort.
min_overlap_frac <- 0.3
n_permutations   <- 1000L
seed             <- 1L

# Across-sub-cluster correlation heatmap color scale (Brendan, 2026-08-11):
# fixed to this range rather than the full [-1, 1] because observed rho across
# all resolutions/configs falls in this band (see F-021,
# .living/findings/dispersion-robustness.md) -- a full [-1,1] scale would wash
# out the real within-vs-across-cell-type structure.
heatmap_rho_range <- c(0.8, 1)

set.seed(seed)

# Logging logic lives in scripts/utils.R (shared with 01) to avoid duplicate
# definitions drifting apart (scilintr R025).
source(file.path(repo_dir, "analysis/subcluster_dispersion_robustness/scripts/utils.R"))
log_msg <- make_logger("03_robustness")
log_msg("Config: %s (reading/writing %s)", config_label, out_dir)

# ---- Helpers ----------------------------------------------------------------

# Load a Memento residual-variance CSV: first column (or a "gene" column) is the
# gene id, remaining columns are per-group estimates.
# ANALYSIS_OK[within-file-helper]: called in the per-resolution loop below; R026's
# cross-file rule cannot see intra-file callers.
load_memento_csv <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Missing Memento output: %s. Run 02_memento_per_subcluster.py.", path))
  }
  df <- read.csv(path, row.names = 1, check.names = FALSE)
  if ("gene" %in% colnames(df)) {
    rownames(df) <- df$gene
    df$gene <- NULL
  }
  # ANALYSIS_OK[complete-cases]: keep only genes with a finite Memento estimate
  # in every group (required for pairwise Spearman across sub-clusters); the
  # retained gene count is logged by the caller ("%d genes x %d groups") below.
  df <- df[stats::complete.cases(df), , drop = FALSE]
  as.matrix(df)
}

# Parse Memento group-column names into (individual, subcluster, parent_ct).
# We do NOT rely on Memento's exact separator: individuals match ^NA[0-9]+ and
# sub-cluster labels are "{cell_type}::{cluster}" (set in step 01). Fail loudly
# on any column we cannot parse.
# ANALYSIS_OK[within-file-helper]: called in the per-resolution loop below; R026's
# cross-file rule cannot see intra-file callers.
parse_group_columns <- function(cols) {
  ind   <- regmatches(cols, regexpr("NA[0-9]+", cols))
  sub_m <- regmatches(cols, regexpr("[A-Za-z_]+::[0-9]+", cols))
  if (length(ind) != length(cols) || length(sub_m) != length(cols)) {
    bad <- cols[!grepl("NA[0-9]+", cols) | !grepl("[A-Za-z_]+::[0-9]+", cols)]
    stop(sprintf("Could not parse individual/sub-cluster from columns: %s",
                 paste(head(bad, 5), collapse = ", ")))
  }
  parent <- sub("::.*$", "", sub_m)
  if (!all(parent %in% cell_types)) {
    stop(sprintf("Parsed parent cell types not in expected set: %s",
                 paste(setdiff(unique(parent), cell_types), collapse = ", ")))
  }
  data.frame(col = cols, individual = ind, subcluster = sub_m,
             parent_ct = parent, stringsAsFactors = FALSE)
}

# Collapse a genes x group matrix to genes x sub-cluster by averaging across
# individuals (one dispersion vector per sub-cluster).
# TODO(assumption): aggregating per-individual estimates by simple mean before
# correlation. An alternative keeps individual-level columns; confirm.
# ANALYSIS_OK[within-file-helper]: called in the per-resolution loop below; R026's
# cross-file rule cannot see intra-file callers.
aggregate_by_subcluster <- function(mat, meta) {
  subclusters <- unique(meta$subcluster)
  agg <- sapply(subclusters, function(sc) {
    cols <- meta$col[meta$subcluster == sc]
    rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)
  })
  colnames(agg) <- subclusters
  agg
}

# Drop ^RP ribosomal and ^MT- mitochondrial genes (manuscript convention,
# main_figs_1-4_related_supp.Rmd lines 52-60 for ribo; mito added per
# basics_dispersion_reanalysis's 01_prepare_basics_inputs.R). Returns a logical
# keep-mask aligned to `genes`.
# ANALYSIS_OK[within-file-helper]: called for both the original cell-type-level
# matrices and each resolution's sub-cluster matrices below; R026's cross-file
# rule cannot see intra-file callers.
keep_nonribo_nonmito <- function(genes) {
  !(grepl(ribo_pattern, genes) | grepl(mito_pattern, genes))
}

# Alexander's -6 Memento-mean-expression exclusion, any-individual rule
# (dd_excl_topexpressed/scripts/01_exclude_topexpressed_genes.R): a gene is
# excluded from a group (cell type or sub-cluster) if ANY of that group's
# `cols` mean values exceeds `threshold`. Returns the excluded gene names.
# ANALYSIS_OK[within-file-helper]: called once per cell type (original data)
# and once per sub-cluster (each resolution) below; R026's cross-file rule
# cannot see intra-file callers.
mean_exclude_group <- function(mean_mat, cols, threshold) {
  present_cols <- cols[cols %in% colnames(mean_mat)]
  if (length(present_cols) == 0) {
    stop(sprintf("None of the expected mean-matrix columns are present: %s",
                 paste(cols, collapse = ", ")))
  }
  sub_mean <- mean_mat[, present_cols, drop = FALSE]
  rownames(sub_mean)[apply(sub_mean, 1, function(x) any(x > threshold, na.rm = TRUE))]
}

# ---- Load original cell-type-level dispersion, aggregate across individuals --
log_msg("Loading original cell-type-level dispersion: %s", original_disp_csv)
orig <- read.csv(original_disp_csv, row.names = 1, check.names = FALSE)
# ANALYSIS_OK[complete-cases]: drop genes without a finite estimate in all 18
# individual-x-cell-type columns; retained gene count logged on the next line.
orig <- as.matrix(orig[stats::complete.cases(orig), , drop = FALSE])
log_msg("Original matrix: %d genes x %d columns.", nrow(orig), ncol(orig))

log_msg("Loading original cell-type-level mean expression: %s", original_mean_csv)
orig_mean <- read.csv(original_mean_csv, row.names = 1, check.names = FALSE)

# Align dispersion and mean matrices to their shared gene universe before
# filtering (the two CSVs are independently complete-cases-filterable and can
# diverge slightly), then drop ribosomal/mitochondrial genes from both.
common_orig <- intersect(rownames(orig), rownames(orig_mean))
log_msg("Dispersion/mean gene overlap: %d shared (dispersion %d, mean %d).",
        length(common_orig), nrow(orig), nrow(orig_mean))
# ANALYSIS_OK[align-to-shared-ids]: restrict both matrices to their shared gene
# universe; dropped-row counts logged on the log_msg call immediately above.
orig      <- orig[common_orig, , drop = FALSE]
# ANALYSIS_OK[align-to-shared-ids]: as above -- same narrowing, same log line.
orig_mean <- orig_mean[common_orig, , drop = FALSE]

keep_orig <- keep_nonribo_nonmito(rownames(orig))
log_msg("Dropping %d ribosomal/mitochondrial genes from original data (%d retained).",
        sum(!keep_orig), sum(keep_orig))
# ANALYSIS_OK[filtering]: drop ^RP ribosomal / ^MT- mitochondrial genes (repo
# convention, see keep_nonribo_nonmito()); dropped count logged immediately above.
orig      <- orig[keep_orig, , drop = FALSE]
# ANALYSIS_OK[filtering]: as above -- same narrowing, same log line.
orig_mean <- orig_mean[keep_orig, , drop = FALSE]

# Aggregate original to genes x cell_type by averaging the 3 individuals, with
# Alexander's -6 mean-threshold exclusion (any-individual rule) applied per
# cell type: a gene excluded from cell type X is NA'd out of X's column only
# (not dropped from other cell types' columns here) before averaging.
n_mean_excluded_by_ct <- integer(length(cell_types))
names(n_mean_excluded_by_ct) <- cell_types
orig_ct <- sapply(cell_types, function(ct) {
  cols <- paste(individuals, ct, sep = "_")
  missing <- setdiff(cols, colnames(orig))
  if (length(missing) > 0) {
    stop(sprintf("Original CSV missing expected columns: %s",
                 paste(missing, collapse = ", ")))
  }
  disp_ct <- orig[, cols, drop = FALSE]
  excluded_ct <- mean_exclude_group(orig_mean, cols, mean_excl_threshold)
  n_mean_excluded_by_ct[[ct]] <<- length(excluded_ct)
  vals <- rowMeans(disp_ct, na.rm = TRUE)
  vals[rownames(disp_ct) %in% excluded_ct] <- NA_real_
  vals
})
colnames(orig_ct) <- cell_types
for (ct in cell_types) {
  log_msg("  [%s] mean-threshold excl (Memento mean > %s in >=1 individual): %d genes.",
          ct, mean_excl_threshold, n_mean_excluded_by_ct[[ct]])
}

# Restrict to genes retained (ribo/mito-clean and below the mean threshold) in
# ALL 6 cell types -- "union-retained", matching dd_excl_topexpressed's Part B
# convention for analyses that need one shared gene universe across cell types.
n_before_union <- nrow(orig_ct)
orig_ct <- orig_ct[stats::complete.cases(orig_ct), , drop = FALSE]
log_msg("Aggregated original to %d genes x %d cell types (%d dropped: excluded in >=1 cell type).",
        nrow(orig_ct), ncol(orig_ct), n_before_union - nrow(orig_ct))

# ---- Derive consistently-low / consistently-high gene sets from original -----
# Quintile-rank genes within each cell type; low = bottom quintile in ALL cell
# types, high = top quintile in ALL cell types.
# ANALYSIS_OK[within-file-helper]: called for orig_ct and each sub-cluster below;
# R026's cross-file rule cannot see intra-file callers.
quintile_rank <- function(v) {
  # na.last = "keep" leaves NA (mean-threshold-excluded genes) as NA rather
  # than assigning them a rank -- otherwise rank()'s default na.last = TRUE
  # would push them to the top quintile.
  r <- rank(v, ties.method = "average", na.last = "keep")
  as.integer(cut(r,
                 breaks = quantile(r, probs = seq(0, 1, length.out = n_quintiles + 1),
                                   na.rm = TRUE),
                 include.lowest = TRUE, labels = FALSE))
}
orig_quint <- apply(orig_ct, 2, quintile_rank)
rownames(orig_quint) <- rownames(orig_ct)

consistently_low  <- rownames(orig_quint)[apply(orig_quint, 1, function(x) all(x == 1L))]
consistently_high <- rownames(orig_quint)[apply(orig_quint, 1, function(x) all(x == n_quintiles))]
log_msg("Gene sets: %d consistently-low, %d consistently-high (of %d genes).",
        length(consistently_low), length(consistently_high), nrow(orig_quint))
if (length(consistently_low) == 0 || length(consistently_high) == 0) {
  stop("Empty consistently-low/high gene set derived from original data.")
}

# ---- Per-resolution analysis ------------------------------------------------
vs_original_rows   <- list()
across_summary_rows <- list()
geneset_rows       <- list()

for (res_label in names(resolution_files)) {
  log_msg("=== Resolution %s ===", res_label)
  mat      <- load_memento_csv(resolution_files[[res_label]])
  mean_mat <- load_memento_csv(mean_resolution_files[[res_label]])

  # Align dispersion and mean matrices (independently complete-cases-filtered,
  # can diverge slightly) before ribo/mito drop, same pattern as the original
  # cell-type-level data above.
  common_sc_genes <- intersect(rownames(mat), rownames(mean_mat))
  log_msg("  Dispersion/mean gene overlap: %d shared (dispersion %d, mean %d).",
          length(common_sc_genes), nrow(mat), nrow(mean_mat))
  # ANALYSIS_OK[align-to-shared-ids]: restrict both matrices to their shared
  # gene universe; dropped-row counts logged on the log_msg call immediately above.
  mat      <- mat[common_sc_genes, , drop = FALSE]
  # ANALYSIS_OK[align-to-shared-ids]: as above -- same narrowing, same log line.
  mean_mat <- mean_mat[common_sc_genes, , drop = FALSE]

  keep_sc <- keep_nonribo_nonmito(rownames(mat))
  log_msg("  Dropping %d ribosomal/mitochondrial genes (%d retained).",
          sum(!keep_sc), sum(keep_sc))
  # ANALYSIS_OK[filtering]: drop ^RP ribosomal / ^MT- mitochondrial genes (repo
  # convention, see keep_nonribo_nonmito()); dropped count logged immediately above.
  mat      <- mat[keep_sc, , drop = FALSE]
  # ANALYSIS_OK[filtering]: as above -- same narrowing, same log line.
  mean_mat <- mean_mat[keep_sc, , drop = FALSE]

  meta <- parse_group_columns(colnames(mat))
  log_msg("  %d genes x %d groups; %d sub-clusters across %d parent cell types.",
          nrow(mat), ncol(mat), length(unique(meta$subcluster)),
          length(unique(meta$parent_ct)))

  # Alexander's -6 mean-threshold exclusion, any-individual rule, applied per
  # sub-cluster (analogous to the per-cell-type application above): NA out a
  # sub-cluster's entry for any gene where >=1 of its individuals in that
  # sub-cluster has a Memento mean > mean_excl_threshold.
  n_mean_excluded_total <- 0L
  for (sc in unique(meta$subcluster)) {
    sc_cols <- meta$col[meta$subcluster == sc]
    excluded_sc <- mean_exclude_group(mean_mat, sc_cols, mean_excl_threshold)
    if (length(excluded_sc) > 0) {
      mat[rownames(mat) %in% excluded_sc, sc_cols] <- NA_real_
      n_mean_excluded_total <- n_mean_excluded_total + length(excluded_sc)
    }
  }
  log_msg("  Mean-threshold excl (Memento mean > %s in >=1 individual): %d (sub-cluster, gene) entries NA'd across %d sub-clusters.",
          mean_excl_threshold, n_mean_excluded_total, length(unique(meta$subcluster)))

  # ID overlap check vs original (fail loudly on mismatch).
  common <- intersect(rownames(mat), rownames(orig_ct))
  overlap_frac <- length(common) / min(nrow(mat), nrow(orig_ct))
  log_msg("  Shared genes with original: %d (overlap fraction %.2f).",
          length(common), overlap_frac)
  if (overlap_frac < min_overlap_frac) {
    stop(sprintf(paste0("Gene-ID overlap between sub-cluster estimates and the ",
                        "original (%.2f) is below %.2f. Likely an identifier ",
                        "mismatch (ENSG vs symbol). Harmonize IDs before ",
                        "correlating."), overlap_frac, min_overlap_frac))
  }

  sub_agg <- aggregate_by_subcluster(mat, meta)
  # ANALYSIS_OK[align-to-shared-ids]: restrict to genes shared with the original
  # (common); overlap count/fraction is logged just above and the run aborts if
  # it falls below min_overlap_frac.
  sub_agg <- sub_agg[common, , drop = FALSE]
  parent_of <- setNames(sub("::.*$", "", colnames(sub_agg)), colnames(sub_agg))

  # ---- (a) across-sub-cluster correlation within this resolution -----------
  cor_mat <- cor(sub_agg, method = cor_method, use = "pairwise.complete.obs")
  write.csv(cor_mat,
            file.path(out_dir, sprintf("across_subcluster_cor_res_%s.csv", res_label)))

  cor_long <- reshape2::melt(cor_mat, varnames = c("subcluster_a", "subcluster_b"),
                             value.name = "rho")
  # ANALYSIS_OK[drop-self-pairs]: removes the diagonal (self-correlation, rho=1)
  # entries of the melted correlation matrix, not analysis genes/samples.
  cor_long <- cor_long[cor_long$subcluster_a != cor_long$subcluster_b, ]
  cor_long$same_parent <- parent_of[as.character(cor_long$subcluster_a)] ==
                          parent_of[as.character(cor_long$subcluster_b)]
  across_summary_rows[[res_label]] <- data.frame(
    resolution = res_label,
    n_subclusters = ncol(sub_agg),
    mean_rho_within_celltype  = mean(cor_long$rho[cor_long$same_parent]),
    mean_rho_across_celltype  = mean(cor_long$rho[!cor_long$same_parent]),
    stringsAsFactors = FALSE
  )

  # Heatmap via ComplexHeatmap (.living/conventions.md, 2026-08-11 -- all
  # heatmaps in this project use ComplexHeatmap, not ggplot geom_tile()).
  # Color scale fixed to heatmap_rho_range (0.8-1), not the full [-1,1].
  hm_col_fun <- circlize::colorRamp2(
    c(heatmap_rho_range[1], mean(heatmap_rho_range), heatmap_rho_range[2]),
    c("blue", "white", "red")
  )
  hm <- ComplexHeatmap::Heatmap(
    cor_mat,
    name = "Spearman rho",
    col = hm_col_fun,
    column_title = sprintf("Across-sub-cluster dispersion correlation (res %s)", res_label),
    row_names_gp = grid::gpar(fontsize = 6),
    column_names_gp = grid::gpar(fontsize = 6)
  )
  png(file.path(out_dir, sprintf("across_subcluster_cor_heatmap_res_%s.png", res_label)),
      width = 10, height = 9, units = "in", res = 200)
  ComplexHeatmap::draw(hm)
  dev.off()

  # ---- (b) each sub-cluster vs original parent cell-type estimate -----------
  for (sc in colnames(sub_agg)) {
    parent <- parent_of[[sc]]
    x <- sub_agg[, sc]
    y <- orig_ct[common, parent]
    ok <- is.finite(x) & is.finite(y)
    rho <- cor(x[ok], y[ok], method = cor_method)
    vs_original_rows[[length(vs_original_rows) + 1]] <- data.frame(
      resolution = res_label, subcluster = sc, parent_ct = parent,
      n_genes = sum(ok), rho = rho, stringsAsFactors = FALSE
    )
  }

  vs_df_res <- do.call(rbind, vs_original_rows[
    sapply(vs_original_rows, function(d) d$resolution[1] == res_label)])
  sp <- ggplot(vs_df_res, aes(parent_ct, rho)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, height = 0, size = 1.5, alpha = 0.7) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = sprintf("Sub-cluster vs original cell-type dispersion (res %s)", res_label),
         x = "Parent cell type", y = "Spearman rho (sub-cluster vs original)")
  ggsave(file.path(out_dir, sprintf("subcluster_vs_original_rho_res_%s.png", res_label)),
         sp, width = 9, height = 6, dpi = 200)

  # ---- (c) gene-set robustness ---------------------------------------------
  # Quintile-rank genes within each sub-cluster; check where the original
  # consistently-low/high genes land.
  low_in_res  <- intersect(consistently_low,  rownames(sub_agg))
  high_in_res <- intersect(consistently_high, rownames(sub_agg))

  for (sc in colnames(sub_agg)) {
    q <- quintile_rank(sub_agg[, sc])
    names(q) <- rownames(sub_agg)
    # na.rm = TRUE: q is NA for genes the mean-threshold exclusion NA'd out of
    # this sub-cluster (see the per-sub-cluster masking above); those genes are
    # excluded from the retention fraction, not counted as failures.
    frac_low_bottom  <- mean(q[low_in_res]  == 1L, na.rm = TRUE)
    frac_high_top    <- mean(q[high_in_res] == n_quintiles, na.rm = TRUE)

    # Permutation null: random gene sets of the same sizes.
    # ANALYSIS_OK[permutation-null]: sample() draws random gene subsets to build
    # the null distribution of quintile retention, not synthetic data; the RNG is
    # seeded once at top level (set.seed(seed) near the top of the script).
    perm_low  <- replicate(n_permutations,
                           mean(sample(q, length(low_in_res))  == 1L, na.rm = TRUE))
    # ANALYSIS_OK[permutation-null]: as above -- null-distribution draw, top-level seed.
    perm_high <- replicate(n_permutations,
                           mean(sample(q, length(high_in_res)) == n_quintiles, na.rm = TRUE))
    p_low  <- (sum(perm_low  >= frac_low_bottom) + 1) / (n_permutations + 1)
    p_high <- (sum(perm_high >= frac_high_top)   + 1) / (n_permutations + 1)

    geneset_rows[[length(geneset_rows) + 1]] <- data.frame(
      resolution = res_label, subcluster = sc, parent_ct = parent_of[[sc]],
      n_low = length(low_in_res), n_high = length(high_in_res),
      frac_low_in_bottom_quintile = frac_low_bottom,
      frac_high_in_top_quintile   = frac_high_top,
      perm_p_low  = p_low,
      perm_p_high = p_high,
      stringsAsFactors = FALSE
    )
  }
}

# ---- Write summary tables ---------------------------------------------------
across_summary <- do.call(rbind, across_summary_rows)
write.csv(across_summary,
          file.path(out_dir, "across_subcluster_correlation_summary.csv"),
          row.names = FALSE)
log_msg("Wrote across-sub-cluster correlation summary.")

vs_original <- do.call(rbind, vs_original_rows)
write.csv(vs_original,
          file.path(out_dir, "subcluster_vs_original_correlation.csv"),
          row.names = FALSE)
log_msg("Wrote sub-cluster vs original correlation table (%d rows).", nrow(vs_original))

geneset_robustness <- do.call(rbind, geneset_rows)
write.csv(geneset_robustness,
          file.path(out_dir, "geneset_robustness_summary.csv"),
          row.names = FALSE)
log_msg("Wrote gene-set robustness summary (%d rows).", nrow(geneset_robustness))

# Boxplot: gene-set quintile-retention fractions across sub-clusters, per res.
gl <- reshape2::melt(
  geneset_robustness,
  id.vars = c("resolution", "subcluster", "parent_ct"),
  measure.vars = c("frac_low_in_bottom_quintile", "frac_high_in_top_quintile"),
  variable.name = "gene_set", value.name = "fraction"
)
# Legend labels renamed from the raw column names to reader-facing text per
# Brendan's request, 2026-09-02 -- named-vector lookup (not dplyr::recode)
# since this script has no dplyr dependency elsewhere.
gene_set_labels <- c(frac_low_in_bottom_quintile = "low-dispersion genes",
                      frac_high_in_top_quintile = "high-dispersion genes")
# ANALYSIS_OK[relabeling]: display-only label remap for the legend, not a row
# filter -- no rows dropped or reordered; asserted below.
gl$gene_set <- gene_set_labels[as.character(gl$gene_set)]
stopifnot("Every gene_set value must have a display label -- unmapped value(s) would silently become NA" =
            !anyNA(gl$gene_set))
gl$gene_set <- factor(gl$gene_set, levels = gene_set_labels)
gp <- ggplot(gl, aes(resolution, fraction, fill = gene_set)) +
  geom_boxplot(position = position_dodge(0.8)) +
  theme_bw() +
  labs(title = "Gene-set quintile retention across sub-clusters",
       x = "Clustering resolution",
       y = "Fraction of set retained in extreme quintile",
       fill = NULL)
# PDF (not PNG) per Brendan's stated preference for editable vector figures
# (2026-09-02); cairo_pdf matches the repo's device convention.
ggsave(file.path(out_dir, "geneset_robustness_boxplot.pdf"), gp,
       width = 8, height = 6, device = cairo_pdf)

log_msg("DONE. Outputs written to %s", out_dir)
