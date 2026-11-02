# ---- Paths ------------------------------------------------------------------
repo_dir     <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
src_dir      <- file.path(repo_dir, "cellranger_cluster-mode_trial")
analysis_dir <- file.path(repo_dir, "analysis/mean_matched_regulatory_revision")
# PATCHED: this fork's own output tree, not the original analysis's.
out_dir      <- file.path(repo_dir, "figure_reproduction/outputs/mean_matched_regulatory_revision")

# Input matrices (genes x 18; columns {ind}_{cell_type}).
# PATCHED: Epithelial_cells sub-cluster-24-excluded canonical matrices (only
# Epithelial_cells' 3 columns differ from the published matrix; see
# figure_reproduction/FIGURE_REPRODUCTION.md).
dispersion_matrix <- file.path(repo_dir, "figure_reproduction/data/res_var_excl_cluster24.csv")
mean_matrix       <- file.path(repo_dir, "figure_reproduction/data/mean_expression_excl_cluster24.csv")
# Per-cell-type -6 top-expressed-gene exclusion indicator (1 = retained), from
# figure_reproduction/scripts/01_exclude_topexpressed_genes.R.
retained_matrix   <- file.path(repo_dir, "figure_reproduction/outputs/retained_gene_matrix.csv")

# s_het per gene (Pritchard). UNCHANGED -- not read by 01/02/05/07 (only by 03,
# not copied here).
shet_xlsx <- file.path(
  src_dir, "analysis/datasets/external_data",
  "Pritchard_shet_supp_data.xlsx"
)

# eQTL results (Erik McIntire). UNCHANGED -- not read by 01/02/05/07 (only 03).
eqtl_result_template <- paste0(
  "/project/gilad/emcintire/chemo_project_2023/pipeline/map_eqtls/",
  "archive/no_cvrt_condox_only/eqtl_results/Control/",
  "%s/steadyState_%s_Control_AllChroms.txt"
)

# DE / DD contrast tables. UNCHANGED -- not read by 01/02/05/07 (only 03).
de_contrast_template <- file.path(
  src_dir, "analysis/first_pub_datasets",
  "%s_mean_var_new_pred_id_contrasts.csv"
)

# Extension (a) housekeeping list path. UNCHANGED -- not read by 01/02/05/07
# (only 04, not copied here).
hk_gene_list_path <- file.path(analysis_dir, "inputs/housekeeping_genes.txt")

# Optional biomaRt cache. UNCHANGED -- not read by 01/02/05/07 (only 03/04).
exon_counts_cache <- file.path(analysis_dir, "outputs", "regulatory", "exon_counts_mane.rds")
bm_mane_cache     <- file.path(analysis_dir, "outputs", "regulatory", "bm_mane_unique.rds")

# ---- Biology ----------------------------------------------------------------
cell_types <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)
individuals  <- c("NA19099", "NA19093", "NA18522")
ribo_pattern <- "^RP"

# ---- Thresholds (named constants) -------------------------------------------
n_quantiles            <- 5
n_mean_percentile_bins <- 100
quintile_membership_min <- length(cell_types)
lowest_quintile  <- 1L
highest_quintile <- n_quantiles

egene_padj      <- 0.05
egene_min_count <- 49

match_seed_primary <- 1L
match_seeds        <- c(1L, 2L, 3L)

network_retry_attempts <- 6L
network_retry_wait_sec <- 10

lfc_not_de <- 1.0
dd_adjp <- 0.05

# ---- Small robust-analysis helpers -------------------------------------------
log_step <- function(msg, n = NULL) {
  if (is.null(n)) {
    message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg))
  } else {
    message(sprintf("[%s] %s  (n = %d)", format(Sys.time(), "%H:%M:%S"), msg, n))
  }
}

assert_that <- function(condition, msg) {
  if (!isTRUE(condition)) stop(sprintf("ASSERTION FAILED: %s", msg), call. = FALSE)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

# ---- Shared statistics / labelling helpers (used by 05) ---------------------
wilcox_effect <- function(values, group, level_a, level_b, feature) {
  x <- values[group == level_a]; y <- values[group == level_b]
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  assert_that(length(x) > 0 && length(y) > 0,
              sprintf("empty group in %s: %s vs %s", feature, level_a, level_b))
  wt <- wilcox.test(x, y, conf.int = TRUE)
  data.frame(
    feature = feature, comparison = sprintf("%s_vs_%s", level_a, level_b),
    n_a = length(x), n_b = length(y),
    median_a = median(x), median_b = median(y),
    hl_shift = unname(wt$estimate),
    ci_lower = wt$conf.int[1], ci_upper = wt$conf.int[2],
    p_value = wt$p.value, row.names = NULL)
}

format_pvalue_label <- function(p) {
  ifelse(p < 2e-16, "p < 2e-16", paste0("p = ", formatC(p, format = "e", digits = 2)))
}

label_group <- function(genes_list) {
  dplyr::bind_rows(
    data.frame(gene = genes_list$first, rank_category = "Low-dispersion"),
    data.frame(gene = genes_list$fifth, rank_category = "High-dispersion"),
    data.frame(gene = genes_list$never, rank_category = "Neither")
  )
}
