#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(clusterProfiler)  # enricher, read.gmt, dotplot via enrichplot
  library(enrichplot)
})

repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
dir <- file.path(repo_dir, "cellranger_cluster-mode_trial")
memento_dir <- file.path(dir, "analysis", "datasets", "memento")
# ANALYSIS_OK[external-path]: read-only MSigDB C5 GO:BP GMT used by every GO analysis here; existence asserted below
gmt_file <- file.path(dir, "analysis", "datasets", "gsea_gene_sets", "background_genes",
                      "c5.go.bp.v2024.1.Hs.symbols.gmt")
pub_res_dir <- file.path(dir, "analysis", "first_pub_datasets", "pathway_analysis", "alt_cluster_results")
out_dir <- file.path(repo_dir, "figure_reproduction", "outputs", "figS_alt_cluster_go")

cluster_resolutions <- c("res_05", "res_08")
ribo_pattern <- "^RP"
n_quantiles <- 5
min_low_celltypes <- 6     # source: rowSums(. == 1) >= 6
padj_cutoff <- 0.05        # source
genecount_cutoff <- 5      # source: Count > 5
n_show <- 8                # source: showCategory = 8
match_tol <- 1e-10
mean_excl_threshold <- -6  # script 01: exclude if any individual's mean > -6

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(gmt_file)) stop("Missing GMT: ", gmt_file)
gmt <- read.gmt(gmt_file)

# Consistently-low genes and GO background for one resolution. With
# apply_filter = FALSE this is the source code path; with TRUE, the -6 rule.
low_genes_and_universe <- function(cluster_res, apply_filter) {
  in_csv <- file.path(memento_dir, sprintf("new_res_var_by_cell_type_&_ind_subset_8_cts_%s.csv", cluster_res))
  mean_csv <- file.path(memento_dir, sprintf("new_mean_expression_by_cell_type_&_ind_subset_8_%s.csv", cluster_res))
  for (f in c(in_csv, mean_csv)) if (!file.exists(f)) stop("Missing input: ", f)
  res_var_df <- read.csv(in_csv, row.names = 1)
  # ANALYSIS_OK[filtering]: ^RP removal as in the source (ribo_genes_full); count logged below
  res_var_df <- res_var_df[!grepl(ribo_pattern, rownames(res_var_df)), ]
  mean_df <- read.csv(mean_csv, row.names = 1)
  if (!all(rownames(res_var_df) %in% rownames(mean_df)) || !setequal(colnames(mean_df), colnames(res_var_df))) {
    stop("Mean matrix ", basename(mean_csv), " does not match the dispersion matrix's genes/samples")
  }
  if (anyNA(res_var_df) || anyNA(mean_df[rownames(res_var_df), ])) stop("NA values in ", cluster_res, " matrices")

  split_parts <- strsplit(colnames(res_var_df), "_")
  cell_types <- unique(sapply(split_parts, function(x) paste(x[-1], collapse = "_")))
  genes <- rownames(res_var_df)

  # Per cell type: retained genes, median dispersion over the 3 individuals,
  # rank and quintile over the retained genes (NA for excluded genes).
  quantile_bins <- matrix(NA_integer_, nrow = length(genes), ncol = length(cell_types),
                          dimnames = list(genes, cell_types))
  retained <- matrix(TRUE, nrow = length(genes), ncol = length(cell_types), dimnames = list(genes, cell_types))
  for (ct in cell_types) {
    cols <- paste0(c("NA19099_", "NA19093_", "NA18522_"), ct)
    if (!all(cols %in% colnames(res_var_df))) stop("Missing columns for ", ct, " in ", in_csv)
    if (apply_filter) {
      retained[, ct] <- !apply(mean_df[genes, cols] > mean_excl_threshold, 1, any)
    }
    keep <- genes[retained[, ct]]
    ranks <- rank(apply(res_var_df[keep, cols], 1, median), ties.method = "average")
    quantile_bins[keep, ct] <- cut(ranks, breaks = quantile(ranks, probs = seq(0, 1, length.out = n_quantiles + 1)),
                                   labels = FALSE, include.lowest = TRUE)
  }
  low_genes <- genes[rowSums(quantile_bins == 1, na.rm = TRUE) >= min_low_celltypes]
  universe <- if (apply_filter) genes[rowSums(retained) == length(cell_types)] else genes
  if (!all(low_genes %in% universe)) stop("Consistently-low genes outside the GO universe for ", cluster_res)
  message(sprintf("[47] %s (%s): %d non-ribo genes, %d cell types, %d excluded by -6 in >= 1 cell type, universe %d, %d consistently-low genes",
                  cluster_res, if (apply_filter) "-6 filter" else "no filter", length(genes), length(cell_types),
                  sum(rowSums(!retained) > 0), length(universe), length(low_genes)))
  list(low_genes = low_genes, universe = universe)
}

# enricher() + the source's post-filter and Description cleanup.
run_go <- function(low_genes, universe) {
  bg_genes <- gmt[gmt$gene %in% universe, ]
  res <- enricher(low_genes, TERM2GENE = bg_genes)
  if (is.null(res)) stop("No enrichment result")
  res_df <- res@result %>% mutate(minuslog10padj = -log10(p.adjust),
                                  pathway_name = gsub("\\.GOBP.*$", "", rownames(res@result)))
  target_pws <- unique(res_df$ID[res_df$p.adjust < padj_cutoff & res_df$Count > genecount_cutoff])
  res_df <- res_df[res_df$ID %in% target_pws, ]
  rownames(res_df) <- res_df$ID
  res_df$Description <- gsub("GOBP_", "", res_df$Description)
  res_df$Description <- str_wrap(res_df$Description, width = 10)
  list(res_df = res_df, bg_genes = bg_genes)
}

for (cluster_res in cluster_resolutions) {
  # Check: unfiltered run reproduces the published result table.
  chk <- low_genes_and_universe(cluster_res, apply_filter = FALSE)
  chk_df <- run_go(chk$low_genes, chk$universe)$res_df
  pub_csv <- file.path(pub_res_dir, sprintf("%s_lowest_quintile_ranked_dispersion_%s_GO_quintile_ranked_bp_res_df.csv",
                                            cluster_res, cluster_res))
  if (!file.exists(pub_csv)) stop("Missing published result: ", pub_csv)
  pub <- read.csv(pub_csv)
  if (!identical(sort(chk_df$ID), sort(pub$ID)) ||
      !(max(abs(chk_df[pub$ID, "p.adjust"] - pub$p.adjust)) < match_tol)) {
    stop("Check failed: ", cluster_res, " unfiltered enrichment does not reproduce ", basename(pub_csv))
  }
  message(sprintf("[47] %s: check passed, unfiltered run reproduces the published table (%d terms)", cluster_res, nrow(chk_df)))

  # Main result: -6 exclusion applied.
  sets <- low_genes_and_universe(cluster_res, apply_filter = TRUE)
  go <- run_go(sets$low_genes, sets$universe)
  res_df <- go$res_df
  message(sprintf("[47] %s: %d terms with the -6 filter (%d unfiltered; %d shared)",
                  cluster_res, nrow(res_df), nrow(chk_df), length(intersect(res_df$ID, chk_df$ID))))
  write.csv(res_df, file.path(out_dir, sprintf("%s_go_res_df.csv", cluster_res)), row.names = FALSE)

  enrichres <- new("enrichResult", readable = FALSE, result = res_df, pvalueCutoff = 0.05,
                   pAdjustMethod = "BH", qvalueCutoff = 0.2, organism = "human",
                   ontology = "UNKNOWN", gene = sets$low_genes, keytype = "UNKNOWN",
                   universe = unique(go$bg_genes$gene), gene2Symbol = character(0), geneSets = go$bg_genes)
  p <- dotplot(enrichres, showCategory = n_show) +
    scale_color_gradient(low = "grey85", high = "#8B0000") +
    labs(title = paste0(cluster_res, "_consistently low-dispersed genes"),
         x = "Gene ratio", y = NULL, color = expression(-log[10](adjusted~p))) +
    theme_classic(base_family = "Arial") +
    theme(
      text = element_text(size = 6),
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 6),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      plot.title = element_text(size = 6, hjust = 0.5),
      legend.key.height = unit(3, "mm"),
      legend.key.width = unit(3, "mm"),
      plot.margin = margin(2, 2, 2, 2, unit = "mm"),
      aspect.ratio = 2)
  ggsave(file.path(out_dir, paste0(cluster_res, "_consistently_low_disp_GO_enrich.pdf")), p,
         device = cairo_pdf, width = 4, height = 5, units = "in")
}
