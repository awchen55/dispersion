#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(limma)              # normalizeQuantiles
  library(variancePartition)  # dream, makeContrastsDream, eBayes, topTable
  library(BiocParallel)       # SnowParam
  library(clusterProfiler)    # enricher, read.gmt
  library(enrichplot)         # dotplot method for enrichResult
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
fig_repro_dir <- file.path(repo_dir, "figure_reproduction")
data_dir <- file.path(fig_repro_dir, "data")
cluster_dir <- file.path(repo_dir, "cellranger_cluster-mode_trial")

mean_csv <- file.path(data_dir, "mean_expression_excl_cluster24.csv")

# ANALYSIS_OK[external-path]: read-only path to the MSigDB C5 GO:BP GMT file
# already used by analysis/pathway_analysis/clusterProfiler.Rmd (the source
# of the GOBP_* terms in GO_enrichments_results.xlsx); existence asserted below.
gmt_file <- file.path(cluster_dir, "analysis", "datasets", "gsea_gene_sets",
                     "background_genes", "c5.go.bp.v2024.1.Hs.symbols.gmt")

out_dir <- file.path(fig_repro_dir, "outputs")
ranked_dir <- file.path(out_dir, "ranked")
refit_dir <- file.path(out_dir, "dd_refit")
fig2_dir <- file.path(out_dir, "fig2_replication")
retained_matrix_csv <- file.path(out_dir, "retained_gene_matrix.csv")

ribo_pattern <- "^RP"

# ANALYSIS_OK[model-formula]: same dream mixed-model design as
# 03_rerun_dd_rank_analyses.R (cell_type fixed, cell_line random), applied here
# to the mean matrix so DE and DD are fit identically.
dream_formula <- ~ 0 + cell_type
n_dream_cores <- 4
min_universe_genes <- 100

individuals <- c("NA19099", "NA19093", "NA18522")
cell_types <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)
n_expected_cols <- length(individuals) * length(cell_types)

# Fig 2D's example contrast (source column suffix "_K"). Confirmed via
# Supplementary_table_s1.xlsx's contrast_reference sheet: K =
# compare_Foregut_endodermal_cells_Cardiomyocytes. Named identically to
# 03_rerun_dd_rank_analyses.R's contrast_matrix (paste0("compare_", x[2], "_", x[1])).
example_contrast <- "compare_Foregut_endodermal_cells_Cardiomyocytes"
dd_sig_thresh <- 0.05  # source: adj.P.Val <= 0.05 throughout

# Fig 2A GO enrichment (source: clusterProfiler.Rmd padj_cutoff/genecount_cutoff).
go_padj_cutoff <- 0.05
go_genecount_cutoff <- 5
n_top_terms_fig2a <- 8   # main_figs_1-4_related_supp.Rmd line 481: head(8)
n_top_terms_part3 <- 15  # clusterProfiler.Rmd's dotplot(..., showCategory = 15)

# ------------------------------- HELPERS -------------------------------------
source(file.path(fig_repro_dir, "scripts", "utils.R"))

dir.create(fig2_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(gmt_file)) stop("Missing GMT background file: ", gmt_file)
gmt_df <- read.gmt(gmt_file)
message(sprintf("[load] GMT: %d gene-term pairs, %d unique terms",
                nrow(gmt_df), length(unique(gmt_df$term))))

# Run enricher() against the GMT background restricted to universe_genes, then
# apply the source notebook's post-hoc filter (p.adjust < go_padj_cutoff & Count
# > go_genecount_cutoff), on top of enricher()'s own default pvalueCutoff/
# qvalueCutoff. Returns the enrichResult object (nrow(result) already filtered).
# ANALYSIS_OK[local-helper]: used 3x below (Fig 2A, Part 3 full, Part 3 retained); R026 cross-file rule cannot see intra-file callers.
run_enricher_gmt <- function(gene_list, universe_genes, label) {
  bg <- gmt_df[gmt_df$gene %in% universe_genes, ]
  message(sprintf("[go:%s] %d query genes; background %d genes, %d gene-term pairs",
                  label, length(gene_list), length(unique(bg$gene)), nrow(bg)))
  res <- enricher(gene_list, TERM2GENE = bg)
  if (is.null(res)) {
    message(sprintf("[go:%s] no enriched terms", label))
    return(res)
  }
  res@result <- res@result[res@result$p.adjust < go_padj_cutoff &
                            res@result$Count > go_genecount_cutoff, ]
  message(sprintf("[go:%s] %d terms pass p.adjust<%.2f & Count>%d",
                  label, nrow(res@result), go_padj_cutoff, go_genecount_cutoff))
  res
}

# --- Load mean matrix, drop ^RP genes (mirrors 03's dispersion loading) ------
if (!file.exists(mean_csv)) stop("Missing mean matrix: ", mean_csv)
mean_df <- read.csv(mean_csv, row.names = 1, check.names = FALSE)
if (ncol(mean_df) != n_expected_cols) {
  stop("Mean matrix has ", ncol(mean_df), " cols; expected ", n_expected_cols)
}
mean_mat <- as.matrix(mean_df)
# ANALYSIS_OK[filtering]: drop ^RP ribosomal genes (repo convention); retained count logged immediately below
mean_mat <- mean_mat[!grepl(ribo_pattern, rownames(mean_mat)), , drop = FALSE]
all_genes <- rownames(mean_mat)
message(sprintf("[load] %d non-ribo genes in mean matrix", length(all_genes)))

# --- Load per-CT retained indicator from script 01 ---------------------------
if (!file.exists(retained_matrix_csv)) {
  stop("Missing retained_gene_matrix.csv — run 01_exclude_topexpressed_genes.R first")
}
retained_ind <- as.matrix(read.csv(retained_matrix_csv, row.names = 1, check.names = FALSE))
# ANALYSIS_OK[column-reorder]: reorder to canonical cell-type order only; no genes/rows dropped
retained_ind <- retained_ind[, cell_types, drop = FALSE]
union_retained_genes <- all_genes[rowSums(retained_ind[all_genes, , drop = FALSE] == 1) == length(cell_types)]
message(sprintf("[load] union-retained genes (excluded in no cell type): %d of %d",
                length(union_retained_genes), length(all_genes)))
if (length(union_retained_genes) < min_universe_genes) {
  stop("Union-retained gene set too small (", length(union_retained_genes), ")")
}

# =============================================================================
# PART 1 — Fig 2D: DE (mean) vs DD (dispersion) scatter, contrast K
# =============================================================================

sample_cols <- colnames(mean_mat)
col_ind <- sub("_.*$", "", sample_cols)
col_ct  <- sub("^[^_]*_", "", sample_cols)
stopifnot(all(col_ind %in% individuals), all(col_ct %in% cell_types),
          length(sample_cols) == n_expected_cols)
metadata <- data.frame(
  cell_type = factor(col_ct, levels = cell_types),
  cell_line = factor(col_ind, levels = individuals),
  row.names = sample_cols)

# 15 pairwise contrasts named compare_{B}_{A} — identical construction to
# 03_rerun_dd_rank_analyses.R's contrast_matrix.
pairs <- combn(cell_types, 2, simplify = FALSE)
contrast_matrix <- setNames(
  vapply(pairs, function(x) paste0("cell_type", x[2], " - ", "cell_type", x[1]),
         character(1)),
  vapply(pairs, function(x) paste0("compare_", x[2], "_", x[1]), character(1)))
stopifnot(example_contrast %in% names(contrast_matrix))

# One dream + eBayes fit on the MEAN matrix for a given gene set — mirrors
# 03_rerun_dd_rank_analyses.R's fit_dream_dd(), applied to mean instead of
# dispersion so DE and DD are fit under the identical pipeline (paired
# full-vs-retained design: only the gene set differs between the two fits).
# ANALYSIS_OK[local-helper]: used within this standalone pipeline-stage script (run via run.sh); not a shared-library export
fit_dream_de <- function(gene_set, tag) {
  mat <- normalizeQuantiles(mean_mat[gene_set, sample_cols, drop = FALSE])
  contrast_dream <- makeContrastsDream(dream_formula, metadata,
                                       contrasts = contrast_matrix)
  param <- SnowParam(n_dream_cores, type = "SOCK",
                     progressbar = FALSE, exportglobals = FALSE)
  message(sprintf("[deB] dream fit '%s': %d genes x %d samples, %d contrasts",
                  tag, length(gene_set), length(sample_cols),
                  length(contrast_matrix)))
  fit <- dream(mat, dream_formula, metadata, contrast_dream, BPPARAM = param)
  eBayes(fit)
}

fit_de_full     <- fit_dream_de(all_genes, "full")
fit_de_retained <- fit_dream_de(union_retained_genes, "retained")
saveRDS(fit_de_full,     file.path(refit_dir, "dream_fit_mean_full.rds"))
saveRDS(fit_de_retained, file.path(refit_dir, "dream_fit_mean_retained.rds"))

# DD fits already produced by 03_rerun_dd_rank_analyses.R Part B.
dd_full_rds     <- file.path(refit_dir, "dream_fit_full.rds")
dd_retained_rds <- file.path(refit_dir, "dream_fit_retained.rds")
if (!file.exists(dd_full_rds) || !file.exists(dd_retained_rds)) {
  stop("Missing dream_fit_full.rds / dream_fit_retained.rds — run 03_rerun_dd_rank_analyses.R first")
}
fit_dd_full     <- readRDS(dd_full_rds)
fit_dd_retained <- readRDS(dd_retained_rds)

# topTable for the example contrast from a DD and a DE fit sharing the same
# gene universe. Returns a list(dd_tt=, de_tt=) restricted to shared genes.
# ANALYSIS_OK[local-helper]: used twice below (full, retained); R026 cross-file rule cannot see intra-file callers.
get_contrast_toptables <- function(fit_de, fit_dd, tag) {
  de_tt <- topTable(fit_de, coef = example_contrast, number = Inf, sort.by = "none")
  dd_tt <- topTable(fit_dd, coef = example_contrast, number = Inf, sort.by = "none")
  shared <- intersect(rownames(de_tt), rownames(dd_tt))
  message(sprintf("[fig2d:%s] %d genes shared between DE and DD fits", tag, length(shared)))
  list(de_tt = de_tt[shared, ], dd_tt = dd_tt[shared, ])
}

# Fig 2D plot_df, mirroring source lines 731-738 exactly (including the sign
# flip on both axes).
# ANALYSIS_OK[local-helper]: used twice below (full, retained); R026 cross-file rule cannot see intra-file callers.
build_de_dd_df <- function(tt, tag) {
  data.frame(
    gene = rownames(tt$de_tt),
    de_logFC = -tt$de_tt$logFC,
    dv_logFC = -tt$dd_tt$logFC,
    dv_sig = factor(tt$dd_tt$adj.P.Val <= dd_sig_thresh, levels = c(FALSE, TRUE)),
    de_sig = tt$de_tt$adj.P.Val <= dd_sig_thresh,
    stringsAsFactors = FALSE)
}

# Source lines 744-799, byte-for-byte aesthetic (labels/colors/theme).
# ANALYSIS_OK[local-helper]: used twice below (full, retained); R026 cross-file rule cannot see intra-file callers.
make_de_dd_scatter <- function(plot_df) {
  ggplot(plot_df, aes(x = de_logFC, y = dv_logFC, color = dv_sig)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey40") +
    geom_point(size = 0.8, alpha = 0.8) +
    scale_color_manual(
      values = c("FALSE" = "grey70", "TRUE" = "#B22222"),
      labels = c("FALSE" = "Non DD", "TRUE" = "DD")
    ) +
    labs(x = "DE log2 fold change", y = "DD log2 fold change", color = NULL) +
    theme_classic(base_family = "Arial") +
    theme(
      text = element_text(size = 6),
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 6),
      axis.line = element_line(linewidth = 0.4),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
      legend.position = "right",
      legend.text = element_text(size = 6),
      legend.key.size = unit(6, "pt"),
      plot.margin = margin(2, 2, 2, 2, unit = "mm")
    )
}

tt_full     <- get_contrast_toptables(fit_de_full, fit_dd_full, "full")
tt_retained <- get_contrast_toptables(fit_de_retained, fit_dd_retained, "retained")
df_full     <- build_de_dd_df(tt_full, "full")
df_retained <- build_de_dd_df(tt_retained, "retained")

write.csv(df_full,     file.path(fig2_dir, "fig2d_scatter_full.csv"), row.names = FALSE)
write.csv(df_retained, file.path(fig2_dir, "fig2d_scatter_retained.csv"), row.names = FALSE)

p_full     <- make_de_dd_scatter(df_full)
p_retained <- make_de_dd_scatter(df_retained)

# Dimensions match the published panel (pdfinfo on reference_material/
# paper_figure_scripts' sibling first_pub_figures/Fig2c_cardvsforend_DD_vs_DE_logFC.pdf:
# 252x144pt = 3.5x2in) -- the source Rmd chunk never calls ggsave() for this
# panel itself, so there's no dimension to copy from the notebook directly.
ggsave(file.path(fig2_dir, "fig2d_scatter_full.pdf"), p_full,
       width = 3.5, height = 2, device = cairo_pdf)
ggsave(file.path(fig2_dir, "fig2d_scatter_retained.pdf"), p_retained,
       width = 3.5, height = 2, device = cairo_pdf)

message(sprintf(
  "[fig2d] full: n=%d, Spearman(DE,DD)=%.3f | retained: n=%d, Spearman(DE,DD)=%.3f",
  nrow(df_full), cor(df_full$de_logFC, df_full$dv_logFC, method = "spearman"),
  nrow(df_retained), cor(df_retained$de_logFC, df_retained$dv_logFC, method = "spearman")))

# =============================================================================
# PART 2 — Fig 2A: GO enrichment of consistently-low-dispersion genes
# =============================================================================

low_genes_file <- file.path(ranked_dir, "consistently_low_retained.txt")
if (!file.exists(low_genes_file)) {
  stop("Missing ", low_genes_file, " — run 02_rerun_ranked_dispersion.R first")
}
low_genes <- readLines(low_genes_file)
message(sprintf("[fig2a] %d consistently-low-dispersion (retained) genes", length(low_genes)))

ego <- run_enricher_gmt(low_genes, union_retained_genes, "fig2a_low_dispersion")
ego_df <- if (is.null(ego)) data.frame() else as.data.frame(ego)
write.csv(ego_df, file.path(fig2_dir, "fig2a_go_enrichment_retained_full.csv"), row.names = FALSE)

if (nrow(ego_df) > 0) {
  # Source lines 477-536, byte-for-byte (dotplot construction + aesthetic).
  # ANALYSIS_OK[plot-filter]: head(n_top_terms_fig2a) restricts the DOTPLOT to
  # the top 8 terms by p.adjust (matches source line 481 exactly); the full
  # enrichment table (all significant terms) is already written unfiltered to
  # fig2a_go_enrichment_retained_full.csv above -- this filter never touches
  # the analysis population, only what's plotted.
  plot_df <- ego_df %>%
    arrange(p.adjust) %>%
    head(n_top_terms_fig2a) %>%
    mutate(
      # ANALYSIS_OK[fixed-width-split]: GeneRatio is always clusterProfiler's
      # "k/n" string; str_split_fixed(..., 2) is guaranteed exactly 2 columns,
      # so [,1]/[,2] index the numerator/denominator, not an arbitrary
      # dataframe column position (matches source lines 483-485 exactly).
      GeneRatio_num = as.numeric(str_split_fixed(GeneRatio, "/", 2)[, 1]) /
                      as.numeric(str_split_fixed(GeneRatio, "/", 2)[, 2]),
      log10_padj = -log10(p.adjust),
      Description = factor(ID, levels = rev(ID))
    )
  write.csv(plot_df, file.path(fig2_dir, "fig2a_go_enrichment_retained_top8.csv"), row.names = FALSE)

  p_go <- ggplot(plot_df, aes(x = GeneRatio_num, y = Description)) +
    geom_point(aes(size = Count, color = log10_padj)) +
    scale_color_gradient(low = "grey85", high = "#8B0000",
                         name = expression(-log[10](adjusted~p))) +
    scale_size(range = c(2, 6)) +
    labs(
      title = "GO enrichment: consistently low-dispersed genes",
      x = "Gene ratio", y = NULL
    ) +
    theme_classic(base_family = "Arial") +
    theme(
      text = element_text(size = 6),
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 6),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      plot.title = element_text(size = 6, hjust = 0.5),
      legend.key.height = unit(3, "mm"),
      legend.key.width  = unit(3, "mm"),
      plot.margin = margin(2, 2, 2, 2, unit = "mm"),
      aspect.ratio = 2
    )

  # The published panel (first_pub_figures/Fig1d_consistently_low_disp_GO_
  # enrich.pdf, 288x360pt = 4x5in) is portrait, but at that size here the
  # theme's aspect.ratio=2 + 8 long GO-term labels collapse the panel to a
  # sliver and all y-axis text overlaps on one line (Brendan flagged as
  # squished again after the first page-size-matched fix). Landscape 5x4in,
  # found by rendering and eyeballing several candidates, is the smallest
  # size at which every label is legible with no overlap; not a page-size
  # match since the source Rmd never calls ggsave() for this panel.
  ggsave(file.path(fig2_dir, "fig2a_go_enrichment_retained.pdf"), p_go,
         width = 5, height = 4, device = cairo_pdf)
}

# =============================================================================
# PART 3 — Cardiomyocytes vs Foregut_endodermal_cells: "less variable AND
# higher mean in Cardiomyocytes" GO enrichment (GO_enrichments_results.xlsx
# "DD_DE_GO_enrichment" sheet; source: pathway_analysis/clusterProfiler.Rmd
# lines 64-160). Background restricted to genes upregulated in Cardiomyocytes
# for this contrast (source's "up_bg" design), independent of variance
# direction, to control for mean-expression-driven detection bias. Computed
# for both the full and retained gene universes for comparison.
# =============================================================================

# ANALYSIS_OK[local-helper]: used twice below (full, retained); R026 cross-file rule cannot see intra-file callers.
run_part3 <- function(tt, tag) {
  de_tt <- tt$de_tt
  dd_tt <- tt$dd_tt
  # Contrast is Foregut - Cardiomyocytes: logFC < 0 means higher mean in Card;
  # dispersion logFC > 0 means less variable in Card (source lines 72-83).
  mean_upreg_in_card <- de_tt$logFC < 0 & de_tt$adj.P.Val <= dd_sig_thresh
  less_var_in_card    <- dd_tt$logFC > 0 & dd_tt$adj.P.Val <= dd_sig_thresh
  background_genes <- rownames(de_tt)[mean_upreg_in_card]
  goi <- rownames(de_tt)[mean_upreg_in_card & less_var_in_card]
  message(sprintf(
    "[part3:%s] %d genes upregulated in Card (background); %d of those also less variable in Card (foreground)",
    tag, length(background_genes), length(goi)))

  res <- run_enricher_gmt(goi, background_genes, paste0("part3_", tag))
  res_df <- if (is.null(res)) data.frame() else as.data.frame(res)
  write.csv(res_df, file.path(fig2_dir, sprintf("fig2_part3_card_vs_foregut_%s.csv", tag)),
            row.names = FALSE)
  if (nrow(res_df) > 0) {
    p <- dotplot(res, showCategory = n_top_terms_part3) +
      ggtitle(sprintf("Mean increased, variability decreased (%s)", tag))
    # 5x5in (matching Fig2d_cardvsforegut_GO_enrich.pdf's page size) still
    # squishes this dotplot's longer term labels (Brendan flagged again after
    # the page-size-matched fix). Reverted to 7x7in, the pre-2026-07-23-fix
    # size, which rendered legibly.
    ggsave(file.path(fig2_dir, sprintf("fig2_part3_card_vs_foregut_%s.pdf", tag)),
           p, width = 7, height = 7)
  }
  res_df
}

part3_full     <- run_part3(tt_full, "full")
part3_retained <- run_part3(tt_retained, "retained")

message("[done] Parts 1-3 written to ", fig2_dir)
