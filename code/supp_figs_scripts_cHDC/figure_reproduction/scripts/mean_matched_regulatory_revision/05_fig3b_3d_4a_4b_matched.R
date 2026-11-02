#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(readxl)
  library(variancePartition)  # topTable
})

args       <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
source(file.path(script_dir, "config.R"))
panel_dir <- ensure_dir(file.path(out_dir, "fig3_fig4_panels"))

TWN_S4_XLSX  <- file.path(src_dir, "analysis/datasets/external_data/saha_coexp_net_TWN_2017/Supplemental_Data_S4.xlsx")
PPI_TSV      <- file.path(src_dir, "analysis/datasets/external_data/protein_interaction_data/interaction_consensus.tsv.zip")
PPI_TISSUE_TSV <- file.path(src_dir, "analysis/datasets/external_data/protein_interaction_data/has_data_Tissue.tsv")
# Script 04's DE fit (Epithelial sub-cluster 24 excluded, union-retained genes).
DE_DREAM_RDS <- file.path(repo_dir, "figure_reproduction/outputs/dd_refit/dream_fit_mean_retained.rds")
EPD_TATA_BED <- file.path(src_dir, "analysis/datasets/external_data/epd/human_epdnew_3sa0y.bed")

# Fixed at 0 in the source Rmd (Rmd 2226) -- a gene counts as DE in a contrast
# if adj.P.Val <= DE_ADJP, regardless of logFC magnitude (any nonzero direction).
DE_FC_THRESH <- 0
DE_ADJP      <- 0.05  # Rmd 2240's literal 0.05, named here for clarity

# ---- Inputs -------------------------------------------------------------
gene_sets <- readRDS(file.path(out_dir, "gene_sets", "dispersion_gene_sets.rds"))
mm        <- readRDS(file.path(out_dir, "mean_matched", "matched_gene_sets.rds"))
all_genes_bg <- rownames(gene_sets$var_quintile_dfs)  # ribo-filtered, -6-retained universe, unmatched

# Low-vs-High labelling of a matched (first/fifth/never) triple, restricted
# to Low/High (3B/3D/4A/4B never use "never" -- only the eGene/s_het/exon
# panels in 03 do).
# ANALYSIS_OK[within-file-helper]: called below for every seed; local to this
# script, not shared, so no cross-file caller by design.
label_low_high <- function(genes_list) {
  bind_rows(
    data.frame(gene = genes_list$first, rank_category = "Low-dispersion"),
    data.frame(gene = genes_list$fifth, rank_category = "High-dispersion"))
}

# Shared boxplot styling, faithful to the source Rmd's exact theme for these
# 4 panels (Rmd 1477-1531 / 1676-1732 / 2256-2311 -- all identical).
# ANALYSIS_OK[within-file-helper]: called below for every panel; local to
# this script, not shared, so no cross-file caller by design.
matched_boxplot <- function(df, y_col, y_lab, title, p_value) {
  # Explicit factor order -- without this, ggplot's default alphabetical
  # ordering puts "High-dispersion" before "Low-dispersion" on the x-axis,
  # inconsistent with the source Rmd's convention (Low always first) and
  # with 06's Fig3C panels (which do set this explicitly).
  df$rank_category <- factor(df$rank_category, levels = c("Low-dispersion", "High-dispersion"))
  ggplot(df, aes(x = rank_category, y = .data[[y_col]], fill = rank_category)) +
    geom_boxplot(width = 0.6, linewidth = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 0.6, alpha = 0.4) +
    # p-value annotation, same style/wording as 4B's forest-plot labels
    # (format_pvalue_label(), config.R) -- Brendan's request, 2026-09-03.
    # y = Inf + vjust pins it just inside the top of the panel regardless of
    # this feature's own y-range, matching how 4B's label sits above its points.
    annotate("text", x = 1.5, y = Inf, label = format_pvalue_label(p_value),
             vjust = 1.3, size = 2, family = "Arial") +
    scale_fill_manual(values = c("Low-dispersion" = "steelblue", "High-dispersion" = "coral")) +
    scale_x_discrete(labels = c("Low-dispersion" = "Low-dispersion\ngenes",
                                 "High-dispersion" = "High-dispersion\ngenes")) +
    labs(x = NULL, y = y_lab, title = title) +
    theme_classic(base_family = "Arial") +
    theme(text = element_text(size = 6), axis.title = element_text(size = 6),
          axis.text = element_text(size = 6), axis.line = element_line(linewidth = 0.4),
          axis.text.x = element_text(hjust = 0.5), legend.position = "none",
          plot.margin = margin(2, 2, 2, 2, unit = "mm"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4))
}

# Seed-stability check shared by all 4 panels: recompute the Low-vs-High
# Wilcoxon effect under every seed, assert the sign of the shift agrees.
# ANALYSIS_OK[within-file-helper]: called below for 3B/3D/4A; local to this
# script, not shared, so no cross-file caller by design.
seed_stability_wilcox <- function(feature_name, value_lookup_fn) {
  bind_rows(lapply(names(mm$matched_by_seed), function(s) {
    lab <- label_low_high(mm$matched_by_seed[[s]])
    v <- value_lookup_fn(lab$gene)
    keep <- !is.na(v)
    r <- wilcox_effect(v[keep], lab$rank_category[keep],
                       "Low-dispersion", "High-dispersion", feature_name)
    r$seed <- s; r
  }))
}

primary <- mm$matched_by_seed[[as.character(mm$match_seed_primary)]]
primary_labelled <- label_low_high(primary)

# =============================================================================
# Fig3B: TWN shared-hub connectivity rank, Low vs High (matched)
# =============================================================================
log_step("=== Fig3B: TWN shared-hub rank ===")
assert_that(file.exists(TWN_S4_XLSX), paste("missing:", TWN_S4_XLSX))
s4_sheets <- excel_sheets(TWN_S4_XLSX)
assert_that(grepl("^Group3", s4_sheets[10]),
            sprintf("Supplemental_Data_S4.xlsx sheet 10 is '%s', expected 'Group3 (TE-TE)' -- file layout changed",
                    s4_sheets[10]))
shared_te_te <- read_excel(TWN_S4_XLSX, sheet = s4_sheets[10], col_names = TRUE)
shared_te_te$rank <- as.integer(shared_te_te$rank)

twn_lookup <- function(genes) shared_te_te$rank[match(genes, shared_te_te$hub)]
twn_frame <- primary_labelled %>%
  filter(gene %in% shared_te_te$hub) %>%
  mutate(hub_rank = twn_lookup(gene))
log_step("matched genes that are TWN shared hubs", nrow(twn_frame))

twn_res <- wilcox_effect(twn_frame$hub_rank, twn_frame$rank_category,
                          "Low-dispersion", "High-dispersion", "twn_hub_rank")
write.csv(twn_res, file.path(panel_dir, "fig3b_twn_hub_rank_wilcox_matched.csv"), row.names = FALSE)

twn_stability <- seed_stability_wilcox("twn_hub_rank", twn_lookup)
write.csv(twn_stability, file.path(panel_dir, "fig3b_twn_hub_rank_seed_stability.csv"), row.names = FALSE)
assert_that(length(unique(sign(twn_stability$hl_shift))) == 1,
            "Fig3B Low-vs-High effect changes SIGN across matching seeds -- not stable")
log_step("Fig3B Low-vs-High effect direction is stable across all matching seeds")

twn_plot <- matched_boxplot(twn_frame, "hub_rank", "Shared hub genes ranking",
                             "Mean-matched TWN shared-hub rank by dispersion class",
                             p_value = twn_res$p_value)
ggsave(file.path(panel_dir, "fig3b_twn_hub_rank_matched_boxplot.pdf"), twn_plot,
       width = 3.2, height = 2.8, units = "in", device = cairo_pdf)
saveRDS(twn_plot, file.path(panel_dir, "fig3b_twn_hub_rank_plot.rds"))  # for 07's combined-figure patchwork

# =============================================================================
# Fig3D: protein-protein interaction connectivity, Low vs High (matched)
# =============================================================================
log_step("=== Fig3D: PPI connectivity ===")
assert_that(file.exists(PPI_TSV), paste("missing:", PPI_TSV))
assert_that(file.exists(PPI_TISSUE_TSV), paste("missing:", PPI_TISSUE_TSV))
ppi <- read_tsv(PPI_TSV, show_col_types = FALSE)
exp_evidence <- read_tsv(PPI_TISSUE_TSV, show_col_types = FALSE)
# ANALYSIS_OK[column-add]: adds a gene_1 symbol column via lookup join; no rows
# dropped or reordered (row count logged via connectivity's nrow() below).
ppi <- ppi %>%
  mutate(gene_1 = exp_evidence$Gene[match(ensembl_gene_id_1, exp_evidence$Ensembl)])
connectivity <- as.data.frame(table(ppi$gene_1))
colnames(connectivity) <- c("Gene", "Freq")
log_step("PPI connectivity table built (genes with >=1 interaction)", nrow(connectivity))

ppi_lookup <- function(genes) connectivity$Freq[match(genes, connectivity$Gene)]
ppi_frame <- primary_labelled %>%
  filter(gene %in% connectivity$Gene) %>%
  mutate(int_count = ppi_lookup(gene))
log_step("matched genes with a PPI connectivity count", nrow(ppi_frame))

ppi_res <- wilcox_effect(ppi_frame$int_count, ppi_frame$rank_category,
                          "Low-dispersion", "High-dispersion", "ppi_int_count")
write.csv(ppi_res, file.path(panel_dir, "fig3d_ppi_connectivity_wilcox_matched.csv"), row.names = FALSE)

ppi_stability <- seed_stability_wilcox("ppi_int_count", ppi_lookup)
write.csv(ppi_stability, file.path(panel_dir, "fig3d_ppi_connectivity_seed_stability.csv"), row.names = FALSE)
assert_that(length(unique(sign(ppi_stability$hl_shift))) == 1,
            "Fig3D Low-vs-High effect changes SIGN across matching seeds -- not stable")
log_step("Fig3D Low-vs-High effect direction is stable across all matching seeds")

ppi_plot <- matched_boxplot(ppi_frame, "int_count", "Protein connections per gene",
                             "Mean-matched PPI connectivity by dispersion class",
                             p_value = ppi_res$p_value)
ggsave(file.path(panel_dir, "fig3d_ppi_connectivity_matched_boxplot.pdf"), ppi_plot,
       width = 3.2, height = 2.8, units = "in", device = cairo_pdf)
saveRDS(ppi_plot, file.path(panel_dir, "fig3d_ppi_connectivity_plot.rds"))  # for 07's combined-figure patchwork

# =============================================================================
# Fig4A: DE (mean-expression) count across 15 pairwise contrasts, Low vs High
# =============================================================================
log_step("=== Fig4A: DE count ===")
assert_that(file.exists(DE_DREAM_RDS), paste("missing:", DE_DREAM_RDS))
de_dream_obj <- readRDS(DE_DREAM_RDS)
contrast_coefs <- colnames(de_dream_obj$coefficients)[1:15]
assert_that(all(grepl("^compare_", contrast_coefs)),
            "first 15 DE dream coefficients are not all 'compare_*' contrasts -- fit layout changed")

contrast_de <- matrix(NA, nrow = nrow(de_dream_obj$coefficients), ncol = length(contrast_coefs),
                       dimnames = list(rownames(de_dream_obj$coefficients), contrast_coefs))
for (contrast in contrast_coefs) {
  tt <- topTable(de_dream_obj, coef = contrast, number = Inf, sort.by = "none")
  contrast_de[, contrast] <- tt$adj.P.Val <= DE_ADJP & abs(tt$logFC) >= DE_FC_THRESH
}
de_count_lookup <- function(genes) {
  idx <- match(genes, rownames(contrast_de))
  vapply(idx, function(i) if (is.na(i)) NA_real_ else sum(contrast_de[i, ]), numeric(1))
}
de_frame <- primary_labelled %>%
  filter(gene %in% rownames(contrast_de)) %>%
  mutate(DE_count = de_count_lookup(gene))
log_step("matched genes with a DE-count value", nrow(de_frame))

de_res <- wilcox_effect(de_frame$DE_count, de_frame$rank_category,
                         "Low-dispersion", "High-dispersion", "DE_count")
write.csv(de_res, file.path(panel_dir, "fig4a_de_count_wilcox_matched.csv"), row.names = FALSE)

de_stability <- seed_stability_wilcox("DE_count", de_count_lookup)
write.csv(de_stability, file.path(panel_dir, "fig4a_de_count_seed_stability.csv"), row.names = FALSE)
assert_that(length(unique(sign(de_stability$hl_shift))) == 1,
            "Fig4A Low-vs-High effect changes SIGN across matching seeds -- not stable")
log_step("Fig4A Low-vs-High effect direction is stable across all matching seeds")

de_plot <- matched_boxplot(de_frame, "DE_count", "DE count per gene",
                            "Mean-matched DE count by dispersion class",
                            p_value = de_res$p_value)
ggsave(file.path(panel_dir, "fig4a_de_count_matched_boxplot.pdf"), de_plot,
       width = 3.2, height = 2.8, units = "in", device = cairo_pdf)
saveRDS(de_plot, file.path(panel_dir, "fig4a_de_count_plot.rds"))  # for 07's combined-figure patchwork

# =============================================================================
# Fig4B: TATA-box promoter enrichment, each of Low/High vs full background
# =============================================================================
log_step("=== Fig4B: TATA enrichment ===")
assert_that(file.exists(EPD_TATA_BED), paste("missing:", EPD_TATA_BED))
epd_tata <- read.table(EPD_TATA_BED, sep = "\t", header = FALSE,
                        col.names = c("chr", "start", "end", "promoter_id", "score", "strand"))
epd_tata$gene <- sub("_[0-9]+$", "", epd_tata$promoter_id)
tata_genes <- intersect(unique(epd_tata$gene), all_genes_bg)
log_step("TATA-annotated genes in the ribo-filtered background", length(tata_genes))

# Fisher's exact test of one target gene set against TATA annotation, vs the
# SAME full (unmatched) background the source Rmd used (Rmd 2331-2351) --
# mean-matching is about Low vs High comparability with each other, not about
# changing what background TATA enrichment is assessed against.
# ANALYSIS_OK[within-file-helper]: called below for Low/High and every seed;
# local to this script, not shared, so no cross-file caller by design.
tata_enrich <- function(target_genes) {
  a <- length(intersect(tata_genes, target_genes))
  b <- length(setdiff(tata_genes, target_genes))
  c <- length(setdiff(target_genes, tata_genes))
  d <- length(setdiff(all_genes_bg, union(tata_genes, target_genes)))
  ct <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
               dimnames = list(TATA = c("Yes", "No"), Target = c("In_Set", "Not_In_Set")))
  fisher.test(ct)
}

low_ft  <- tata_enrich(primary$first)
high_ft <- tata_enrich(primary$fifth)
tata_df <- data.frame(
  gene_set = factor(c("Low-dispersion", "High-dispersion"), levels = c("Low-dispersion", "High-dispersion")),
  odds_ratio = c(unname(low_ft$estimate), unname(high_ft$estimate)),
  ci_lower   = c(low_ft$conf.int[1], high_ft$conf.int[1]),
  ci_upper   = c(low_ft$conf.int[2], high_ft$conf.int[2]),
  p_value    = c(low_ft$p.value, high_ft$p.value)
)
tata_df$log2_or       <- log2(tata_df$odds_ratio)
tata_df$log2_ci_lower <- log2(tata_df$ci_lower)
tata_df$log2_ci_upper <- log2(tata_df$ci_upper)
write.csv(tata_df, file.path(panel_dir, "fig4b_tata_fisher_matched.csv"), row.names = FALSE)
print(tata_df)

tata_stability <- bind_rows(lapply(names(mm$matched_by_seed), function(s) {
  m <- mm$matched_by_seed[[s]]
  lo <- tata_enrich(m$first); hi <- tata_enrich(m$fifth)
  data.frame(seed = s,
             log2_or_low  = log2(unname(lo$estimate)), p_low  = lo$p.value,
             log2_or_high = log2(unname(hi$estimate)), p_high = hi$p.value)
}))
write.csv(tata_stability, file.path(panel_dir, "fig4b_tata_seed_stability.csv"), row.names = FALSE)
assert_that(length(unique(sign(tata_stability$log2_or_low))) == 1,
            "Fig4B Low-dispersion TATA enrichment direction changes SIGN across matching seeds -- not stable")
assert_that(length(unique(sign(tata_stability$log2_or_high))) == 1,
            "Fig4B High-dispersion TATA enrichment direction changes SIGN across matching seeds -- not stable")
log_step("Fig4B Low/High TATA enrichment directions are each stable across all matching seeds")

tata_df$p_label <- format_pvalue_label(tata_df$p_value)
tata_plot <- ggplot(tata_df, aes(x = gene_set, y = log2_or)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey40") +
  geom_errorbar(aes(ymin = log2_ci_lower, ymax = log2_ci_upper), width = 0.15, linewidth = 0.4) +
  geom_point(size = 1.6) +
  geom_text(aes(y = log2_ci_upper, label = p_label), vjust = -0.6, size = 2, family = "Arial") +
  scale_x_discrete(labels = c("Low-dispersion" = "Low-dispersion\ngenes",
                               "High-dispersion" = "High-dispersion\ngenes")) +
  labs(x = NULL, y = expression(log[2]("Odds ratio")),
       title = "Mean-matched TATA-box enrichment by dispersion class") +
  theme_classic(base_family = "Arial") +
  theme(text = element_text(size = 6), axis.title = element_text(size = 6),
        axis.text = element_text(size = 6), axis.line = element_line(linewidth = 0.4),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        legend.position = "none", plot.margin = margin(2, 2, 2, 2, unit = "mm"))
ggsave(file.path(panel_dir, "fig4b_tata_matched_forestplot.pdf"), tata_plot,
       width = 3.2, height = 2.8, units = "in", device = cairo_pdf)
saveRDS(tata_plot, file.path(panel_dir, "fig4b_tata_plot.rds"))  # for 07's combined-figure patchwork

log_step(paste("Fig3B/3D/4A/4B mean-matched results written to", panel_dir))
