#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(readr)
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
cluster_dir <- file.path(repo_dir, "cellranger_cluster-mode_trial")
analysis_dir <- file.path(repo_dir, "figure_reproduction")
out_dir <- file.path(analysis_dir, "outputs")
ranked_dir <- file.path(out_dir, "ranked")
dd_rank_dir <- file.path(out_dir, "dd_rank")
filt_mat_dir <- file.path(out_dir, "filtered_matrices")
fig3_dir <- file.path(out_dir, "fig3_replication")
dir.create(fig3_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(analysis_dir, "scripts", "utils.R"))

individuals <- c("NA19099", "NA19093", "NA18522")
cell_types <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)

twn_s2_xlsx <- file.path(cluster_dir, "analysis", "datasets", "external_data",
                          "saha_coexp_net_TWN_2017", "Supplemental_Data_S2.xlsx")
twn_s4_xlsx <- file.path(cluster_dir, "analysis", "datasets", "external_data",
                          "saha_coexp_net_TWN_2017", "Supplemental_Data_S4.xlsx")
ppi_tsv <- file.path(cluster_dir, "analysis", "datasets", "external_data",
                      "protein_interaction_data", "interaction_consensus.tsv.zip")
ppi_exp_tsv <- file.path(cluster_dir, "analysis", "datasets", "external_data",
                          "protein_interaction_data", "has_data_Tissue.tsv")

# =============================================================================
# PART A -- Fig3A: eGene enrichment across global dispersion quintiles
# =============================================================================
egene_csv <- file.path(dd_rank_dir, "egene_quintile_fisher_retained.csv")
if (!file.exists(egene_csv)) {
  stop("Missing ", egene_csv, " -- run 03_rerun_dd_rank_analyses.R Part A first")
}
plot_df <- read.csv(egene_csv)
plot_df$quintile <- factor(plot_df$quintile, levels = sort(unique(plot_df$quintile)))
plot_df <- plot_df %>%
  mutate(label = ifelse(
    quintile == levels(quintile)[1] & p.value < 0.05,
    paste0("p = ", format(p.value, scientific = TRUE, digits = 2)),
    NA
  ))

p_egene <- ggplot(plot_df, aes(x = quintile, y = log2OR)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_errorbar(aes(ymin = log2_ci_lower, ymax = log2_ci_upper), width = 0.12, linewidth = 0.4) +
  geom_point(aes(color = p.value < 0.05), size = 1.8) +
  geom_text(aes(label = label), position = position_nudge(y = 1), size = 2, family = "Arial") +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "#B22222"), guide = "none") +
  labs(
    x = "Dispersion (global quintile rank)",
    y = expression(log[2]("Odds ratio")),
    title = "Enrichment of eGenes across dispersion quintiles"
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    text = element_text(size = 6),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    plot.title = element_text(size = 6, hjust = 0.5),
    axis.line = element_line(linewidth = 0.4),
    plot.margin = margin(2, 2, 2, 2, unit = "mm")
  )

ggsave(file.path(fig3_dir, "Fig3a_eQTL_fisher_global_ranks_excl_topexpressed.pdf"),
       p_egene, device = cairo_pdf, width = 3.54, height = 3.2, units = "in")
saveRDS(p_egene, file.path(fig3_dir, "fig3a_panel.rds"))
message("[fig3a] done")

# =============================================================================
# PART B -- Fig3B: TWN shared-hub gene rank, global consistently-low/-high
# =============================================================================
low_genes <- readLines(file.path(ranked_dir, "consistently_low_retained.txt"))
high_genes <- readLines(file.path(ranked_dir, "consistently_high_retained.txt"))
message(sprintf("[global] consistently-low: %d genes, consistently-high: %d genes",
                length(low_genes), length(high_genes)))

s4_sheets <- excel_sheets(twn_s4_xlsx)
stopifnot("Supplemental_Data_S4.xlsx sheet 10 layout changed -- expected 'Group3 (TE-TE)'" =
            grepl("^Group3", s4_sheets[10]))
shared_te_te <- read_excel(twn_s4_xlsx, sheet = s4_sheets[10], col_names = TRUE)
shared_te_te$rank <- as.integer(shared_te_te$rank)

df_b <- data.frame(genes = c(low_genes, high_genes),
                    rank = c(rep("1", length(low_genes)), rep("5", length(high_genes)))) %>%
  dplyr::filter(genes %in% shared_te_te$hub) %>%
  dplyr::mutate(hub_rank = shared_te_te$rank[match(genes, shared_te_te$hub)])
df_b$rank <- factor(df_b$rank, levels = c("1", "5"))
message(sprintf("[fig3b] %d/%d low, %d/%d high genes found as TWN shared hubs",
                sum(df_b$rank == "1"), length(low_genes),
                sum(df_b$rank == "5"), length(high_genes)))

p_twn_shared <- make_low_high_boxplot(df_b, "hub_rank", "Shared hub genes ranking")
ggsave(file.path(fig3_dir, "Fig3b_heart_skeletal_TWN_excl_topexpressed.pdf"),
       p_twn_shared, device = cairo_pdf, width = 2.2, height = 2.54, units = "in")
saveRDS(p_twn_shared, file.path(fig3_dir, "fig3b_panel.rds"))
message("[fig3b] done")

# =============================================================================
# PART C -- Fig3C: cell-type-specific TWN connections (own retained ranking)
# =============================================================================
s2_sheets <- excel_sheets(twn_s2_xlsx)

# Bottom/top `percentage` of a cell type's OWN retained (-6-excluded) median
# dispersion ranking (mirrors the source's ct_ranks_df, but computed on
# script 01's per-cell-type retained matrix instead of the full universe).
ct_specific_twn_percentage_retained <- function(cell_type, sheet, percentage) {
  disp_path <- file.path(filt_mat_dir, paste0(cell_type, "_dispersion_retained.csv"))
  if (!file.exists(disp_path)) {
    stop("Missing ", disp_path, " -- run 01_exclude_topexpressed_genes.R first")
  }
  disp_mat <- as.matrix(read.csv(disp_path, row.names = 1, check.names = FALSE))
  med_disp <- apply(disp_mat, 1, median)
  ct_ranks <- rank(med_disp, ties.method = "average")

  n_genes <- length(ct_ranks)
  low_cut <- percentage * n_genes
  high_cut <- n_genes - percentage * n_genes
  low_disp <- names(ct_ranks)[ct_ranks <= low_cut]
  high_disp <- names(ct_ranks)[ct_ranks >= high_cut]

  te_te <- read_excel(twn_s2_xlsx, sheet = sheet, col_names = FALSE)
  colnames(te_te) <- c("hub_name", "n_genes", "genes")

  df <- data.frame(genes = c(low_disp, high_disp),
                    rank_category = c(rep("low", length(low_disp)), rep("high", length(high_disp)))) %>%
    dplyr::filter(genes %in% te_te$hub_name) %>%
    dplyr::mutate(n_hub_genes = as.integer(te_te$n_genes[match(genes, te_te$hub_name)]))
  df$rank_category <- factor(df$rank_category, levels = c("low", "high"), ordered = TRUE)
  message(sprintf("[fig3c:%s] %d retained genes ranked; %d/%d low, %d/%d high found as TWN hubs",
                  cell_type, n_genes,
                  sum(df$rank_category == "low"), length(low_disp),
                  sum(df$rank_category == "high"), length(high_disp)))

  ggplot(df, aes(x = rank_category, y = n_hub_genes, fill = rank_category)) +
    geom_boxplot(width = 0.6, linewidth = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 0.6, alpha = 0.4) +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 2) +
    scale_fill_manual(values = c("low" = "steelblue", "high" = "coral")) +
    scale_x_discrete(labels = c("low" = "Low-dispersion\ngenes", "high" = "High-dispersion\ngenes")) +
    labs(x = NULL, y = "Connections per gene", title = gsub("_", " ", cell_type)) +
    theme_classic(base_family = "Arial") +
    theme(
      text = element_text(size = 6),
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 6),
      plot.title = element_text(size = 6, hjust = 0.5),
      axis.line = element_line(linewidth = 0.4),
      axis.text.x = element_text(hjust = 0.5),
      legend.position = "none",
      plot.margin = margin(2, 2, 2, 2, unit = "mm"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4)
    )
}

p_twn_fib <- ct_specific_twn_percentage_retained("Cardiac_fibroblasts", s2_sheets[30], percentage = 0.3)
p_twn_card <- ct_specific_twn_percentage_retained("Cardiomyocytes", s2_sheets[34], percentage = 0.3)

ggsave(file.path(fig3_dir, "Fig3c_fibro_card_fibro_TWN_excl_topexpressed.pdf"),
       p_twn_fib, device = cairo_pdf, width = 2.2, height = 2.54, units = "in")
ggsave(file.path(fig3_dir, "Fig3c_heart_card_TWN_excl_topexpressed.pdf"),
       p_twn_card, device = cairo_pdf, width = 2.2, height = 2.54, units = "in")
saveRDS(p_twn_fib, file.path(fig3_dir, "fig3c_fib_panel.rds"))
saveRDS(p_twn_card, file.path(fig3_dir, "fig3c_card_panel.rds"))
message("[fig3c] done")

# =============================================================================
# PART D -- Fig3D: protein-protein interaction connectivity, global Q1/Q5
# =============================================================================
ppi <- read_tsv(ppi_tsv, show_col_types = FALSE)
exp_evidence <- read_tsv(ppi_exp_tsv, show_col_types = FALSE)
ppi <- ppi %>%
  dplyr::mutate(gene_1 = exp_evidence$Gene[match(ensembl_gene_id_1, exp_evidence$Ensembl)])

# ANALYSIS_OK[data-assembly-fix]: source Rmd chunk references connectivity$Gene/
# $Freq but as.data.frame(table(x)) actually produces Var1/Freq -- a source-
# notebook state-dependent quirk (see script header); using the real column
# names here, not the source's literal (broken outside its own session) code.
connectivity <- as.data.frame(table(ppi$gene_1))
colnames(connectivity) <- c("Gene", "Freq")

df_d <- data.frame(genes = c(low_genes, high_genes),
                    rank = c(rep("1", length(low_genes)), rep("5", length(high_genes)))) %>%
  dplyr::filter(genes %in% connectivity$Gene) %>%
  dplyr::mutate(int_count = connectivity$Freq[match(genes, connectivity$Gene)])
df_d$rank <- factor(df_d$rank, levels = c("1", "5"))
message(sprintf("[fig3d] %d/%d low, %d/%d high genes found in PPI connectivity table",
                sum(df_d$rank == "1"), length(low_genes),
                sum(df_d$rank == "5"), length(high_genes)))

p_ppi <- make_low_high_boxplot(df_d, "int_count", "Protein connections per gene")
ggsave(file.path(fig3_dir, "Fig3d_protein_ints_excl_topexpressed.pdf"),
       p_ppi, device = cairo_pdf, width = 2.2, height = 2.54, units = "in")
saveRDS(p_ppi, file.path(fig3_dir, "fig3d_panel.rds"))
message("[fig3d] done")

# ------------------------------- REPORTABLE VALUES ------------------------------
wilcox_p <- function(df, y_var) wilcox.test(df[[y_var]] ~ df$rank)$p.value

write.csv(
  data.frame(
    key = c("fig3a_quintile1_p_value", "fig3a_quintile1_log2OR",
            "fig3b_twn_shared_hub_wilcox_p",
            "fig3c_cardiac_fibroblasts_wilcox_p", "fig3c_cardiomyocytes_wilcox_p",
            "fig3d_ppi_connectivity_wilcox_p"),
    value = c(plot_df$p.value[plot_df$quintile == levels(plot_df$quintile)[1]],
              plot_df$log2OR[plot_df$quintile == levels(plot_df$quintile)[1]],
              wilcox_p(df_b, "hub_rank"),
              wilcox.test(n_hub_genes ~ rank_category, data = p_twn_fib$data)$p.value,
              wilcox.test(n_hub_genes ~ rank_category, data = p_twn_card$data)$p.value,
              wilcox_p(df_d, "int_count")),
    provenance = "figure_reproduction/scripts/27_fig3_panels_excl_topexpressed.R"
  ),
  file.path(out_dir, "reportable_values_27.csv"),
  row.names = FALSE
)

message("[done] Figure 3 panels written to ", fig3_dir)
