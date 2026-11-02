#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(variancePartition)  # topTable (re-exported from limma)
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
cluster_dir <- file.path(repo_dir, "cellranger_cluster-mode_trial")
analysis_dir <- file.path(repo_dir, "figure_reproduction")
out_dir <- file.path(analysis_dir, "outputs")
ranked_dir <- file.path(out_dir, "ranked")
fig4_dir <- file.path(out_dir, "fig4_replication")
dir.create(fig4_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(analysis_dir, "scripts", "utils.R"))

cell_types <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)

dream_de_rds <- file.path(out_dir, "dd_refit", "dream_fit_mean_retained.rds")
epd_bed <- file.path(cluster_dir, "analysis", "datasets", "external_data",
                      "epd", "human_epdnew_3sa0y.bed")
fantom_tss_txt <- "/project2/gilad/wenhe/cox_baboon_diet/combined_analysis/gene_annotation/fantom/short_descrip_split_gene_only.txt"
abc_txt <- "/project2/gilad/wenhe/cox_baboon_diet/combined_analysis/gene_annotation/ABC_prediction_score005.txt"
mane_summary_txt <- file.path(repo_dir, "data", "raw", "mane_v1_5", "MANE.GRCh38.v1.5.summary.txt.gz")
mane_gff_gz <- file.path(repo_dir, "data", "raw", "mane_v1_5", "MANE.GRCh38.v1.5.ensembl_genomic.gff.gz")
retained_matrix_csv <- file.path(out_dir, "retained_gene_matrix.csv")

low_genes <- readLines(file.path(ranked_dir, "consistently_low_retained.txt"))
high_genes <- readLines(file.path(ranked_dir, "consistently_high_retained.txt"))
message(sprintf("[global] consistently-low: %d genes, consistently-high: %d genes",
                length(low_genes), length(high_genes)))

low_high_df <- function(matched_var, values_by_gene) {
  df <- data.frame(genes = c(low_genes, high_genes),
                    rank = c(rep("1", length(low_genes)), rep("5", length(high_genes))))
  df[[matched_var]] <- values_by_gene[df$genes]
  df <- df[!is.na(df[[matched_var]]), , drop = FALSE]
  df$rank <- factor(df$rank, levels = c("1", "5"))
  df
}

# =============================================================================
# PART A -- Fig4A: DE count across 15 contrasts, global low vs high
# =============================================================================
de_dream_obj <- readRDS(dream_de_rds)
fc_thresh <- 0
contrast_coefs <- colnames(de_dream_obj$coefficients)[1:15]

contrast_DE <- as.data.frame(matrix(
  NA, nrow = nrow(de_dream_obj$coefficients), ncol = length(contrast_coefs),
  dimnames = list(rownames(de_dream_obj$coefficients), contrast_coefs)))

for (contrast in contrast_coefs) {
  res_sum_table <- topTable(de_dream_obj, coef = contrast, number = Inf, sort.by = "none")
  contrast_DE[[contrast]] <- res_sum_table$adj.P.Val <= 0.05 & abs(res_sum_table$logFC) >= fc_thresh
}

de_count_by_gene <- setNames(rowSums(contrast_DE), rownames(contrast_DE))
df_a <- low_high_df("DE_count", de_count_by_gene)
message(sprintf("[fig4a] %d/%d low, %d/%d high genes found in the dream DE fit",
                sum(df_a$rank == "1"), length(low_genes),
                sum(df_a$rank == "5"), length(high_genes)))

p_de <- make_low_high_boxplot(df_a, "DE_count", "DE count per gene")
ggsave(file.path(fig4_dir, "Fig4a_DE_count_low_vs_high_excl_topexpressed.pdf"),
       p_de, device = cairo_pdf, width = 2.2, height = 2.54, units = "in")
saveRDS(p_de, file.path(fig4_dir, "fig4a_panel.rds"))
message("[fig4a] done")

# =============================================================================
# PART B -- Fig4B: TATA (EPD) Fisher enrichment, global low vs high
# =============================================================================
retained_ind <- as.matrix(read.csv(retained_matrix_csv, row.names = 1, check.names = FALSE))
retained_ind <- retained_ind[, cell_types, drop = FALSE]
union_retained_genes <- rownames(retained_ind)[rowSums(retained_ind == 1) == length(cell_types)]
message(sprintf("[fig4b] union-retained background: %d genes", length(union_retained_genes)))

epd_tata <- read.table(epd_bed, sep = "\t", header = FALSE)
colnames(epd_tata) <- c("chr", "start", "end", "promoter_id", "score", "strand")
epd_tata <- epd_tata %>%
  dplyr::mutate(gene = sub("_[0-9]+$", "", promoter_id))
tata_genes <- intersect(unique(epd_tata$gene), union_retained_genes)

tata_enrich_fisher <- function(target_genes) {
  target_genes <- intersect(target_genes, union_retained_genes)
  a <- length(intersect(tata_genes, target_genes))
  b <- length(setdiff(tata_genes, target_genes))
  c <- length(setdiff(target_genes, tata_genes))
  d <- length(setdiff(union_retained_genes, union(tata_genes, target_genes)))
  contingency <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                        dimnames = list(TATA = c("Yes", "No"), Target = c("In_Set", "Not_In_Set")))
  fisher.test(contingency)
}

low_ft <- tata_enrich_fisher(low_genes)
high_ft <- tata_enrich_fisher(high_genes)

plot_df_b <- data.frame(
  gene_set = factor(c("low_disp", "high_disp"), levels = c("low_disp", "high_disp")),
  odd_ratio = c(as.numeric(low_ft$estimate), as.numeric(high_ft$estimate)),
  ci_lower = c(low_ft$conf.int[1], high_ft$conf.int[1]),
  ci_upper = c(low_ft$conf.int[2], high_ft$conf.int[2]),
  p_val = c(low_ft$p.value, high_ft$p.value)
) %>%
  dplyr::mutate(
    p_label = ifelse(p_val < 2e-16, "p = <2e-16 ", paste0("p = ", formatC(p_val, format = "e", digits = 2))),
    log2_or = log2(odd_ratio),
    log2_ci_lower = log2(ci_lower),
    log2_ci_upper = log2(ci_upper)
  )
message(sprintf("[fig4b] low OR=%.3f p=%.3g; high OR=%.3f p=%.3g",
                plot_df_b$odd_ratio[1], plot_df_b$p_val[1],
                plot_df_b$odd_ratio[2], plot_df_b$p_val[2]))

p_tata <- ggplot(plot_df_b, aes(x = gene_set, y = log2_or)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey40") +
  geom_errorbar(aes(ymin = log2_ci_lower, ymax = log2_ci_upper), width = 0.15, linewidth = 0.4) +
  geom_point(size = 1.6) +
  geom_text(aes(y = log2_ci_upper, label = p_label), vjust = -0.6, size = 2, family = "Arial") +
  scale_x_discrete(labels = c("low_disp" = "Low-dispersion\ngenes", "high_disp" = "High-dispersion\ngenes")) +
  labs(x = NULL, y = expression(log[2]("Odds ratio"))) +
  theme_classic(base_family = "Arial") +
  theme(
    text = element_text(size = 6),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.line = element_line(linewidth = 0.4),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2, unit = "mm")
  )
ggsave(file.path(fig4_dir, "Fig4b_tata_enrich_low_vs_high_excl_topexpressed.pdf"),
       p_tata, device = cairo_pdf, width = 2.2, height = 2.54, units = "in")
saveRDS(p_tata, file.path(fig4_dir, "fig4b_panel.rds"))
message("[fig4b] done")

# =============================================================================
# PART C -- Fig4C: FANTOM TSS count, global low vs high
# =============================================================================
tss_usage <- read.csv(fantom_tss_txt, header = FALSE)
tss_df <- as.data.frame(table(tss_usage$V1))
tss_by_gene <- setNames(tss_df$Freq, as.character(tss_df$Var1))

df_c <- low_high_df("tss_count", tss_by_gene)
message(sprintf("[fig4c] %d/%d low, %d/%d high genes found in FANTOM TSS table",
                sum(df_c$rank == "1"), length(low_genes),
                sum(df_c$rank == "5"), length(high_genes)))

p_tss <- make_low_high_boxplot(df_c, "tss_count", "TSS count")
ggsave(file.path(fig4_dir, "Fig4c_TSS_count_low_vs_high_excl_topexpressed.pdf"),
       p_tss, device = cairo_pdf, width = 2.2, height = 2.54, units = "in")
saveRDS(p_tss, file.path(fig4_dir, "fig4c_panel.rds"))
message("[fig4c] done")

# =============================================================================
# PART D -- Fig4D: ABC accumulated enhancer length, global low vs high
# =============================================================================
abc_pred <- read_table(abc_txt, show_col_types = FALSE)
enh_by_gene <- setNames(abc_pred$avg_length_across_celltype, abc_pred$TargetGene)

df_d <- low_high_df("enh_length", enh_by_gene)
message(sprintf("[fig4d] %d/%d low, %d/%d high genes found in ABC enhancer table",
                sum(df_d$rank == "1"), length(low_genes),
                sum(df_d$rank == "5"), length(high_genes)))

p_enh <- make_low_high_boxplot(df_d, "enh_length", "Accum. enhancer length\n(per tissue/cell type)")
ggsave(file.path(fig4_dir, "Fig4d_enh_length_low_vs_high_excl_topexpressed.pdf"),
       p_enh, device = cairo_pdf, width = 2.2, height = 2.54, units = "in")
saveRDS(p_enh, file.path(fig4_dir, "fig4d_panel.rds"))
message("[fig4d] done")

# =============================================================================
# PART E -- Fig4E: MANE Select intron count, global low vs high
# =============================================================================
strip_version <- function(x) sub("\\..*$", "", x)
mane_summary <- read.delim(mane_summary_txt, stringsAsFactors = FALSE, comment.char = "")
colnames(mane_summary)[1] <- sub("^#", "", colnames(mane_summary)[1])
mane_select <- mane_summary[mane_summary$MANE_status == "MANE Select", , drop = FALSE]
mane_select$ensembl_transcript_id <- strip_version(mane_select$Ensembl_nuc)

mane_gff <- read.delim(mane_gff_gz, header = FALSE, comment.char = "#", stringsAsFactors = FALSE,
                       col.names = c("seqid", "source", "feature", "start", "end",
                                     "score", "strand", "frame", "attributes"))
exon_gff <- mane_gff[mane_gff$feature == "exon", , drop = FALSE]
extract_attr <- function(attributes, key) {
  m <- regmatches(attributes, regexpr(paste0(key, "=[^;]+"), attributes))
  ifelse(nchar(m) == 0, NA_character_, sub(paste0("^", key, "="), "", m))
}
exon_gff$ensembl_transcript_id <- strip_version(extract_attr(exon_gff$attributes, "transcript_id"))
stopifnot("failed to parse transcript_id from some GFF exon rows" =
            !any(is.na(exon_gff$ensembl_transcript_id)))
exon_counts_tbl <- as.data.frame(table(exon_gff$ensembl_transcript_id), stringsAsFactors = FALSE)
colnames(exon_counts_tbl) <- c("ensembl_transcript_id", "num_exons")
mane_select$num_introns <- exon_counts_tbl$num_exons[
  match(mane_select$ensembl_transcript_id, exon_counts_tbl$ensembl_transcript_id)] - 1L
stopifnot("some MANE Select transcripts have no exon rows in the GFF" =
            !any(is.na(mane_select$num_introns)))
intron_by_gene <- setNames(mane_select$num_introns, mane_select$symbol)
message(sprintf("[fig4e] MANE Select intron-count table: %d genes (genome-wide)", length(intron_by_gene)))

df_e <- low_high_df("intron_count", intron_by_gene)
message(sprintf("[fig4e] %d/%d low, %d/%d high genes found in MANE Select",
                sum(df_e$rank == "1"), length(low_genes),
                sum(df_e$rank == "5"), length(high_genes)))

p_intron <- make_low_high_boxplot(df_e, "intron_count", "Intron count")
ggsave(file.path(fig4_dir, "Fig4e_intron_count_low_vs_high_excl_topexpressed.pdf"),
       p_intron, device = cairo_pdf, width = 2.2, height = 2.54, units = "in")
saveRDS(p_intron, file.path(fig4_dir, "fig4e_panel.rds"))
message("[fig4e] done")

# ------------------------------- REPORTABLE VALUES ------------------------------
wilcox_p <- function(df, y_var) wilcox.test(df[[y_var]] ~ df$rank)$p.value

write.csv(
  data.frame(
    key = c("fig4a_de_count_wilcox_p",
            "fig4b_tata_low_log2OR", "fig4b_tata_low_p_value",
            "fig4b_tata_high_log2OR", "fig4b_tata_high_p_value",
            "fig4c_tss_count_wilcox_p",
            "fig4d_enh_length_wilcox_p",
            "fig4e_intron_count_wilcox_p"),
    value = c(wilcox_p(df_a, "DE_count"),
              plot_df_b$log2_or[1], plot_df_b$p_val[1],
              plot_df_b$log2_or[2], plot_df_b$p_val[2],
              wilcox_p(df_c, "tss_count"),
              wilcox_p(df_d, "enh_length"),
              wilcox_p(df_e, "intron_count")),
    provenance = "figure_reproduction/scripts/29_fig4_panels_excl_topexpressed.R"
  ),
  file.path(out_dir, "reportable_values_29.csv"),
  row.names = FALSE
)

message("[done] Figure 4 panels written to ", fig4_dir)
