#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(readr)
  library(variancePartition)  # topTable (re-exported from limma)
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
cluster_dir <- file.path(repo_dir, "cellranger_cluster-mode_trial")
fig_repro_dir <- file.path(repo_dir, "figure_reproduction")
out_root <- file.path(fig_repro_dir, "outputs")
out_dir <- file.path(out_root, "regulatory_network_excl_cellstate")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dispersion_csv <- file.path(fig_repro_dir, "data", "res_var_excl_cluster24.csv")
retained_matrix_csv <- file.path(out_root, "retained_gene_matrix.csv")
filt_mat_dir <- file.path(out_root, "filtered_matrices")
ranked_dir <- file.path(out_root, "ranked")
excluded_csv <- file.path(out_root, "dd_refit", "excl_cellstate_genes", "excluded_genes.csv")
baseline_values_csvs <- file.path(out_root, c("reportable_values_27.csv", "reportable_values_29.csv"))

# External inputs -- identical paths to 03 / 27 / 29.
# ANALYSIS_OK[external-path]: read-only path to Erik's frozen eQTL results, identical to 03 Part A
eqtl_dir <- file.path(
  "/project/gilad/emcintire/chemo_project_2023/pipeline/map_eqtls",
  "archive", "no_cvrt_condox_only", "eqtl_results", "Control", "%s")
eqtl_file_tmpl <- "steadyState_%s_Control_AllChroms.txt"
ext_dir <- file.path(cluster_dir, "analysis", "datasets", "external_data")
twn_s2_xlsx <- file.path(ext_dir, "saha_coexp_net_TWN_2017", "Supplemental_Data_S2.xlsx")
twn_s4_xlsx <- file.path(ext_dir, "saha_coexp_net_TWN_2017", "Supplemental_Data_S4.xlsx")
ppi_tsv <- file.path(ext_dir, "protein_interaction_data", "interaction_consensus.tsv.zip")
ppi_exp_tsv <- file.path(ext_dir, "protein_interaction_data", "has_data_Tissue.tsv")
# Script 04's DE fit (Epithelial sub-cluster 24 excluded, union-retained genes).
dream_de_rds <- file.path(out_root, "dd_refit", "dream_fit_mean_retained.rds")
epd_bed <- file.path(ext_dir, "epd", "human_epdnew_3sa0y.bed")
fantom_tss_txt <- "/project2/gilad/wenhe/cox_baboon_diet/combined_analysis/gene_annotation/fantom/short_descrip_split_gene_only.txt"
abc_txt <- "/project2/gilad/wenhe/cox_baboon_diet/combined_analysis/gene_annotation/ABC_prediction_score005.txt"
mane_summary_txt <- file.path(repo_dir, "data", "raw", "mane_v1_5", "MANE.GRCh38.v1.5.summary.txt.gz")
mane_gff_gz <- file.path(repo_dir, "data", "raw", "mane_v1_5", "MANE.GRCh38.v1.5.ensembl_genomic.gff.gz")

# Parameters -- identical to 02 / 03 Part A / 27 / 29.
ribo_pattern <- "^RP"
n_quantiles <- 5
tolerance_frac <- 0.15        # 03 Part A
egene_pvalue_thresh <- 0.05   # 03 Part A
egene_min_count <- 49         # 03 Part A
ct_percentage <- 0.3          # 27 Part C (bottom/top 30%)
de_adjp_thresh <- 0.05        # 29 Part A
n_excluded_expected <- 146L   # 97 cell cycle + 44 markers + 6 chamber, TBX5 in two sets

# ANALYSIS_OK[experimental-design]: the 3 pooled iPSC lines (column prefixes), not an exclusion list; same as 02/03
individuals <- c("NA19099", "NA19093", "NA18522")
cell_types <- c(
  "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
  "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"
)

source(file.path(fig_repro_dir, "scripts", "utils.R"))  # ct_columns()

for (f in c(dispersion_csv, retained_matrix_csv, excluded_csv, baseline_values_csvs,
            twn_s2_xlsx, twn_s4_xlsx, ppi_tsv, ppi_exp_tsv, dream_de_rds, epd_bed,
            fantom_tss_txt, abc_txt, mane_summary_txt, mane_gff_gz)) {
  if (!file.exists(f)) stop("Missing input: ", f)
}

# ------------------------------- GENE UNIVERSE --------------------------------
excl_tab <- read.csv(excluded_csv, stringsAsFactors = FALSE)
excluded_genes <- unique(excl_tab$gene)
if (length(excluded_genes) != n_excluded_expected) {
  stop("excluded_genes.csv has ", length(excluded_genes), " genes; expected ", n_excluded_expected)
}

res_var_mat <- as.matrix(read.csv(dispersion_csv, row.names = 1, check.names = FALSE))
if (anyNA(res_var_mat)) stop("Dispersion matrix contains NA values")
# ANALYSIS_OK[filtering]: drop ^RP ribosomal genes, identical to 02/03; count logged below
res_var_mat <- res_var_mat[!grepl(ribo_pattern, rownames(res_var_mat)), , drop = FALSE]
all_genes <- rownames(res_var_mat)

retained_ind <- as.matrix(read.csv(retained_matrix_csv, row.names = 1, check.names = FALSE))
if (!setequal(colnames(retained_ind), cell_types)) stop("retained_gene_matrix.csv columns != cell types")
# ANALYSIS_OK[column-reorder]: canonical cell-type column order only; no rows dropped
retained_ind <- retained_ind[, cell_types, drop = FALSE]
if (!all(all_genes %in% rownames(retained_ind))) stop("retained_gene_matrix.csv missing non-ribo genes")

union_retained <- all_genes[rowSums(retained_ind[all_genes, , drop = FALSE] == 1) == length(cell_types)]
# ANALYSIS_OK[filtering]: the requested cell-state exclusion; removed count logged, genes listed in 39's excluded_genes.csv
union_excl <- setdiff(union_retained, excluded_genes)
message(sprintf("[universe] %d non-ribo genes; union-retained %d -> %d after removing %d cell-state genes",
                length(all_genes), length(union_retained), length(union_excl),
                length(union_retained) - length(union_excl)))

# ---------------------- CONSISTENTLY-LOW / -HIGH (02 logic) -------------------
quint <- matrix(NA_integer_, nrow = length(all_genes), ncol = length(cell_types),
                dimnames = list(all_genes, cell_types))
for (ct in cell_types) {
  # ANALYSIS_OK[filtering]: per-CT retained genes (01) minus the cell-state genes; count logged
  ct_genes <- setdiff(all_genes[retained_ind[all_genes, ct] == 1], excluded_genes)
  med <- apply(res_var_mat[ct_genes, ct_columns(ct), drop = FALSE], 1, median)
  r <- rank(med, ties.method = "average")
  quint[ct_genes, ct] <- cut(r, breaks = quantile(r, probs = seq(0, 1, length.out = n_quantiles + 1)),
                             labels = FALSE, include.lowest = TRUE)
  message(sprintf("[rank] %-25s %d genes ranked", ct, length(ct_genes)))
}
low_genes <- rownames(quint)[rowSums(quint == 1, na.rm = TRUE) == length(cell_types)]
high_genes <- rownames(quint)[rowSums(quint == n_quantiles, na.rm = TRUE) == length(cell_types)]
stopifnot(!any(c(low_genes, high_genes) %in% excluded_genes))
writeLines(low_genes, file.path(out_dir, "consistently_low_excl_cellstate.txt"))
writeLines(high_genes, file.path(out_dir, "consistently_high_excl_cellstate.txt"))

base_low <- readLines(file.path(ranked_dir, "consistently_low_retained.txt"))
base_high <- readLines(file.path(ranked_dir, "consistently_high_retained.txt"))
set_cmp <- data.frame(
  gene_set = c("consistently_low", "consistently_high"),
  n_baseline = c(length(base_low), length(base_high)),
  n_excl_cellstate = c(length(low_genes), length(high_genes)),
  n_removed_cellstate = c(sum(base_low %in% excluded_genes), sum(base_high %in% excluded_genes)),
  n_other_lost = c(length(setdiff(setdiff(base_low, excluded_genes), low_genes)),
                   length(setdiff(setdiff(base_high, excluded_genes), high_genes))),
  n_gained = c(length(setdiff(low_genes, base_low)), length(setdiff(high_genes, base_high))))
write.csv(set_cmp, file.path(out_dir, "consistent_sets_excl_cellstate.csv"), row.names = FALSE)
print(set_cmp)

# Low/high data frame for one per-gene feature (29's low_high_df, rank 1 = low, 5 = high).
low_high_df <- function(var, values_by_gene) {
  df <- data.frame(genes = c(low_genes, high_genes),
                   rank = c(rep("1", length(low_genes)), rep("5", length(high_genes))))
  df[[var]] <- unname(values_by_gene[df$genes])
  # ANALYSIS_OK[filtering]: genes without a value in the external feature table are dropped (as in 27/29); found counts logged by log_found()
  df <- df[!is.na(df[[var]]), , drop = FALSE]
  df$rank <- factor(df$rank, levels = c("1", "5"))
  df
}
log_found <- function(tag, df) {
  message(sprintf("[%s] %d/%d low, %d/%d high genes have a value", tag,
                  sum(df$rank == "1"), length(low_genes), sum(df$rank == "5"), length(high_genes)))
}

# ------------------------- Fig3A: eGene quintile Fisher -----------------------
# ANALYSIS_OK[within-file-helper]: called below in this script only; local by design
load_egenes <- function(ct) {
  path <- file.path(sprintf(eqtl_dir, ct), sprintf(eqtl_file_tmpl, ct))
  if (!file.exists(path)) stop("Missing eQTL results file: ", path)
  tab <- readr::read_table(path, show_col_types = FALSE)
  if (!all(c("gene", "bh_adj_pvalue") %in% colnames(tab))) stop("eQTL file missing columns: ", path)
  unique(tab$gene[tab$bh_adj_pvalue <= egene_pvalue_thresh])
}
egenes_by_ct <- list()
for (ct in cell_types) {
  eg <- load_egenes(ct)
  message(sprintf("[egene] %-25s %d eGenes", ct, length(eg)))
  if (length(eg) >= egene_min_count) egenes_by_ct[[ct]] <- eg
}
all_egenes <- unique(unlist(egenes_by_ct, use.names = FALSE))

ct_ranks <- sapply(cell_types, function(ct) {
  rank(apply(res_var_mat[union_excl, ct_columns(ct), drop = FALSE], 1, median), ties.method = "average")
})
avg_rank <- apply(ct_ranks, 1, median)
# ANALYSIS_OK[filtering]: 03 Part A's rank-consistency tolerance filter (15% of universe size); kept count logged
within_tol <- apply(abs(ct_ranks - avg_rank) <= tolerance_frac * nrow(ct_ranks), 1, all)
filtered_ranks <- ct_ranks[within_tol, , drop = FALSE]
message(sprintf("[fig3a] %d genes within rank tolerance", nrow(filtered_ranks)))
global_rank <- rank(apply(filtered_ranks, 1, median), ties.method = "average")
global_quintile <- setNames(
  cut(global_rank, breaks = quantile(global_rank, probs = seq(0, 1, length.out = n_quantiles + 1)),
      labels = FALSE, include.lowest = TRUE),
  rownames(filtered_ranks))
background <- names(global_quintile)
egenes_bg <- intersect(all_egenes, background)
fig3a_df <- do.call(rbind, lapply(seq_len(n_quantiles), function(i) {
  q_genes <- background[global_quintile == i]
  a <- length(intersect(egenes_bg, q_genes))
  ft <- fisher.test(matrix(c(a, length(egenes_bg) - a, length(q_genes) - a,
                             length(background) - length(union(egenes_bg, q_genes))),
                           nrow = 2, byrow = TRUE))
  data.frame(quintile = i, egene_count = a, p.value = ft$p.value, OR = unname(ft$estimate),
             log2OR = log2(unname(ft$estimate)),
             log2_ci_lower = log2(ft$conf.int[1]), log2_ci_upper = log2(ft$conf.int[2]))
}))
fig3a_df$padj <- p.adjust(fig3a_df$p.value, method = "BH")
print(fig3a_df)

# ------------------------- Fig3B: TWN shared-hub rank -------------------------
s4_sheets <- excel_sheets(twn_s4_xlsx)
stopifnot("Supplemental_Data_S4.xlsx sheet 10 layout changed" = grepl("^Group3", s4_sheets[10]))
shared_te_te <- read_excel(twn_s4_xlsx, sheet = s4_sheets[10], col_names = TRUE)
fig3b_df <- low_high_df("hub_rank", setNames(as.integer(shared_te_te$rank), shared_te_te$hub))
log_found("fig3b", fig3b_df)

# ------------------ Fig3C: cell-type-specific TWN connections ------------------
s2_sheets <- excel_sheets(twn_s2_xlsx)
ct_twn_df <- function(ct, sheet) {
  disp_mat <- as.matrix(read.csv(file.path(filt_mat_dir, paste0(ct, "_dispersion_retained.csv")),
                                 row.names = 1, check.names = FALSE))
  # ANALYSIS_OK[filtering]: the requested cell-state exclusion before per-CT ranking; count logged
  disp_mat <- disp_mat[!rownames(disp_mat) %in% excluded_genes, , drop = FALSE]
  r <- rank(apply(disp_mat, 1, median), ties.method = "average")
  n <- length(r)
  low_disp <- names(r)[r <= ct_percentage * n]
  high_disp <- names(r)[r >= n - ct_percentage * n]
  te_te <- read_excel(twn_s2_xlsx, sheet = sheet, col_names = FALSE)
  colnames(te_te) <- c("hub_name", "n_genes", "genes")
  df <- data.frame(genes = c(low_disp, high_disp),
                   rank_category = c(rep("low", length(low_disp)), rep("high", length(high_disp)))) %>%
    dplyr::filter(genes %in% te_te$hub_name) %>%
    dplyr::mutate(n_hub_genes = as.integer(te_te$n_genes[match(genes, te_te$hub_name)]))
  df$rank_category <- factor(df$rank_category, levels = c("low", "high"), ordered = TRUE)
  message(sprintf("[fig3c:%s] %d genes ranked; %d/%d low, %d/%d high are TWN hubs", ct, n,
                  sum(df$rank_category == "low"), length(low_disp),
                  sum(df$rank_category == "high"), length(high_disp)))
  df
}
fig3c_fib_df <- ct_twn_df("Cardiac_fibroblasts", s2_sheets[30])   # same sheets as 27
fig3c_card_df <- ct_twn_df("Cardiomyocytes", s2_sheets[34])

# ------------------------- Fig3D: PPI connectivity ----------------------------
ppi <- read_tsv(ppi_tsv, show_col_types = FALSE)
exp_evidence <- read_tsv(ppi_exp_tsv, show_col_types = FALSE)
ppi_gene <- exp_evidence$Gene[match(ppi$ensembl_gene_id_1, exp_evidence$Ensembl)]
connectivity <- table(ppi_gene)
fig3d_df <- low_high_df("int_count", setNames(as.integer(connectivity), names(connectivity)))
log_found("fig3d", fig3d_df)

# ------------------------- Fig4A: DE count ------------------------------------
de_dream_obj <- readRDS(dream_de_rds)
contrast_coefs <- colnames(de_dream_obj$coefficients)[1:15]
de_calls <- sapply(contrast_coefs, function(k) {
  topTable(de_dream_obj, coef = k, number = Inf, sort.by = "none")$adj.P.Val <= de_adjp_thresh
})
fig4a_df <- low_high_df("DE_count", setNames(rowSums(de_calls), rownames(de_dream_obj$coefficients)))
log_found("fig4a", fig4a_df)

# ------------------------- Fig4B: TATA Fisher ---------------------------------
epd_tata <- read.table(epd_bed, sep = "\t", header = FALSE)
tata_genes <- intersect(unique(sub("_[0-9]+$", "", epd_tata$V4)), union_excl)
# ANALYSIS_OK[within-file-helper]: called below in this script only; local by design
tata_fisher <- function(target) {
  target <- intersect(target, union_excl)
  a <- length(intersect(tata_genes, target))
  ft <- fisher.test(matrix(c(a, length(tata_genes) - a, length(target) - a,
                             length(union_excl) - length(union(tata_genes, target))),
                           nrow = 2, byrow = TRUE))
  data.frame(odd_ratio = unname(ft$estimate), p_val = ft$p.value,
             log2_or = log2(unname(ft$estimate)),
             log2_ci_lower = log2(ft$conf.int[1]), log2_ci_upper = log2(ft$conf.int[2]))
}
fig4b_df <- cbind(gene_set = factor(c("low_disp", "high_disp"), levels = c("low_disp", "high_disp")),
                  rbind(tata_fisher(low_genes), tata_fisher(high_genes)))
message(sprintf("[fig4b] background %d genes; low OR=%.3f p=%.3g; high OR=%.3f p=%.3g",
                length(union_excl), fig4b_df$odd_ratio[1], fig4b_df$p_val[1],
                fig4b_df$odd_ratio[2], fig4b_df$p_val[2]))

# ------------------------- Fig4C: FANTOM TSS count ----------------------------
tss_tab <- table(read.csv(fantom_tss_txt, header = FALSE)$V1)
fig4c_df <- low_high_df("tss_count", setNames(as.integer(tss_tab), names(tss_tab)))
log_found("fig4c", fig4c_df)

# ------------------------- Fig4D: ABC enhancer length -------------------------
abc_pred <- read_table(abc_txt, show_col_types = FALSE)
fig4d_df <- low_high_df("enh_length", setNames(abc_pred$avg_length_across_celltype, abc_pred$TargetGene))
log_found("fig4d", fig4d_df)

# ------------------------- Fig4E: MANE Select intron count --------------------
# ANALYSIS_OK[within-file-helper]: called below in this script only; local by design
strip_version <- function(x) sub("\\..*$", "", x)
mane_summary <- read.delim(mane_summary_txt, stringsAsFactors = FALSE, comment.char = "")
colnames(mane_summary)[1] <- sub("^#", "", colnames(mane_summary)[1])
mane_select <- mane_summary[mane_summary$MANE_status == "MANE Select", , drop = FALSE]
mane_gff <- read.delim(mane_gff_gz, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
exon_attr <- mane_gff$V9[mane_gff$V3 == "exon"]
exon_tx <- strip_version(sub("^transcript_id=", "", regmatches(exon_attr, regexpr("transcript_id=[^;]+", exon_attr))))
stopifnot("failed to parse transcript_id from some GFF exon rows" = length(exon_tx) == length(exon_attr))
exon_counts <- table(exon_tx)
n_introns <- as.integer(exon_counts[strip_version(mane_select$Ensembl_nuc)]) - 1L
stopifnot("some MANE Select transcripts have no exon rows" = !anyNA(n_introns))
fig4e_df <- low_high_df("intron_count", setNames(n_introns, mane_select$symbol))
log_found("fig4e", fig4e_df)

# ------------------------------- SAVE PANEL DATA -------------------------------
panel_data <- list(fig3a = fig3a_df, fig3b = fig3b_df, fig3c_fib = fig3c_fib_df,
                   fig3c_card = fig3c_card_df, fig3d = fig3d_df, fig4a = fig4a_df,
                   fig4b = fig4b_df, fig4c = fig4c_df, fig4d = fig4d_df, fig4e = fig4e_df)
saveRDS(panel_data, file.path(out_dir, "panel_data.rds"))

# ------------------------- BASELINE vs EXCLUDED STATS ------------------------
# ANALYSIS_OK[within-file-helper]: called below in this script only; local by design
wp <- function(df, y, g) wilcox.test(df[[y]] ~ df[[g]])$p.value
q1 <- fig3a_df$quintile == 1
new_vals <- c(
  fig3a_quintile1_p_value = fig3a_df$p.value[q1],
  fig3a_quintile1_log2OR = fig3a_df$log2OR[q1],
  fig3b_twn_shared_hub_wilcox_p = wp(fig3b_df, "hub_rank", "rank"),
  fig3c_cardiac_fibroblasts_wilcox_p = wp(fig3c_fib_df, "n_hub_genes", "rank_category"),
  fig3c_cardiomyocytes_wilcox_p = wp(fig3c_card_df, "n_hub_genes", "rank_category"),
  fig3d_ppi_connectivity_wilcox_p = wp(fig3d_df, "int_count", "rank"),
  fig4a_de_count_wilcox_p = wp(fig4a_df, "DE_count", "rank"),
  fig4b_tata_low_log2OR = fig4b_df$log2_or[1],
  fig4b_tata_low_p_value = fig4b_df$p_val[1],
  fig4b_tata_high_log2OR = fig4b_df$log2_or[2],
  fig4b_tata_high_p_value = fig4b_df$p_val[2],
  fig4c_tss_count_wilcox_p = wp(fig4c_df, "tss_count", "rank"),
  fig4d_enh_length_wilcox_p = wp(fig4d_df, "enh_length", "rank"),
  fig4e_intron_count_wilcox_p = wp(fig4e_df, "intron_count", "rank"))

baseline <- do.call(rbind, lapply(baseline_values_csvs, read.csv, stringsAsFactors = FALSE))
if (!setequal(baseline$key, names(new_vals))) stop("Baseline reportable-value keys differ from this script's keys")
stats_cmp <- data.frame(key = names(new_vals),
                        baseline = baseline$value[match(names(new_vals), baseline$key)],
                        excl_cellstate = unname(new_vals))
write.csv(stats_cmp, file.path(out_dir, "stats_baseline_vs_excl_cellstate.csv"), row.names = FALSE)
print(stats_cmp)

write.csv(
  rbind(
    data.frame(key = c("n_excluded_genes", "n_union_retained", "n_union_excl_cellstate",
                       "n_low_excl_cellstate", "n_high_excl_cellstate"),
               value = c(length(excluded_genes), length(union_retained), length(union_excl),
                         length(low_genes), length(high_genes))),
    data.frame(key = paste0(names(new_vals), "_excl_cellstate"), value = unname(new_vals))
  ) %>% dplyr::mutate(provenance = "figure_reproduction/scripts/regulatory_network_excl_cellstate/01_panels_excl_cellstate.R"),
  file.path(out_dir, "reportable_values_rnx.csv"), row.names = FALSE)

message("[done] wrote panel data and stats to ", out_dir)
