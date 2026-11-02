#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args       <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
source(file.path(script_dir, "config.R"))
panel_dir <- file.path(out_dir, "fig3_fig4_panels")
comb_dir  <- ensure_dir(file.path(out_dir, "combined"))

# External inputs -- identical paths to regulatory_network_excl_cellstate/01.
# ANALYSIS_OK[external-path]: read-only path to Erik's frozen eQTL results, identical to 03 / rnx 01
eqtl_tmpl <- eqtl_result_template
fantom_tss_txt <- "/project2/gilad/wenhe/cox_baboon_diet/combined_analysis/gene_annotation/fantom/short_descrip_split_gene_only.txt"
abc_txt <- "/project2/gilad/wenhe/cox_baboon_diet/combined_analysis/gene_annotation/ABC_prediction_score005.txt"
mane_summary_txt <- file.path(repo_dir, "data", "raw", "mane_v1_5", "MANE.GRCh38.v1.5.summary.txt.gz")
mane_gff_gz <- file.path(repo_dir, "data", "raw", "mane_v1_5", "MANE.GRCh38.v1.5.ensembl_genomic.gff.gz")
for (f in c(fantom_tss_txt, abc_txt, mane_summary_txt, mane_gff_gz)) {
  assert_that(file.exists(f), paste("missing input:", f))
}

gene_sets <- readRDS(file.path(out_dir, "gene_sets", "dispersion_gene_sets.rds"))
mm        <- readRDS(file.path(out_dir, "mean_matched", "matched_gene_sets.rds"))
all_genes_bg <- rownames(gene_sets$var_quintile_dfs)  # ribo-filtered, -6-retained universe, unmatched
primary <- mm$matched_by_seed[[as.character(mm$match_seed_primary)]]
log_step(sprintf("primary seed %d matched sizes: low=%d neither=%d high=%d", mm$match_seed_primary,
                 length(primary$first), length(primary$never), length(primary$fifth)))

# =============================================================================
# Fig3A: eGene enrichment per matched group vs full background
# =============================================================================
egenes <- character(0)
for (ct in cell_types) {
  path <- sprintf(eqtl_tmpl, ct, ct)
  assert_that(file.exists(path), paste("missing eQTL result:", path))
  tab <- read_table(path, show_col_types = FALSE)
  assert_that(all(c("gene", "bh_adj_pvalue") %in% colnames(tab)), paste("eQTL file missing columns:", path))
  eg <- unique(tab$gene[tab$bh_adj_pvalue <= egene_padj])
  log_step(sprintf("%s eGenes", ct), length(eg))
  if (length(eg) >= egene_min_count) egenes <- c(egenes, eg)
}
egenes_bg <- intersect(unique(egenes), all_genes_bg)
log_step("eGenes in the ribo-filtered background", length(egenes_bg))

# ANALYSIS_OK[within-file-helper]: called below for each group and seed; local by design
egene_fisher <- function(target) {
  a <- length(intersect(egenes_bg, target))
  ft <- fisher.test(matrix(c(a, length(egenes_bg) - a, length(target) - a,
                             length(all_genes_bg) - length(union(egenes_bg, target))),
                           nrow = 2, byrow = TRUE))
  data.frame(n_genes = length(target), egene_count = a, OR = unname(ft$estimate),
             log2OR = log2(unname(ft$estimate)), p.value = ft$p.value,
             log2_ci_lower = log2(ft$conf.int[1]), log2_ci_upper = log2(ft$conf.int[2]))
}
group_levels <- c(first = "Low-dispersion", never = "Neither", fifth = "High-dispersion")
egene_by_seed <- bind_rows(lapply(names(mm$matched_by_seed), function(s) {
  bind_rows(lapply(names(group_levels), function(g) {
    cbind(seed = s, group = group_levels[[g]], egene_fisher(mm$matched_by_seed[[s]][[g]]))
  }))
}))
fig3a_df <- egene_by_seed %>% filter(seed == as.character(mm$match_seed_primary))
fig3a_df$group <- factor(fig3a_df$group, levels = unname(group_levels))
print(fig3a_df)
write.csv(fig3a_df, file.path(comb_dir, "egene_fisher_matched.csv"), row.names = FALSE)

# =============================================================================
# Fig4C-E: per-gene feature lookups (identical parsing to rnx 01)
# =============================================================================
tss_tab <- table(read.csv(fantom_tss_txt, header = FALSE)$V1)
tss_by_gene <- setNames(as.integer(tss_tab), names(tss_tab))

abc_pred <- read_table(abc_txt, show_col_types = FALSE)
enh_by_gene <- setNames(abc_pred$avg_length_across_celltype, abc_pred$TargetGene)

# ANALYSIS_OK[within-file-helper]: called below in this script only; local by design
strip_version <- function(x) sub("\\..*$", "", x)
mane_summary <- read.delim(mane_summary_txt, stringsAsFactors = FALSE, comment.char = "")
colnames(mane_summary)[1] <- sub("^#", "", colnames(mane_summary)[1])
mane_select <- mane_summary[mane_summary$MANE_status == "MANE Select", , drop = FALSE]
mane_gff <- read.delim(mane_gff_gz, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
exon_attr <- mane_gff$V9[mane_gff$V3 == "exon"]
exon_tx <- strip_version(sub("^transcript_id=", "", regmatches(exon_attr, regexpr("transcript_id=[^;]+", exon_attr))))
assert_that(length(exon_tx) == length(exon_attr), "failed to parse transcript_id from some GFF exon rows")
exon_counts <- table(exon_tx)
n_introns <- as.integer(exon_counts[strip_version(mane_select$Ensembl_nuc)]) - 1L
assert_that(!anyNA(n_introns), "some MANE Select transcripts have no exon rows")
intron_by_gene <- setNames(n_introns, mane_select$symbol)

# Low/High frame for one feature, in the rnx 01 layout (rank "1" = low, "5" = high).
# ANALYSIS_OK[within-file-helper]: called below for each feature and seed; local by design
low_high_df <- function(var, values_by_gene, matched) {
  df <- data.frame(genes = c(matched$first, matched$fifth),
                   rank = c(rep("1", length(matched$first)), rep("5", length(matched$fifth))))
  df[[var]] <- unname(values_by_gene[df$genes])
  # ANALYSIS_OK[filtering]: matched genes without a value in the external feature table are dropped (as in rnx 01 / 05); found counts logged below
  df <- df[!is.na(df[[var]]), , drop = FALSE]
  df$rank <- factor(df$rank, levels = c("1", "5"))
  df
}
features <- list(fig4c = list(var = "tss_count", values = tss_by_gene),
                 fig4d = list(var = "enh_length", values = enh_by_gene),
                 fig4e = list(var = "intron_count", values = intron_by_gene))

wilcox_rows <- list()
new_panels <- list()
for (nm in names(features)) {
  ft <- features[[nm]]
  for (s in names(mm$matched_by_seed)) {
    df <- low_high_df(ft$var, ft$values, mm$matched_by_seed[[s]])
    r <- wilcox_effect(df[[ft$var]], as.character(df$rank), "1", "5", ft$var)
    r$panel <- nm; r$seed <- s
    wilcox_rows[[paste(nm, s)]] <- r
    if (s == as.character(mm$match_seed_primary)) {
      log_step(sprintf("[%s] %d/%d low, %d/%d high matched genes have a value", nm,
                       sum(df$rank == "1"), length(primary$first),
                       sum(df$rank == "5"), length(primary$fifth)))
      new_panels[[nm]] <- df
    }
  }
}
wilcox_all <- bind_rows(wilcox_rows)
write.csv(wilcox_all %>% filter(seed == as.character(mm$match_seed_primary)),
          file.path(comb_dir, "wilcox_4c_4d_4e_matched.csv"), row.names = FALSE)

# ---- Seed stability (3A + 4C-E) -------------------------------------------
stab <- bind_rows(
  egene_by_seed %>% transmute(panel = "fig3a", term = group, seed, effect = log2OR, p_value = p.value),
  wilcox_all %>% transmute(panel, term = feature, seed, effect = hl_shift, p_value))
write.csv(stab, file.path(comb_dir, "seed_stability_3a_4c_4d_4e.csv"), row.names = FALSE)
sign_flips <- stab %>% group_by(panel, term) %>% summarise(n_sign = n_distinct(sign(effect)), .groups = "drop") %>%
  filter(n_sign > 1)
# 3A (eGene) is reference-only -- not in the composite since 2026-10-07 -- so a
# sign flip there is logged, not fatal. Its "Neither" log2OR sits near 0
# (|log2OR| < 0.7, mostly p > 0.05), so the sign is not stable after the
# 2026-10-08 -6 exclusion. Plotted panels 4C-E must still be sign-stable.
flips_3a <- sign_flips %>% filter(panel == "fig3a")
if (nrow(flips_3a) > 0) {
  log_step(paste("WARNING: reference-only fig3a effect sign changes across matching seeds for:",
                 paste(flips_3a$term, collapse = ", ")))
}
flips_plotted <- sign_flips %>% filter(panel != "fig3a")
assert_that(nrow(flips_plotted) == 0,
            paste("effect sign changes across matching seeds for:", paste(flips_plotted$panel, flips_plotted$term, collapse = ", ")))
log_step("4C / 4D / 4E effect directions are stable across all matching seeds")

# =============================================================================
# Gather 05's / copied-06's panel data (ggplot $data, no recomputation)
# =============================================================================
# ANALYSIS_OK[within-file-helper]: called below for each saved panel; local by design
panel_data_from_rds <- function(filename, value_col) {
  f <- file.path(panel_dir, filename)
  assert_that(file.exists(f), sprintf("missing %s -- run 05 / mm_copy_fig3c first", f))
  d <- readRDS(f)$data
  assert_that(all(c("gene", "rank_category", value_col) %in% colnames(d)), paste("unexpected columns in", f))
  assert_that(setequal(unique(as.character(d$rank_category)), c("Low-dispersion", "High-dispersion")),
              paste("unexpected groups in", f))
  data.frame(genes = d$gene,
             rank = factor(ifelse(d$rank_category == "Low-dispersion", "1", "5"), levels = c("1", "5")),
             setNames(list(d[[value_col]]), value_col))
}
tata <- readRDS(file.path(panel_dir, "fig4b_tata_plot.rds"))$data
assert_that(nrow(tata) == 2, "fig4b TATA panel data should have 2 rows")

panel_data <- list(
  fig3a      = fig3a_df,
  fig3b      = panel_data_from_rds("fig3b_twn_hub_rank_plot.rds", "hub_rank"),
  fig3c_fib  = panel_data_from_rds("fig3c_Cardiac_fibroblasts_plot.rds", "n_genes"),
  fig3c_card = panel_data_from_rds("fig3c_Cardiomyocytes_plot.rds", "n_genes"),
  fig3d      = panel_data_from_rds("fig3d_ppi_connectivity_plot.rds", "int_count"),
  fig4a      = panel_data_from_rds("fig4a_de_count_plot.rds", "DE_count"),
  fig4b      = tata,
  fig4c      = new_panels$fig4c,
  fig4d      = new_panels$fig4d,
  fig4e      = new_panels$fig4e)
saveRDS(panel_data, file.path(comb_dir, "panel_data.rds"))

q <- function(k) wilcox_all$p_value[wilcox_all$panel == k & wilcox_all$seed == as.character(mm$match_seed_primary)]
vals <- list(
  mm_fig3a_low_log2OR = fig3a_df$log2OR[fig3a_df$group == "Low-dispersion"],
  mm_fig3a_low_p = fig3a_df$p.value[fig3a_df$group == "Low-dispersion"],
  mm_fig3a_high_log2OR = fig3a_df$log2OR[fig3a_df$group == "High-dispersion"],
  mm_fig3a_high_p = fig3a_df$p.value[fig3a_df$group == "High-dispersion"],
  mm_fig4c_tss_wilcox_p = q("fig4c"), mm_fig4d_enh_wilcox_p = q("fig4d"), mm_fig4e_intron_wilcox_p = q("fig4e"))
# Reportable values (R-side register_value shim, per the figure_reproduction precedent).
write.csv(data.frame(key = names(vals), value = signif(unlist(vals), 4),
                     provenance = "figure_reproduction/scripts/mean_matched_regulatory_revision/08_remaining_panels_matched.R"),
          file.path(comb_dir, "reportable_values_mm_combined.csv"), row.names = FALSE)
log_step(paste("wrote combined panel data to", comb_dir))
