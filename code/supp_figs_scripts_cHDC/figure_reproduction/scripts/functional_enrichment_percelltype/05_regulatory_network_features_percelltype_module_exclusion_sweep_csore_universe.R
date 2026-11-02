#!/usr/bin/env Rscript

stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
})

OUT_DIR <- "figure_reproduction/outputs/functional_enrichment_percelltype"
QUINTILE_CSV <- file.path(OUT_DIR, "percelltype_quintile_assignments_csore_universe.csv")
CONSENSUS_MODULES_DIR <- "figure_reproduction/outputs/coexpression_plasticity_lowfidelity/consensus_modules_prototype"
MODULE_SUMMARY_CSV <- "figure_reproduction/outputs/module_dispersion_profile/phase2_consensus/module_dispersion_summary_consensus.csv"
CONSENSUS_QUANTILE_SUFFIX <- "q0.5"

# ANALYSIS_OK[pinned-input]: external, read-only data sources identical to
# the original analysis -- fixed upstream files, untouched by the Epithelial
# sub-cluster-24 exclusion, not something to duplicate.
TSS_TXT <- "/project2/gilad/wenhe/cox_baboon_diet/combined_analysis/gene_annotation/fantom/short_descrip_split_gene_only.txt"
ABC_ENHANCER_TXT <- "/project2/gilad/wenhe/cox_baboon_diet/combined_analysis/gene_annotation/ABC_prediction_score005.txt"
PPI_TSV <- "cellranger_cluster-mode_trial/analysis/datasets/external_data/protein_interaction_data/interaction_consensus.tsv.zip"
PPI_TISSUE_TSV <- "cellranger_cluster-mode_trial/analysis/datasets/external_data/protein_interaction_data/has_data_Tissue.tsv"
TWN_S4_XLSX <- "cellranger_cluster-mode_trial/analysis/datasets/external_data/saha_coexp_net_TWN_2017/Supplemental_Data_S4.xlsx"
MANE_SUMMARY_TXT <- "data/raw/mane_v1_5/MANE.GRCh38.v1.5.summary.txt.gz"
MANE_GFF_GZ      <- "data/raw/mane_v1_5/MANE.GRCh38.v1.5.ensembl_genomic.gff.gz"

for (f in c(QUINTILE_CSV, MODULE_SUMMARY_CSV, TSS_TXT, ABC_ENHANCER_TXT, PPI_TSV, PPI_TISSUE_TSV,
            TWN_S4_XLSX, MANE_SUMMARY_TXT, MANE_GFF_GZ)) {
  if (!file.exists(f)) stop("Missing required input: ", f, " -- run script 03 first for the quintile CSV.")
}

CELL_TYPES <- c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                "Epithelial_cells", "Activated_fibroblasts",
                "Foregut_endodermal_cells")
LOW_QUINTILE  <- 1L
HIGH_QUINTILE <- 5L
EXCLUSION_FRACTIONS <- seq(0.1, 0.7, by = 0.1)

log_msg <- function(...) message("[05_regulatory_network_features_percelltype_module_exclusion_sweep_csore_universe] ", ...)

wilcox_effect_two <- function(x, y, feature) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  stopifnot("empty Low or High group in wilcox_effect_two" = length(x) > 0 && length(y) > 0)
  wt <- wilcox.test(x, y, conf.int = TRUE)
  data.frame(
    feature = feature, n_low = length(x), n_high = length(y),
    median_low = median(x), median_high = median(y),
    hl_shift = unname(wt$estimate), ci_lower = wt$conf.int[1], ci_upper = wt$conf.int[2],
    p_value = wt$p.value, stringsAsFactors = FALSE
  )
}

tss_usage <- read.csv(TSS_TXT, header = FALSE, stringsAsFactors = FALSE)
tss_table <- as.data.frame(table(tss_usage$V1), stringsAsFactors = FALSE)
colnames(tss_table) <- c("gene", "tss_count")
tss_lookup <- setNames(tss_table$tss_count, tss_table$gene)
log_msg("TSS count table: ", nrow(tss_table), " genes (from ", nrow(tss_usage), " CAGE TSS rows)")

abc_pred <- readr::read_table(ABC_ENHANCER_TXT, show_col_types = FALSE)
stopifnot("ABC enhancer file missing expected columns" =
            all(c("TargetGene", "avg_length_across_celltype") %in% colnames(abc_pred)))
stopifnot("ABC enhancer file has duplicate TargetGene rows -- unexpected" = !anyDuplicated(abc_pred$TargetGene))
enh_lookup <- setNames(abc_pred$avg_length_across_celltype, abc_pred$TargetGene)
log_msg("Enhancer length table: ", length(enh_lookup), " genes")

ppi <- readr::read_tsv(PPI_TSV, show_col_types = FALSE)
exp_evidence <- readr::read_tsv(PPI_TISSUE_TSV, show_col_types = FALSE)
ppi$gene_1 <- exp_evidence$Gene[match(ppi$ensembl_gene_id_1, exp_evidence$Ensembl)]
ppi_table <- as.data.frame(table(ppi$gene_1), stringsAsFactors = FALSE)
colnames(ppi_table) <- c("gene", "n_interactions")
ppi_lookup <- setNames(ppi_table$n_interactions, ppi_table$gene)
log_msg("PPI connectivity table: ", nrow(ppi_table), " genes with >=1 interaction")

s4_sheets <- readxl::excel_sheets(TWN_S4_XLSX)
stopifnot("Supplemental_Data_S4.xlsx sheet 10 layout changed -- expected 'Group3 (TE-TE)'" =
            grepl("^Group3", s4_sheets[10]))
shared_te_te <- readxl::read_excel(TWN_S4_XLSX, sheet = s4_sheets[10], col_names = TRUE)
shared_te_te$rank <- as.integer(shared_te_te$rank)
twn_lookup <- setNames(shared_te_te$rank, shared_te_te$hub)
log_msg("TWN shared-hub rank table: ", length(twn_lookup), " hub genes")

strip_version <- function(x) sub("\\..*$", "", x)
mane_summary <- read.delim(MANE_SUMMARY_TXT, stringsAsFactors = FALSE, comment.char = "")
colnames(mane_summary)[1] <- sub("^#", "", colnames(mane_summary)[1])
mane_select <- mane_summary[mane_summary$MANE_status == "MANE Select", , drop = FALSE]
mane_select$ensembl_transcript_id <- strip_version(mane_select$Ensembl_nuc)
log_msg("MANE Select transcripts (genome-wide): ", nrow(mane_select))

mane_gff <- read.delim(MANE_GFF_GZ, header = FALSE, comment.char = "#", stringsAsFactors = FALSE,
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
intron_lookup <- setNames(mane_select$num_introns, mane_select$symbol)
log_msg("Intron count table: ", length(intron_lookup), " genes (MANE Select, genome-wide)")

FEATURE_LOOKUPS <- list(
  tss_count        = tss_lookup,
  enhancer_length  = enh_lookup,
  ppi_connectivity = ppi_lookup,
  twn_hub_rank     = twn_lookup,
  intron_count     = intron_lookup
)

quintile_df <- read.csv(QUINTILE_CSV, stringsAsFactors = FALSE)
stopifnot("percelltype_quintile_assignments_csore_universe.csv missing expected columns" =
            all(c("gene", "cell_type", "quintile") %in% colnames(quintile_df)))

summary_df <- read.csv(MODULE_SUMMARY_CSV, stringsAsFactors = FALSE)
stopifnot("module_dispersion_summary_consensus.csv missing expected columns" =
            all(c("cell_type", "module", "module_rank") %in% colnames(summary_df)))

result_rows <- list()

for (ct in CELL_TYPES) {
  ct_quintiles <- quintile_df[quintile_df$cell_type == ct, ]
  low_genes <- ct_quintiles$gene[ct_quintiles$quintile == LOW_QUINTILE]
  high_genes_full <- ct_quintiles$gene[ct_quintiles$quintile == HIGH_QUINTILE]

  module_csv <- file.path(CONSENSUS_MODULES_DIR, paste0(ct, "_consensus_", CONSENSUS_QUANTILE_SUFFIX, "_module_membership.csv"))
  if (!file.exists(module_csv)) stop("Missing consensus module membership file: ", module_csv)
  module_df <- read.csv(module_csv, stringsAsFactors = FALSE)
  stopifnot("consensus module_membership.csv missing expected columns" =
              all(c("gene", "module", "in_module") %in% colnames(module_df)))
  real_module_df <- module_df[module_df$module != 0 & module_df$in_module, , drop = FALSE]

  ct_summary <- summary_df[summary_df$cell_type == ct, ]
  stopifnot("module_rank" %in% colnames(ct_summary), nrow(ct_summary) > 0)
  n_modules <- nrow(ct_summary)

  log_msg("=== ", ct, " === low=", length(low_genes), " high=", length(high_genes_full),
          " (", n_modules, " consensus modules)")

  for (feat_name in names(FEATURE_LOOKUPS)) {
    lookup <- FEATURE_LOOKUPS[[feat_name]]
    values_low <- lookup[low_genes]

    for (x in c(0, EXCLUSION_FRACTIONS)) {
      if (x == 0) {
        n_exclude <- 0L
        excluded_genes <- character(0)
      } else {
        n_exclude <- ceiling(x * n_modules)
        top_modules <- ct_summary$module[order(ct_summary$module_rank, decreasing = TRUE)][seq_len(n_exclude)]
        excluded_genes <- real_module_df$gene[real_module_df$module %in% top_modules]
      }
      reduced_high <- setdiff(high_genes_full, excluded_genes)
      values_high <- lookup[reduced_high]

      res <- wilcox_effect_two(values_low, values_high, feat_name)
      log_msg("  ", feat_name, " x=", x, ": n_low_with_value=", res$n_low,
              " n_high_with_value=", res$n_high, " (", length(reduced_high), "/", length(high_genes_full),
              " high-dispersion genes remain) HL_shift=", round(res$hl_shift, 3),
              " p=", signif(res$p_value, 3))

      result_rows[[paste(ct, feat_name, x, sep = "_")]] <- data.frame(
        cell_type = ct, feature = feat_name, exclusion_fraction = x,
        n_modules_total = n_modules, n_modules_excluded = n_exclude,
        n_high_dispersion_original = length(high_genes_full),
        n_high_dispersion_remaining = length(reduced_high),
        res, stringsAsFactors = FALSE
      )
    }
  }
}

results_df <- do.call(rbind, result_rows); rownames(results_df) <- NULL
stopifnot("expected 1 row per cell type x feature x exclusion fraction" =
            nrow(results_df) == length(CELL_TYPES) * length(FEATURE_LOOKUPS) * (length(EXCLUSION_FRACTIONS) + 1))
results_csv <- file.path(OUT_DIR, "regulatory_network_features_percelltype_module_exclusion_sweep_csore_universe_results.csv")
write.csv(results_df, results_csv, row.names = FALSE)
log_msg("Wrote: ", results_csv, " (", nrow(results_df), " rows)")

log_msg("=== x=0 baseline (Low vs High, no exclusion) ===")
print(results_df[results_df$exclusion_fraction == 0,
                  c("cell_type", "feature", "n_low", "n_high", "median_low", "median_high", "hl_shift", "p_value")])

log_msg("DONE. Run script 05b to (re)generate the combined plot.")
