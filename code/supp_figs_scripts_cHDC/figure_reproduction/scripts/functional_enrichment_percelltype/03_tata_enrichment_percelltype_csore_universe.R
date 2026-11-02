#!/usr/bin/env Rscript

stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(ggplot2)
})

COEXP_OUT_DIR <- "figure_reproduction/outputs/coexpression_plasticity_lowfidelity"
RES_VAR_DIR   <- file.path(COEXP_OUT_DIR, "memento_per_celltype")
CONSENSUS_DIR <- file.path(COEXP_OUT_DIR, "consensus_modules_prototype")
CONSENSUS_QUANTILE_SUFFIX <- "q0.5"
EPD_BED <- "cellranger_cluster-mode_trial/analysis/datasets/external_data/epd/human_epdnew_3sa0y.bed"
OUT_DIR <- "figure_reproduction/outputs/functional_enrichment_percelltype"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CELL_TYPES <- c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                "Epithelial_cells", "Activated_fibroblasts",
                "Foregut_endodermal_cells")
INDIVIDUALS <- c("NA19099", "NA19093", "NA18522")
N_QUANTILES <- 5L
BOTTOM_QUINTILE <- 1L
TOP_QUINTILE <- N_QUANTILES
MAX_DROP_FRAC <- 0.05

celltype_colors <- c(
  "Activated fibroblasts"    = "#117733", "Cardiac fibroblasts" = "#4477AA",
  "Cardiomyocytes"           = "#882255", "Epithelial cells"    = "#AA4499",
  "Foregut endodermal cells" = "#CC6677", "Mesenchymal cells"   = "#DDCC77"
)

log_msg <- function(...) message("[03_tata_enrichment_percelltype_csore_universe] ", ...)

values_to_quintiles <- function(values) {
  ranks <- rank(values, ties.method = "average")
  cut(ranks,
      breaks = quantile(ranks, probs = seq(0, 1, length.out = N_QUANTILES + 1), na.rm = TRUE),
      labels = FALSE, include.lowest = TRUE)
}

if (!file.exists(EPD_BED)) stop("Missing ", EPD_BED)
epd_tata <- read.table(EPD_BED, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(epd_tata) <- c("chr", "start", "end", "promoter_id", "score", "strand")
epd_tata$gene <- sub("_[0-9]+$", "", epd_tata$promoter_id)

tata_enrich_func <- function(target_genes, background_genes, tata_genes) {
  a <- length(intersect(tata_genes, target_genes))
  b <- length(setdiff(tata_genes, target_genes))
  c_ <- length(setdiff(target_genes, tata_genes))
  d <- length(setdiff(background_genes, union(tata_genes, target_genes)))
  stopifnot(a >= 0, b >= 0, c_ >= 0, d >= 0, a + b + c_ + d == length(background_genes))
  contingency <- matrix(c(a, b, c_, d), nrow = 2, byrow = TRUE,
                         dimnames = list(TATA = c("Yes", "No"), Target = c("In_Set", "Not_In_Set")))
  ft <- fisher.test(contingency)
  list(ft = ft, a = a)
}

quintile_rows <- list()
result_rows <- list()

for (ct in CELL_TYPES) {
  module_csv <- file.path(CONSENSUS_DIR, paste0(ct, "_consensus_", CONSENSUS_QUANTILE_SUFFIX, "_module_membership.csv"))
  if (!file.exists(module_csv)) stop("Missing consensus module membership file: ", module_csv)
  module_df <- read.csv(module_csv, stringsAsFactors = FALSE)
  stopifnot(all(c("gene", "module", "in_module") %in% colnames(module_df)))
  consensus_universe <- module_df$gene

  res_var_csv <- file.path(RES_VAR_DIR, paste0(ct, "_res_var.csv"))
  if (!file.exists(res_var_csv)) stop("Missing res_var file: ", res_var_csv)
  res_var_df <- read.csv(res_var_csv, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot("res_var.csv missing 'gene' column" = "gene" %in% colnames(res_var_df))
  stopifnot("res_var.csv missing an expected individual column" = all(INDIVIDUALS %in% colnames(res_var_df)))
  if (anyDuplicated(res_var_df$gene)) stop("Duplicate gene rows in ", res_var_csv)

  rv_matched <- res_var_df[match(consensus_universe, res_var_df$gene), c("gene", INDIVIDUALS)]
  has_all_3 <- !is.na(rv_matched$gene) & rowSums(is.na(rv_matched[, INDIVIDUALS, drop = FALSE])) == 0
  drop_frac <- sum(!has_all_3) / length(consensus_universe)
  if (drop_frac > MAX_DROP_FRAC) {
    stop(sprintf("[%s] %.1f%% of the consensus universe lacks a res_var value in all 3 individuals (>%.0f%% cap)",
                  ct, 100 * drop_frac, 100 * MAX_DROP_FRAC))
  }
  background_genes <- consensus_universe[has_all_3]
  n_background <- length(background_genes)
  rv_sub <- rv_matched[has_all_3, , drop = FALSE]
  stopifnot(identical(rv_sub$gene, background_genes))
  median_disp_per_gene <- apply(as.matrix(rv_sub[, INDIVIDUALS, drop = FALSE]), 1, median)

  quintiles <- values_to_quintiles(median_disp_per_gene)
  names(quintiles) <- background_genes
  quintile_rows[[ct]] <- data.frame(gene = background_genes, cell_type = ct, quintile = quintiles, stringsAsFactors = FALSE)

  tata_genes <- intersect(unique(epd_tata$gene), background_genes)
  log_msg("=== ", ct, " === background=", n_background, " (CS-CORE consensus universe); TATA genes in universe=", length(tata_genes))

  for (quint_label in c("low_disp", "high_disp")) {
    q <- if (quint_label == "low_disp") BOTTOM_QUINTILE else TOP_QUINTILE
    target_genes <- background_genes[quintiles == q]
    res <- tata_enrich_func(target_genes, background_genes, tata_genes)
    log_msg("  ", quint_label, ": n_target=", length(target_genes),
            " OR=", round(as.numeric(res$ft$estimate), 3), " p=", signif(res$ft$p.value, 3))
    result_rows[[paste0(ct, "_", quint_label)]] <- data.frame(
      cell_type = ct, gene_set = quint_label, n_target = length(target_genes),
      n_tata_in_target = res$a, odds_ratio = as.numeric(res$ft$estimate),
      ci_lower = res$ft$conf.int[1], ci_upper = res$ft$conf.int[2],
      p_value = res$ft$p.value, stringsAsFactors = FALSE
    )
  }
}

quintile_df <- do.call(rbind, quintile_rows); rownames(quintile_df) <- NULL
quintile_csv <- file.path(OUT_DIR, "percelltype_quintile_assignments_csore_universe.csv")
write.csv(quintile_df, quintile_csv, row.names = FALSE)
log_msg("Wrote: ", quintile_csv)

results_df <- do.call(rbind, result_rows); rownames(results_df) <- NULL
stopifnot(nrow(results_df) == length(CELL_TYPES) * 2)
results_csv <- file.path(OUT_DIR, "tata_enrichment_percelltype_csore_universe_results.csv")
write.csv(results_df, results_csv, row.names = FALSE)
log_msg("Wrote: ", results_csv)
print(results_df[, c("cell_type", "gene_set", "n_target", "odds_ratio", "p_value")])

plot_df <- results_df
plot_df$cell_type_label <- factor(gsub("_", " ", plot_df$cell_type), levels = gsub("_", " ", CELL_TYPES))
plot_df$gene_set <- factor(plot_df$gene_set, levels = c("low_disp", "high_disp"))
plot_df$log2_or <- log2(plot_df$odds_ratio)
plot_df$log2_ci_lower <- log2(plot_df$ci_lower)
plot_df$log2_ci_upper <- log2(plot_df$ci_upper)

p <- ggplot(plot_df, aes(x = gene_set, y = log2_or, color = cell_type_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey40") +
  geom_errorbar(aes(ymin = log2_ci_lower, ymax = log2_ci_upper),
                width = 0.15, linewidth = 0.4, position = position_dodge(width = 0.6)) +
  geom_point(size = 1.6, position = position_dodge(width = 0.6)) +
  scale_color_manual(values = celltype_colors, name = NULL) +
  scale_x_discrete(labels = c(low_disp = "Bottom 20%\n(low dispersion)", high_disp = "Top 20%\n(high dispersion)")) +
  labs(x = NULL, y = expression(log[2]("Odds ratio")),
       title = "TATA-promoter enrichment, per-cell-type top/bottom 20%, CS-CORE consensus universe",
       subtitle = "Same method as script 01, but background/target restricted to the ~5,000-gene-per-individual CS-CORE consensus universe") +
  theme_classic(base_family = "Arial") +
  theme(
    text = element_text(size = 6), axis.title = element_text(size = 6), axis.text = element_text(size = 6),
    axis.line = element_line(linewidth = 0.4), legend.text = element_text(size = 6),
    legend.title = element_text(size = 6), legend.key.size = unit(6, "pt"),
    plot.title = element_text(size = 7, face = "bold"), plot.subtitle = element_text(size = 5.5),
    plot.margin = margin(2, 2, 2, 2, unit = "mm"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4)
  )

out_pdf <- file.path(OUT_DIR, "tata_enrichment_percelltype_csore_universe.pdf")
ggsave(out_pdf, p, width = 4.5, height = 3.2, units = "in", device = cairo_pdf)
log_msg("Wrote: ", out_pdf)

log_msg("DONE.")
