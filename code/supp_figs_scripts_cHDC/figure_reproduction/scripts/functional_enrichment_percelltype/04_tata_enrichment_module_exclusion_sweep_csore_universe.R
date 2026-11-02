#!/usr/bin/env Rscript

stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(ggplot2)
})

OUT_DIR <- "figure_reproduction/outputs/functional_enrichment_percelltype"
QUINTILE_CSV <- file.path(OUT_DIR, "percelltype_quintile_assignments_csore_universe.csv")
BASELINE_CSV <- file.path(OUT_DIR, "tata_enrichment_percelltype_csore_universe_results.csv")
EPD_BED <- "cellranger_cluster-mode_trial/analysis/datasets/external_data/epd/human_epdnew_3sa0y.bed"
CONSENSUS_MODULES_DIR <- "figure_reproduction/outputs/coexpression_plasticity_lowfidelity/consensus_modules_prototype"
MODULE_SUMMARY_CSV <- "figure_reproduction/outputs/module_dispersion_profile/phase2_consensus/module_dispersion_summary_consensus.csv"
CONSENSUS_QUANTILE_SUFFIX <- "q0.5"

for (f in c(QUINTILE_CSV, BASELINE_CSV, EPD_BED, MODULE_SUMMARY_CSV)) {
  if (!file.exists(f)) stop("Missing required input: ", f, " -- run script 03 and module_dispersion_profile scripts 02 first.")
}

CELL_TYPES <- c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                "Epithelial_cells", "Activated_fibroblasts",
                "Foregut_endodermal_cells")
TOP_QUINTILE <- 5L
EXCLUSION_FRACTIONS <- seq(0.1, 0.7, by = 0.1)

celltype_colors <- c(
  "Activated fibroblasts"    = "#117733", "Cardiac fibroblasts" = "#4477AA",
  "Cardiomyocytes"           = "#882255", "Epithelial cells"    = "#AA4499",
  "Foregut endodermal cells" = "#CC6677", "Mesenchymal cells"   = "#DDCC77"
)

log_msg <- function(...) message("[04_tata_enrichment_module_exclusion_sweep_csore_universe] ", ...)

quintile_df <- read.csv(QUINTILE_CSV, stringsAsFactors = FALSE)
stopifnot(all(c("gene", "cell_type", "quintile") %in% colnames(quintile_df)))

baseline_df <- read.csv(BASELINE_CSV, stringsAsFactors = FALSE)
stopifnot(all(c("cell_type", "gene_set", "odds_ratio", "ci_lower", "ci_upper", "p_value", "n_target") %in% colnames(baseline_df)))

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

result_rows <- list()

for (ct in CELL_TYPES) {
  ct_quintiles <- quintile_df[quintile_df$cell_type == ct, ]
  background_genes <- ct_quintiles$gene
  n_background <- length(background_genes)
  high_disp_genes <- ct_quintiles$gene[ct_quintiles$quintile == TOP_QUINTILE]
  n_high_disp <- length(high_disp_genes)
  tata_genes <- intersect(unique(epd_tata$gene), background_genes)

  module_csv <- file.path(CONSENSUS_MODULES_DIR, paste0(ct, "_consensus_", CONSENSUS_QUANTILE_SUFFIX, "_module_membership.csv"))
  if (!file.exists(module_csv)) stop("Missing consensus module membership file: ", module_csv)
  module_df <- read.csv(module_csv, stringsAsFactors = FALSE)
  stopifnot(all(c("gene", "module", "in_module") %in% colnames(module_df)))
  real_module_df <- module_df[module_df$module != 0 & module_df$in_module, , drop = FALSE]

  summary_df <- read.csv(MODULE_SUMMARY_CSV, stringsAsFactors = FALSE)
  ct_summary <- summary_df[summary_df$cell_type == ct, ]
  stopifnot("module_rank" %in% colnames(ct_summary), nrow(ct_summary) > 0)
  n_modules <- nrow(ct_summary)

  n_high_disp_in_consensus_universe <- sum(high_disp_genes %in% module_df$gene)
  n_high_disp_in_real_module <- sum(high_disp_genes %in% real_module_df$gene)
  log_msg("=== ", ct, " === background=", n_background, "; ", n_high_disp,
          " high-dispersion genes; ", n_high_disp_in_consensus_universe,
          " in the ", n_modules, "-module consensus universe (should be ~all, by construction); ",
          n_high_disp_in_real_module, " in a real (non-module-0) module.")

  base_row <- baseline_df[baseline_df$cell_type == ct & baseline_df$gene_set == "high_disp", ]
  stopifnot(nrow(base_row) == 1)
  result_rows[[paste0(ct, "_x0")]] <- data.frame(
    cell_type = ct, exclusion_fraction = 0, n_modules_total = n_modules, n_modules_excluded = 0L,
    n_high_disp_original = n_high_disp, n_high_disp_excluded = 0L, n_high_disp_remaining = n_high_disp,
    n_target = base_row$n_target, odds_ratio = base_row$odds_ratio,
    ci_lower = base_row$ci_lower, ci_upper = base_row$ci_upper, p_value = base_row$p_value,
    stringsAsFactors = FALSE
  )

  for (x in EXCLUSION_FRACTIONS) {
    n_exclude <- ceiling(x * n_modules)
    top_modules <- ct_summary$module[order(ct_summary$module_rank, decreasing = TRUE)][seq_len(n_exclude)]
    excluded_genes <- real_module_df$gene[real_module_df$module %in% top_modules]

    reduced_target <- setdiff(high_disp_genes, excluded_genes)
    n_removed <- n_high_disp - length(reduced_target)
    res <- tata_enrich_func(reduced_target, background_genes, tata_genes)
    log_msg("  x=", x, ": excluded ", n_exclude, "/", n_modules, " modules -> removed ",
            n_removed, "/", n_high_disp, " high-dispersion genes; n_remaining=", length(reduced_target),
            " OR=", round(as.numeric(res$ft$estimate), 3), " p=", signif(res$ft$p.value, 3))

    result_rows[[paste0(ct, "_x", x)]] <- data.frame(
      cell_type = ct, exclusion_fraction = x, n_modules_total = n_modules, n_modules_excluded = n_exclude,
      n_high_disp_original = n_high_disp, n_high_disp_excluded = n_removed,
      n_high_disp_remaining = length(reduced_target),
      n_target = length(reduced_target), odds_ratio = as.numeric(res$ft$estimate),
      ci_lower = res$ft$conf.int[1], ci_upper = res$ft$conf.int[2], p_value = res$ft$p.value,
      stringsAsFactors = FALSE
    )
  }
}

results_df <- do.call(rbind, result_rows); rownames(results_df) <- NULL
stopifnot(nrow(results_df) == length(CELL_TYPES) * (length(EXCLUSION_FRACTIONS) + 1))
results_csv <- file.path(OUT_DIR, "tata_enrichment_module_exclusion_sweep_csore_universe_results.csv")
write.csv(results_df, results_csv, row.names = FALSE)
log_msg("Wrote: ", results_csv)
print(results_df[, c("cell_type", "exclusion_fraction", "n_high_disp_remaining", "odds_ratio", "p_value")])

plot_df <- results_df
plot_df$cell_type_label <- factor(gsub("_", " ", plot_df$cell_type), levels = gsub("_", " ", CELL_TYPES))
plot_df$log2_or <- log2(plot_df$odds_ratio)
plot_df$log2_ci_lower <- log2(plot_df$ci_lower)
plot_df$log2_ci_upper <- log2(plot_df$ci_upper)

p <- ggplot(plot_df, aes(x = exclusion_fraction, y = log2_or, color = cell_type_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey40") +
  geom_errorbar(aes(ymin = log2_ci_lower, ymax = log2_ci_upper), width = 0.01, linewidth = 0.3, alpha = 0.5) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.4) +
  scale_color_manual(values = celltype_colors, name = NULL) +
  scale_x_continuous(breaks = c(0, EXCLUSION_FRACTIONS), labels = scales::percent) +
  labs(x = "Fraction of highest-dispersion consensus modules excluded",
       y = expression(log[2]("Odds ratio")),
       title = "TATA enrichment vs. top-x% consensus module exclusion, CS-CORE consensus universe",
       subtitle = "x=0 is script 03's unmodified baseline; background/target both restricted to the consensus-module universe throughout") +
  theme_classic(base_family = "Arial") +
  theme(
    text = element_text(size = 6), axis.title = element_text(size = 6), axis.text = element_text(size = 6),
    axis.line = element_line(linewidth = 0.4), legend.text = element_text(size = 6),
    legend.title = element_text(size = 6), legend.key.size = unit(6, "pt"),
    plot.title = element_text(size = 7, face = "bold"), plot.subtitle = element_text(size = 5.5),
    plot.margin = margin(2, 2, 2, 2, unit = "mm"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4)
  )

out_pdf <- file.path(OUT_DIR, "tata_enrichment_module_exclusion_sweep_csore_universe.pdf")
ggsave(out_pdf, p, width = 5, height = 3.5, units = "in", device = cairo_pdf)
log_msg("Wrote: ", out_pdf)

log_msg("DONE.")
