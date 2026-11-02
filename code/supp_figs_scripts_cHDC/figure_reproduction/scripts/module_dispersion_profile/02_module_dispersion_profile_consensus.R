#!/usr/bin/env Rscript

stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

# ---------------------------------------------------------------------------
# NAMED CONSTANTS
# ---------------------------------------------------------------------------
COEXP_OUT_DIR <- "figure_reproduction/outputs/coexpression_plasticity_lowfidelity"
RES_VAR_DIR   <- file.path(COEXP_OUT_DIR, "memento_per_celltype")
CONSENSUS_DIR <- file.path(COEXP_OUT_DIR, "consensus_modules_prototype")
OUT_DIR       <- "figure_reproduction/outputs/module_dispersion_profile/phase2_consensus"
if (!dir.exists(CONSENSUS_DIR)) {
  stop("Missing consensus module directory: ", CONSENSUS_DIR,
       " -- run the coexpression_plasticity_lowfidelity rules first.")
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CELL_TYPES <- c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                "Epithelial_cells", "Activated_fibroblasts",
                "Foregut_endodermal_cells")
INDIVIDUALS <- c("NA19099", "NA19093", "NA18522")
CONSENSUS_QUANTILE_SUFFIX <- "q0.5"
N_QUANTILES <- 5L
TARGET_QUINTILE <- N_QUANTILES
MAX_DROP_FRAC <- 0.05
MIN_EXPECTED_MODULE_GENES <- 1000L

log_msg <- function(...) message("[02_module_dispersion_profile_consensus] ", ...)

values_to_quintiles <- function(values) {
  ranks <- rank(values, ties.method = "average")
  cut(ranks,
      breaks = quantile(ranks, probs = seq(0, 1, length.out = N_QUANTILES + 1), na.rm = TRUE),
      labels = FALSE, include.lowest = TRUE)
}

celltype_colors <- c(
  "Activated fibroblasts"    = "#117733",
  "Cardiac fibroblasts"      = "#4477AA",
  "Cardiomyocytes"           = "#882255",
  "Epithelial cells"         = "#AA4499",
  "Foregut endodermal cells" = "#CC6677",
  "Mesenchymal cells"        = "#DDCC77"
)

detail_rows <- list()
summary_rows <- list()

for (ct in CELL_TYPES) {
  log_msg("=== ", ct, " ===")

  module_csv <- file.path(CONSENSUS_DIR, paste0(ct, "_consensus_", CONSENSUS_QUANTILE_SUFFIX, "_module_membership.csv"))
  if (!file.exists(module_csv)) stop("Missing consensus module membership file: ", module_csv)
  module_df <- read.csv(module_csv, stringsAsFactors = FALSE)
  stopifnot(
    "consensus module_membership.csv missing expected columns" =
      all(c("gene", "module", "in_module") %in% colnames(module_df)),
    "consensus module_membership.csv gene count below the expected consensus-universe floor" =
      nrow(module_df) >= MIN_EXPECTED_MODULE_GENES,
    "consensus module_membership.csv has duplicate genes" = !anyDuplicated(module_df$gene)
  )

  res_var_csv <- file.path(RES_VAR_DIR, paste0(ct, "_res_var.csv"))
  if (!file.exists(res_var_csv)) stop("Missing res_var file: ", res_var_csv)
  res_var_df <- read.csv(res_var_csv, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot("res_var.csv missing 'gene' column" = "gene" %in% colnames(res_var_df))
  stopifnot("res_var.csv missing an expected individual column" =
              all(INDIVIDUALS %in% colnames(res_var_df)))
  if (anyDuplicated(res_var_df$gene)) stop("Duplicate gene rows in ", res_var_csv)

  consensus_universe <- module_df$gene
  rv_matched <- res_var_df[match(consensus_universe, res_var_df$gene), c("gene", INDIVIDUALS)]
  has_all_3 <- !is.na(rv_matched$gene) & rowSums(is.na(rv_matched[, INDIVIDUALS, drop = FALSE])) == 0
  n_dropped <- sum(!has_all_3)
  drop_frac <- n_dropped / length(consensus_universe)
  if (drop_frac > MAX_DROP_FRAC) {
    stop(sprintf("[%s] %.1f%% of the %d-gene consensus universe lacks a res_var value in all 3 individuals (>%.0f%% cap) -- investigate before proceeding",
                  ct, 100 * drop_frac, length(consensus_universe), 100 * MAX_DROP_FRAC))
  }
  if (n_dropped > 0) {
    log_msg("  dropping ", n_dropped, " / ", length(consensus_universe),
            " consensus-universe genes missing a res_var value in >=1 individual (",
            round(100 * drop_frac, 2), "%)")
  }
  universe <- consensus_universe[has_all_3]
  n_universe <- length(universe)

  module_sub <- module_df[module_df$gene %in% universe, , drop = FALSE]
  rv_sub <- rv_matched[has_all_3, , drop = FALSE]
  stopifnot(identical(rv_sub$gene, universe))
  median_disp_per_gene <- apply(as.matrix(rv_sub[, INDIVIDUALS, drop = FALSE]), 1, median)
  dispersion_lookup <- setNames(median_disp_per_gene, universe)

  quintiles <- values_to_quintiles(median_disp_per_gene)
  high_disp_genes <- universe[quintiles == TARGET_QUINTILE]
  n_high_disp <- length(high_disp_genes)
  log_msg("  universe=", n_universe, " high_dispersion(top quintile, median-across-individuals)=",
          n_high_disp, " (", round(100 * n_high_disp / n_universe, 1), "%)")

  real_module_sub <- module_sub[module_sub$module != 0 & module_sub$in_module, , drop = FALSE]
  real_modules <- sort(unique(real_module_sub$module))
  n_modules <- length(real_modules)
  stopifnot("at least one real consensus module expected per cell type" = n_modules > 0)

  n_high_disp_in_module0 <- n_high_disp - sum(high_disp_genes %in% real_module_sub$gene)
  log_msg("  ", n_modules, " real consensus module(s); ", n_high_disp - n_high_disp_in_module0,
          " / ", n_high_disp, " high-dispersion genes fall in a real module (",
          n_high_disp_in_module0, " in module 0, excluded).")

  module_summary <- do.call(rbind, lapply(real_modules, function(m) {
    genes_in_module <- real_module_sub$gene[real_module_sub$module == m]
    disp_values <- dispersion_lookup[genes_in_module]
    data.frame(
      cell_type = ct, module = m,
      module_size = length(genes_in_module),
      median_dispersion = median(disp_values),
      n_high_dispersion = sum(genes_in_module %in% high_disp_genes),
      stringsAsFactors = FALSE
    )
  }))
  module_summary$pct_high_dispersion <- 100 * module_summary$n_high_dispersion / module_summary$module_size
  module_summary$module_rank <- rank(module_summary$median_dispersion, ties.method = "first")
  if (anyDuplicated(module_summary$median_dispersion)) {
    log_msg("  NOTE: tied median_dispersion values among modules -- rank order among ties is arbitrary (ties.method='first').")
  }
  summary_rows[[ct]] <- module_summary

  detail_rows[[ct]] <- data.frame(
    gene = real_module_sub$gene, cell_type = ct,
    module = real_module_sub$module,
    median_dispersion_across_individuals = dispersion_lookup[real_module_sub$gene],
    high_dispersion = real_module_sub$gene %in% high_disp_genes,
    stringsAsFactors = FALSE
  )
}

n_expected_celltypes <- length(CELL_TYPES)

detail_df <- do.call(rbind, detail_rows); rownames(detail_df) <- NULL
stopifnot(length(unique(detail_df$cell_type)) == n_expected_celltypes)
detail_csv <- file.path(OUT_DIR, "module_dispersion_gene_detail_consensus.csv")
write.csv(detail_df, detail_csv, row.names = FALSE)
log_msg("Wrote: ", detail_csv, " (", nrow(detail_df), " gene x cell-type rows, module 0 excluded)")

summary_df <- do.call(rbind, summary_rows); rownames(summary_df) <- NULL
stopifnot(length(unique(summary_df$cell_type)) == n_expected_celltypes)
summary_csv <- file.path(OUT_DIR, "module_dispersion_summary_consensus.csv")
write.csv(summary_df, summary_csv, row.names = FALSE)
log_msg("Wrote: ", summary_csv, " (", nrow(summary_df), " module x cell-type rows)")

plot_df <- summary_df
plot_df$module_rank <- ave(plot_df$module_rank, plot_df$cell_type, FUN = function(x) rank(x, ties.method = "first"))
label_nudge <- 0.04 * diff(range(plot_df$median_dispersion, na.rm = TRUE))
plot_df$label_y <- plot_df$median_dispersion + sign(plot_df$median_dispersion) * label_nudge
plot_df$label_count <- as.character(plot_df$n_high_dispersion)
plot_df$label_pct <- sprintf("%.0f%%", plot_df$pct_high_dispersion)

LABEL_SIZE <- 1.7 * 1.5

make_consensus_panel <- function(ct_df, ct_label, fill_color, label_col, panel_title) {
  ggplot(ct_df, aes(x = module_rank, y = median_dispersion)) +
    geom_col(width = 0.7, fill = fill_color, color = "grey20", linewidth = 0.15) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey40") +
    geom_text(aes(y = label_y, label = .data[[label_col]]), size = LABEL_SIZE, family = "Arial") +
    labs(x = "Module rank (ascending median dispersion)",
         y = "Median dispersion (res_var) per consensus module",
         title = panel_title) +
    theme_classic(base_family = "Arial") +
    theme(
      text = element_text(size = 13), axis.text = element_text(size = 11),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.title = element_text(size = 13),
      plot.title = element_text(size = 14, face = "bold"),
      axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
      plot.margin = margin(4, 4, 4, 4, unit = "mm"),
      panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.3)
    )
}

out_pdf <- file.path(OUT_DIR, "module_dispersion_profile_consensus.pdf")
cairo_pdf(out_pdf, width = 11, height = 15, onefile = TRUE)
for (ct in CELL_TYPES) {
  ct_label <- gsub("_", " ", ct)
  ct_df <- plot_df[plot_df$cell_type == ct, ]
  stopifnot("cell type missing from plot_df" = nrow(ct_df) > 0)
  fill_color <- celltype_colors[[ct_label]]

  p_count <- make_consensus_panel(ct_df, ct_label, fill_color, "label_count",
                                   "Raw count of top-quintile-dispersion genes per module")
  p_pct   <- make_consensus_panel(ct_df, ct_label, fill_color, "label_pct",
                                   "Proportion of top-quintile-dispersion genes per module")
  gridExtra::grid.arrange(
    p_count, p_pct, nrow = 2,
    top = grid::textGrob(
      paste0(ct_label, " -- consensus module dispersion profile"),
      gp = grid::gpar(fontsize = 17, fontface = "bold", fontfamily = "Arial")
    )
  )
}
dev.off()
log_msg("Wrote: ", out_pdf, " (", length(CELL_TYPES), " pages, one per cell type, count+proportion stacked)")

log_msg("DONE.")
