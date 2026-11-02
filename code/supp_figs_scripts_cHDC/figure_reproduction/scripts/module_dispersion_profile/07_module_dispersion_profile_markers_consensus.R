#!/usr/bin/env Rscript

stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(ggplot2)
})

CONSENSUS_DIR <- "figure_reproduction/outputs/coexpression_plasticity_lowfidelity/consensus_modules_prototype"
OUT_DIR       <- "figure_reproduction/outputs/module_dispersion_profile/phase2_consensus"
SUMMARY_CSV   <- file.path(OUT_DIR, "module_dispersion_summary_consensus.csv")
if (!file.exists(SUMMARY_CSV)) stop("Missing ", SUMMARY_CSV, " -- run script 02 first.")

CELL_TYPES <- c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                "Epithelial_cells", "Activated_fibroblasts",
                "Foregut_endodermal_cells")
CONSENSUS_QUANTILE_SUFFIX <- "q0.5"

# ANALYSIS_OK[hardcoded-gene-list]: literal HGNC gene symbols, identical to
# the established source (coexpression_plasticity_lowfidelity/scripts/16),
# not sample/individual IDs.
CELLTYPE_MARKERS <- list(
  Cardiomyocytes           = c("MYH6", "NKX2-5", "TNNT2", "ACTC1", "RYR2"),
  Mesenchymal_cells        = c("TWIST1", "PRRX2", "FOXC1"),
  Cardiac_fibroblasts      = c("TCF21", "ITGA8", "PDGFRA", "COL1A1", "COL3A1"),
  Activated_fibroblasts    = c("POSTN", "DCN", "LUM", "PDPN"),
  Epithelial_cells         = c("WFDC2", "GRHL2", "CLDN6"),
  Foregut_endodermal_cells = c("HHEX", "ALDH1A1", "ONECUT1")
)
stopifnot(setequal(names(CELLTYPE_MARKERS), CELL_TYPES))

celltype_colors <- c(
  "Activated fibroblasts"    = "#117733", "Cardiac fibroblasts" = "#4477AA",
  "Cardiomyocytes"           = "#882255", "Epithelial cells"    = "#AA4499",
  "Foregut endodermal cells" = "#CC6677", "Mesenchymal cells"   = "#DDCC77"
)
MARKER_HIGHLIGHT_COLOR <- "#D55E00"

log_msg <- function(...) message("[07_module_dispersion_profile_markers_consensus] ", ...)

summary_df <- read.csv(SUMMARY_CSV, stringsAsFactors = FALSE)
stopifnot(
  "module_dispersion_summary_consensus.csv missing expected columns" =
    all(c("cell_type", "module", "module_size", "median_dispersion",
          "module_rank") %in% colnames(summary_df))
)

lookup_rows <- list()

for (ct in CELL_TYPES) {
  markers <- CELLTYPE_MARKERS[[ct]]
  module_csv <- file.path(CONSENSUS_DIR, paste0(ct, "_consensus_", CONSENSUS_QUANTILE_SUFFIX, "_module_membership.csv"))
  if (!file.exists(module_csv)) stop("Missing consensus module membership file: ", module_csv)
  module_df <- read.csv(module_csv, stringsAsFactors = FALSE)
  stopifnot("consensus module_membership.csv missing expected columns" =
              all(c("gene", "module", "in_module") %in% colnames(module_df)))

  markers_found <- intersect(markers, module_df$gene)
  n_missing <- length(markers) - length(markers_found)
  if (n_missing > 0) {
    log_msg("  ", ct, ": ", n_missing, " / ", length(markers),
            " marker(s) not in the consensus gene universe: ",
            paste(setdiff(markers, module_df$gene), collapse = ", "))
  }

  for (g in markers_found) {
    m <- module_df$module[module_df$gene == g]
    in_mod <- module_df$in_module[module_df$gene == g]
    lookup_rows[[paste0(ct, "_", g)]] <- data.frame(
      cell_type = ct, gene = g, module = m,
      in_real_module = (m != 0 & in_mod), stringsAsFactors = FALSE
    )
    if (m == 0 || !in_mod) {
      log_msg("  ", ct, ": marker ", g, " is in consensus module 0 (not a real module) -- not shown on the plot.")
    }
  }
}

lookup_df <- do.call(rbind, lookup_rows); rownames(lookup_df) <- NULL
lookup_csv <- file.path(OUT_DIR, "marker_module_lookup_consensus.csv")
write.csv(lookup_df, lookup_csv, row.names = FALSE)
log_msg("Wrote: ", lookup_csv, " (", nrow(lookup_df), " marker x cell-type rows found in the universe)")

real_lookup <- lookup_df[lookup_df$in_real_module, , drop = FALSE]
marker_labels <- aggregate(gene ~ cell_type + module, data = real_lookup,
                            FUN = function(x) paste(sort(x), collapse = ",\n"))

plot_df <- summary_df
plot_df$module_rank <- ave(plot_df$module_rank, plot_df$cell_type, FUN = function(x) rank(x, ties.method = "first"))
plot_df <- merge(plot_df, marker_labels, by = c("cell_type", "module"), all.x = TRUE)
plot_df$has_marker <- !is.na(plot_df$gene)
label_nudge <- 0.04 * diff(range(plot_df$median_dispersion, na.rm = TRUE))
plot_df$label_y <- plot_df$median_dispersion + sign(plot_df$median_dispersion) * label_nudge

out_pdf <- file.path(OUT_DIR, "module_dispersion_profile_markers_consensus.pdf")
cairo_pdf(out_pdf, width = 11, height = 8.5, onefile = TRUE)
for (ct in CELL_TYPES) {
  ct_label <- gsub("_", " ", ct)
  ct_df <- plot_df[plot_df$cell_type == ct, ]
  stopifnot("cell type missing from plot_df" = nrow(ct_df) > 0)

  p <- ggplot(ct_df, aes(x = module_rank, y = median_dispersion, fill = has_marker)) +
    geom_col(width = 0.7, color = "grey20", linewidth = 0.15) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey40") +
    geom_text(data = ct_df[ct_df$has_marker, ], aes(y = label_y, label = gene),
              size = 2.55, family = "Arial", lineheight = 0.85) +
    scale_fill_manual(values = c(`TRUE` = MARKER_HIGHLIGHT_COLOR, `FALSE` = celltype_colors[[ct_label]]), guide = "none") +
    labs(x = "Module rank (ascending median dispersion)", y = "Median dispersion (res_var) per consensus module",
         title = paste0(ct_label, " -- consensus modules containing its own marker genes"),
         subtitle = paste0("Module 0 excluded; markers: ", paste(CELLTYPE_MARKERS[[ct]], collapse = ", "))) +
    theme_classic(base_family = "Arial") +
    theme(
      text = element_text(size = 13), axis.text = element_text(size = 11),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.title = element_text(size = 13),
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 9),
      axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
      plot.margin = margin(4, 4, 4, 4, unit = "mm"),
      panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.3)
    )
  print(p)
}
dev.off()
log_msg("Wrote: ", out_pdf, " (", length(CELL_TYPES), " pages, one per cell type)")

log_msg("DONE.")
