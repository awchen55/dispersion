#!/usr/bin/env Rscript

stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(ggplot2)
})

CONSENSUS_MODULE_DIR <- "figure_reproduction/outputs/coexpression_plasticity_lowfidelity/consensus_modules_prototype"
OUT_DIR         <- "figure_reproduction/outputs/module_dispersion_profile/phase2_consensus"
SUMMARY_CSV     <- file.path(OUT_DIR, "module_dispersion_summary_consensus.csv")
MARKER_CSV      <- file.path(OUT_DIR, "marker_module_lookup_consensus.csv")
CC_COMPOSITION_CSV <- file.path(OUT_DIR, "cc_module_phase_composition_consensus.csv")
for (f in c(SUMMARY_CSV, MARKER_CSV, CC_COMPOSITION_CSV)) {
  if (!file.exists(f)) stop("Missing ", f, " -- run scripts 02, 07, and 12 first.")
}

CELL_TYPES <- c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                "Epithelial_cells", "Activated_fibroblasts",
                "Foregut_endodermal_cells")
CONSENSUS_QUANTILE_SUFFIX <- "q0.5"
LARGE_MODULE_THRESHOLD <- 10L

# ANALYSIS_OK[hardcoded-gene-list]: literal HGNC gene symbols, identical to
# the established repo-wide source.
CHAMBER_SPECIFICITY <- c("MYL2", "MYL7", "HEY2", "IRX4", "TBX5", "NR2F2")

log_msg <- function(...) message("[14_module_dispersion_profile_consensus_combined] ", ...)

OTHER_MODULE_COLOR    <- "grey30"
CC_MODULE_COLOR       <- "#009E73"
CELLTYPE_MARKER_COLOR <- "#D55E00"
CHAMBER_MARKER_COLOR  <- "#CC79A7"

summary_df <- read.csv(SUMMARY_CSV, stringsAsFactors = FALSE)
stopifnot(
  "module_dispersion_summary_consensus.csv missing expected columns" =
    all(c("cell_type", "module", "module_size", "median_dispersion", "module_rank",
          "n_high_dispersion", "pct_high_dispersion") %in% colnames(summary_df))
)

marker_lookup <- read.csv(MARKER_CSV, stringsAsFactors = FALSE)
stopifnot("marker_module_lookup_consensus.csv missing expected columns" =
            all(c("cell_type", "gene", "module", "in_real_module") %in% colnames(marker_lookup)))
marker_lookup <- marker_lookup[marker_lookup$in_real_module, , drop = FALSE]

cc_composition <- read.csv(CC_COMPOSITION_CSV, stringsAsFactors = FALSE)
stopifnot("cc_module_phase_composition_consensus.csv missing expected columns" =
            all(c("cell_type", "module", "dominant_phase", "n_total_cc", "purity",
                  "is_large_module") %in% colnames(cc_composition)))

cc_module_df <- cc_composition[cc_composition$is_large_module, ]
overrides <- data.frame(
  cell_type = c("Cardiomyocytes", "Activated_fibroblasts"),
  module = c(11, 4),
  expected_phase = c("S", "S"),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(overrides))) {
  ov <- cc_composition[cc_composition$cell_type == overrides$cell_type[i] &
                          cc_composition$module == overrides$module[i], ]
  stopifnot(
    "override module not found in cc composition data -- investigate" = nrow(ov) == 1,
    "override module's dominant phase doesn't match expectation -- investigate before including it" =
      ov$dominant_phase == overrides$expected_phase[i]
  )
  cc_module_df <- rbind(cc_module_df, ov)
}
cc_module_df$category <- "cell_cycle_module"
cc_module_df$driving_detail <- sprintf("%s-dominant cc module, %d cc genes, purity=%.2f",
                                        cc_module_df$dominant_phase, cc_module_df$n_total_cc, cc_module_df$purity)
log_msg(nrow(cc_module_df), " cell-cycle module(s) selected across ",
        length(unique(cc_module_df$cell_type)), " cell types (Cardiac_fibroblasts has none).")

celltype_marker_df <- aggregate(gene ~ cell_type + module, data = marker_lookup,
                                 FUN = function(x) paste(sort(x), collapse = ", "))
celltype_marker_df$category <- "celltype_marker"
celltype_marker_df$driving_detail <- paste0("cell-type marker(s): ", celltype_marker_df$gene)
celltype_marker_df$gene <- NULL
log_msg(nrow(celltype_marker_df), " cell-type-marker module(s) found across ",
        length(unique(celltype_marker_df$cell_type)), " cell types.")

chamber_rows <- list()
for (ct in CELL_TYPES) {
  module_csv <- file.path(CONSENSUS_MODULE_DIR, paste0(ct, "_consensus_", CONSENSUS_QUANTILE_SUFFIX, "_module_membership.csv"))
  if (!file.exists(module_csv)) stop("Missing consensus module membership file: ", module_csv)
  module_df <- read.csv(module_csv, stringsAsFactors = FALSE)
  stopifnot("consensus module_membership.csv missing expected columns" =
              all(c("gene", "module", "in_module") %in% colnames(module_df)))

  found <- intersect(CHAMBER_SPECIFICITY, module_df$gene)
  n_missing <- length(CHAMBER_SPECIFICITY) - length(found)
  if (n_missing > 0) {
    log_msg("  ", ct, ": ", n_missing, " / ", length(CHAMBER_SPECIFICITY),
            " chamber-specificity marker(s) not in the consensus gene universe: ",
            paste(setdiff(CHAMBER_SPECIFICITY, module_df$gene), collapse = ", "))
  }
  for (g in found) {
    m <- module_df$module[module_df$gene == g]
    in_mod <- module_df$in_module[module_df$gene == g]
    if (m == 0 || !in_mod) {
      log_msg("  ", ct, ": chamber marker ", g, " is in consensus module 0 (not a real module) -- not shown on the plot.")
      next
    }
    chamber_rows[[paste0(ct, "_", g)]] <- data.frame(cell_type = ct, module = m, gene = g, stringsAsFactors = FALSE)
  }
}
chamber_lookup <- do.call(rbind, chamber_rows); rownames(chamber_lookup) <- NULL
chamber_marker_df <- aggregate(gene ~ cell_type + module, data = chamber_lookup,
                                FUN = function(x) paste(sort(x), collapse = ", "))
chamber_marker_df$category <- "chamber_marker"
chamber_marker_df$driving_detail <- paste0("chamber-specificity marker(s): ", chamber_marker_df$gene)
chamber_marker_df$gene <- NULL
log_msg(nrow(chamber_marker_df), " chamber-marker module(s) found across ",
        length(unique(chamber_marker_df$cell_type)), " cell types.")

cols <- c("cell_type", "module", "category", "driving_detail")
category_df <- rbind(cc_module_df[, cols], celltype_marker_df[, cols], chamber_marker_df[, cols])
category_df$priority <- match(category_df$category, c("cell_cycle_module", "celltype_marker", "chamber_marker"))
dup_key <- paste(category_df$cell_type, category_df$module)
n_overlap_rows <- sum(duplicated(dup_key) | duplicated(dup_key, fromLast = TRUE))
if (n_overlap_rows > 0) {
  overlap_detail <- category_df[dup_key %in% dup_key[duplicated(dup_key)], ]
  overlap_detail <- overlap_detail[order(overlap_detail$cell_type, overlap_detail$module), ]
  log_msg(nrow(overlap_detail), " module x cell-type rows involved in a multi-category overlap ",
          "(resolved by priority cell_cycle_module > celltype_marker > chamber_marker): ",
          paste(sprintf("%s module %d [%s]", overlap_detail$cell_type, overlap_detail$module, overlap_detail$category), collapse = "; "))
}
category_df <- category_df[order(category_df$priority), ]
category_df <- category_df[!duplicated(dup_key[order(category_df$priority)]), ]
category_df$priority <- NULL

category_csv <- file.path(OUT_DIR, "module_category_lookup_consensus.csv")
write.csv(category_df, category_csv, row.names = FALSE)
log_msg("Wrote: ", category_csv, " (", nrow(category_df),
        " non-'other' module x cell-type rows, one category each after priority resolution)")

plot_df <- summary_df
plot_df$module_rank <- ave(plot_df$module_rank, plot_df$cell_type, FUN = function(x) rank(x, ties.method = "first"))
plot_df <- merge(plot_df, category_df[, c("cell_type", "module", "category")],
                  by = c("cell_type", "module"), all.x = TRUE)
plot_df$category[is.na(plot_df$category)] <- "other"
label_nudge <- 0.04 * diff(range(plot_df$median_dispersion, na.rm = TRUE))
plot_df$count_label_y <- plot_df$median_dispersion + sign(plot_df$median_dispersion) * label_nudge
plot_df$count_label <- as.character(plot_df$n_high_dispersion)
proportion_str <- sprintf("%.2f", plot_df$pct_high_dispersion / 100)
plot_df$pct_label <- ifelse(substr(proportion_str, 1, 2) == "0.", substring(proportion_str, 2), proportion_str)

LABEL_SIZE <- 1.7 * 1.5

out_pdf <- file.path(OUT_DIR, "module_dispersion_profile_consensus_combined.pdf")
cairo_pdf(out_pdf, width = 11, height = 8.5, onefile = TRUE)
for (ct in CELL_TYPES) {
  ct_label <- gsub("_", " ", ct)
  ct_df <- plot_df[plot_df$cell_type == ct, ]
  stopifnot("cell type missing from plot_df" = nrow(ct_df) > 0)
  ct_df$category <- factor(ct_df$category,
                            levels = c("other", "cell_cycle_module", "celltype_marker", "chamber_marker"))
  ct_range <- range(ct_df$median_dispersion, na.rm = TRUE)
  ct_df$pct_label_y <- ct_range[1] - 0.14 * diff(ct_range)

  fill_values <- c(other = OTHER_MODULE_COLOR, cell_cycle_module = CC_MODULE_COLOR,
                    celltype_marker = CELLTYPE_MARKER_COLOR, chamber_marker = CHAMBER_MARKER_COLOR)
  fill_labels <- c(other = "Other module", cell_cycle_module = "Cell cycle module (F-060, large S/G2M)",
                    celltype_marker = "Cell-type marker gene", chamber_marker = "Chamber-specificity marker gene")

  p <- ggplot(ct_df, aes(x = module_rank, y = median_dispersion, fill = category)) +
    geom_col(width = 0.7, color = "grey20", linewidth = 0.15) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey40") +
    geom_text(aes(y = count_label_y, label = count_label), size = LABEL_SIZE, family = "Arial") +
    geom_text(aes(y = pct_label_y, label = pct_label), size = LABEL_SIZE, family = "Arial",
              angle = 90, hjust = 0) +
    scale_fill_manual(values = fill_values, labels = fill_labels, name = NULL, drop = TRUE) +
    labs(x = "Module rank (ascending median dispersion)",
         y = "Median dispersion",
         title = ct_label) +
    theme_classic(base_family = "Arial") +
    theme(
      text = element_text(size = 13), axis.text = element_text(size = 11),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.title = element_text(size = 13),
      plot.title = element_text(size = 16, face = "bold"),
      legend.position = "bottom",
      axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
      plot.margin = margin(4, 4, 4, 4, unit = "mm"),
      panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.3)
    )
  print(p)
}
dev.off()
log_msg("Wrote: ", out_pdf, " (", length(CELL_TYPES), " pages, one per cell type)")

log_msg("DONE.")
