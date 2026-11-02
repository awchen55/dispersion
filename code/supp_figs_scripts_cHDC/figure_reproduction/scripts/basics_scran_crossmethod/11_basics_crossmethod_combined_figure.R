#!/usr/bin/env Rscript

stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

OUT_DIR <- "figure_reproduction/outputs/basics_scran_crossmethod"
RHO_CSV <- file.path(OUT_DIR, "comparison/celltype_specificity/rho_matrix_matched_1000cell.csv")
RET_CSV <- file.path(OUT_DIR, "consistent_geneset_retention/retention_summary.csv")
PAIR_CSV <- file.path(OUT_DIR, "crosscheck_memento/matched_1000cell/pairwise_crosscheck_summary.csv")
OUT_PDF <- file.path(OUT_DIR, "basics_crossmethod_combined.pdf")
N_PAIRS <- choose(6, 2)  # 6 cell types -> 15 pairwise contrasts
for (f in c(RHO_CSV, RET_CSV, PAIR_CSV)) if (!file.exists(f)) stop("Missing ", f)

# --- A: heatmap (same scale/limits as script 08) -----------------------------
rho_mat <- as.matrix(read.csv(RHO_CSV, row.names = 1, check.names = FALSE))
stopifnot("rho matrix must be square 6x6" = all(dim(rho_mat) == 6L),
          "rho matrix rows/cols must match" = identical(rownames(rho_mat), colnames(rho_mat)))
hm_df <- as.data.frame(as.table(rho_mat))
colnames(hm_df) <- c("basics_cell_type", "memento_cell_type", "rho")
hm_df$is_matched <- hm_df$basics_cell_type == hm_df$memento_cell_type
p_a <- ggplot(hm_df, aes(memento_cell_type, basics_cell_type, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", rho), fontface = ifelse(is_matched, "bold", "plain")), size = 3) +
  scale_fill_gradient(low = "white", high = "red", limits = c(0.8, 1), oob = scales::squish) +
  labs(x = "Memento cell type", y = "BASiCS cell type", fill = "Spearman\nrho") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom", legend.key.width = grid::unit(1.2, "cm"))

# --- B: retention bars, scran dropped ----------------------------------------
ret <- read.csv(RET_CSV, stringsAsFactors = FALSE)
ret <- ret[ret$method != "scran_fullcell", ]
stopifnot("expected basics + memento-subset rows x low/high" = nrow(ret) == 4L,
          "scran rows not removed" = !any(grepl("scran", ret$method)))
ret$geneset <- factor(ret$geneset, levels = c("consistently_low", "consistently_high"),
                      labels = c("Consistently low", "Consistently high"))
p_b <- ggplot(ret, aes(method, prop_retained_of_memento_geneset, fill = geneset)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", prop_retained_of_memento_geneset)),
            position = position_dodge(width = 0.7), vjust = -0.4, size = 3) +
  labs(x = NULL, y = "Proportion of Memento's gene set retained", fill = NULL) +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "bottom")

# --- C: per-pair logFC Spearman rho, Memento-significant genes ---------------
pair_df <- read.csv(PAIR_CSV, stringsAsFactors = FALSE)
stopifnot("expected 15 pairwise rows" = nrow(pair_df) == N_PAIRS)
rho_df <- data.frame(
  pair = gsub("_", " ", paste0(pair_df$cell_type_1, " vs ", pair_df$cell_type_2)),
  rho = pair_df$logfc_spearman_rho_sig_memento,
  lo = pair_df$logfc_spearman_ci_lo_sig_memento,
  hi = pair_df$logfc_spearman_ci_hi_sig_memento)
stopifnot("non-finite rho/CI in pairwise summary" = all(is.finite(c(rho_df$rho, rho_df$lo, rho_df$hi))))
# ANALYSIS_OK[plot-ordering]: re-levels the factor for display order only -- same 15 pairs/values.
rho_df$pair <- factor(rho_df$pair, levels = rho_df$pair[order(rho_df$rho)])
p_c <- ggplot(rho_df, aes(x = rho, y = pair)) +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.4) +
  labs(x = "Spearman rho (BASiCS ResDispDistance vs. Memento logFC)", y = NULL) +
  theme_bw()

# Layout: A (heatmap) and B (retention bars) on the top row; C (per-pair logFC
# rho) spans the full width underneath.
# wrap_elements(full = ) keeps C out of patchwork's panel alignment with A/B;
# otherwise A's plot area is pushed right to line up with C's (long y labels).
combined <- ((p_a + p_b) + plot_layout(widths = c(1.15, 1))) / wrap_elements(full = p_c) +
  plot_layout(heights = c(1, 1.1)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))
ggsave(OUT_PDF, combined, width = 11, height = 9.5, device = cairo_pdf)
message("[11_basics_crossmethod_combined_figure] wrote ", OUT_PDF)
