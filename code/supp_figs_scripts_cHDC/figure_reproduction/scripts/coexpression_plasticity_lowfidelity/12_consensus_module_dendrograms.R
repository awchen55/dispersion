#!/usr/bin/env Rscript

stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(WGCNA)
  library(ragg)
})

BASE_DIR <- "figure_reproduction/outputs/coexpression_plasticity_lowfidelity"
TOM_DIR  <- file.path(BASE_DIR, "consensus_tom")
MEMB_DIR <- file.path(BASE_DIR, "consensus_modules_prototype")
OUT_PDF  <- file.path(BASE_DIR, "consensus_module_qc", "consensus_module_dendrograms.pdf")
OUT_PNG  <- file.path(BASE_DIR, "consensus_module_qc", "consensus_module_dendrograms.png")
dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)

CELL_TYPES <- c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells")

log_msg <- function(...) message("[12_consensus_module_dendrograms] ", ...)

# Build every tree and colour vector first, so any missing input fails before
# the PDF device is opened.
panels <- list()
for (ct in CELL_TYPES) {
  tom_rds <- file.path(TOM_DIR, paste0(ct, "_consensus_tom.rds"))
  memb_csv <- file.path(MEMB_DIR, paste0(ct, "_consensus_q0.5_module_membership.csv"))
  for (f in c(tom_rds, memb_csv)) if (!file.exists(f)) stop("Missing ", f, " -- run script 08 first.")
  tom <- readRDS(tom_rds)
  genes <- rownames(tom)
  memb_df <- read.csv(memb_csv, stringsAsFactors = FALSE)
  if (!setequal(memb_df$gene, genes)) stop(ct, ": membership genes differ from the saved TOM's genes")
  memb <- setNames(memb_df$module, memb_df$gene)[genes]
  stopifnot("unassigned (module 0) must be the only non-module label" = all(memb >= 0))
  tree <- stats::hclust(stats::as.dist(1 - tom), method = "average")
  n_mod <- length(unique(memb[memb != 0]))
  cols <- WGCNA::labels2colors(memb)
  stopifnot("module-0 genes must be grey" = all(cols[memb == 0] == "grey"))
  panels[[ct]] <- list(tree = tree, colors = cols, n_mod = n_mod, n_genes = length(genes),
                       n_unassigned = sum(memb == 0))
  log_msg(ct, ": ", length(genes), " genes, ", n_mod, " modules, ", sum(memb == 0), " unassigned")
  rm(tom); gc(FALSE)
}

# Layout: 2 x 3 grid of panels; each panel = dendrogram region + colour-bar
# region, numbered in plot order (dendrogram then colours, panel by panel).
layout_mat <- rbind(c(1, 3, 5), c(2, 4, 6), c(7, 9, 11), c(8, 10, 12))

draw_panels <- function(cex_main, cex_axis, cex_lab, cex_color, two_line_title) {
  layout(layout_mat, heights = c(3, 0.5, 3, 0.5))
  for (ct in CELL_TYPES) {
    p <- panels[[ct]]
    ttl <- if (two_line_title) {
      # Combined-figure version: title is the cell type only (module/unassigned
      # counts dropped -- module counts are in panel C), broken after its first
      # word onto two rows; single-word names get a leading blank row so every
      # title is bottom-aligned over its dendrogram.
      sub(" ", "\n", if (grepl("_", ct)) gsub("_", " ", ct) else paste0(" ", ct))
    } else {
      sprintf("%s (%d modules; %d/%d genes unassigned)", gsub("_", " ", ct), p$n_mod, p$n_unassigned, p$n_genes)
    }
    WGCNA::plotDendroAndColors(
      p$tree, p$colors, groupLabels = "Module", dendroLabels = FALSE, hang = 0.03,
      addGuide = FALSE, guideHang = 0.05, setLayout = FALSE,
      marAll = c(0.5, if (two_line_title) 5 else 4, if (two_line_title) 3.8 else 2.5, 0.5), cex.colorLabels = cex_color,
      main = ttl, cex.main = cex_main, cex.axis = cex_axis,
      ylab = "1 - consensus TOM", cex.lab = cex_lab
    )
  }
}

grDevices::cairo_pdf(OUT_PDF, width = 15, height = 9)
draw_panels(cex_main = 1.1, cex_axis = 1, cex_lab = 0.9, cex_color = 0.8, two_line_title = FALSE)
dev.off()
log_msg("Wrote: ", OUT_PDF)

# Sized for embedding as panel B of script 11's figure (8 x 5.6 in at 400 dpi).
ragg::agg_png(OUT_PNG, width = 8, height = 5.6, units = "in", res = 400, background = "white")
draw_panels(cex_main = 1.05, cex_axis = 1, cex_lab = 1, cex_color = 0.6, two_line_title = TRUE)
dev.off()
log_msg("Wrote: ", OUT_PNG, "\n[12_consensus_module_dendrograms] DONE.")
