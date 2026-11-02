#!/usr/bin/env Rscript

stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(png)
})

BASE_DIR <- "figure_reproduction/outputs/coexpression_plasticity_lowfidelity"
QC_DIR   <- file.path(BASE_DIR, "consensus_module_qc")
R2_CSV   <- file.path(BASE_DIR, "diagnostics", "effect_size_sweep_all_celltypes.csv")
DENDRO_PNG <- file.path(QC_DIR, "consensus_module_dendrograms.png")
CNT_CSV  <- file.path(QC_DIR, "consensus_module_count.csv")
for (f in c(R2_CSV, DENDRO_PNG, CNT_CSV)) if (!file.exists(f)) stop("Missing ", f)

CELL_TYPES <- c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells")
INDIVIDUALS <- c("NA19099", "NA19093", "NA18522")
FLOOR <- 0.25
RSQUARED_CUT <- 0.85
R2 <- "R²"
ct_label <- gsub("_", " ", CELL_TYPES)

log_msg <- function(...) message("[11_consensus_qc_combined_figure] ", ...)

celltype_colors <- c(
  "Activated fibroblasts" = "#117733", "Cardiac fibroblasts" = "#4477AA",
  "Cardiomyocytes" = "#882255", "Epithelial cells" = "#AA4499",
  "Foregut endodermal cells" = "#CC6677", "Mesenchymal cells" = "#DDCC77"
)
base_theme <- theme_classic(base_family = "Arial") +
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        plot.title = element_text(size = 11, face = "bold"),
        plot.margin = margin(3, 3, 3, 3, unit = "mm"))

# --- A: R^2 at chosen power, floor 0.25 -------------------------------------
r2 <- read.csv(R2_CSV, stringsAsFactors = FALSE)
r2 <- r2[abs(r2$floor - FLOOR) < 1e-9, ]
stopifnot("expected 18 rows at floor 0.25" = nrow(r2) == length(CELL_TYPES) * length(INDIVIDUALS),
          "unexpected cell type" = all(r2$cell_type %in% CELL_TYPES))
r2$cell_type_label <- factor(gsub("_", " ", r2$cell_type), levels = ct_label)
r2$sample_label <- factor(paste0(r2$cell_type_label, "\n", r2$individual),
                          levels = paste0(rep(ct_label, each = 3), "\n", rep(INDIVIDUALS, times = 6)))
p_a <- ggplot(r2, aes(x = sample_label, y = chosen_r2, fill = cell_type_label)) +
  geom_col(color = "black", linewidth = 0.3) +
  geom_hline(yintercept = RSQUARED_CUT, linetype = "dashed", linewidth = 0.3, color = "grey30") +
  geom_text(aes(label = sprintf("p=%d%s", chosen_power, ifelse(crosses_at_all, "", "*"))),
            vjust = -0.4, size = 2.6, family = "Arial") +
  scale_fill_manual(values = celltype_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 1.12), breaks = c(0, 0.5, 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = paste0("Scale-free fit (", R2, ")"),
       title = paste0(R2, " at the chosen power, per sample (floor = 0.25)"),
       subtitle = paste0("p = chosen power; * = never reached ", R2, "=", RSQUARED_CUT,
                         " (global-best shown). Dashed line: ", RSQUARED_CUT, ".")) +
  base_theme +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1), plot.subtitle = element_text(size = 8))

# --- B: gene dendrograms with module colour bars (rendered by script 12) ---------
dendro_img <- png::readPNG(DENDRO_PNG)
stopifnot("dendrogram PNG is empty" = length(dim(dendro_img)) == 3 && all(dim(dendro_img)[1:2] > 100))
p_b <- ggdraw() +
  draw_grob(grid::rasterGrob(dendro_img, x = 0.5, y = 0.465, width = 1, height = 0.93, interpolate = TRUE)) +
  draw_label("Consensus module gene dendrograms", x = 0.02, y = 0.97, hjust = 0, vjust = 1,
             fontface = "bold", size = 11)

# --- C: modules per cell type -------------------------------------------------
cnt <- read.csv(CNT_CSV, stringsAsFactors = FALSE)
stopifnot("unexpected cell type in count CSV" = all(cnt$cell_type %in% CELL_TYPES), nrow(cnt) == length(CELL_TYPES))
cnt$cell_type_label <- factor(gsub("_", " ", cnt$cell_type), levels = ct_label)
p_c <- ggplot(cnt, aes(x = cell_type_label, y = n_modules, fill = cell_type_label)) +
  geom_col(color = "black", linewidth = 0.3) +
  geom_text(aes(label = n_modules), vjust = -0.4, size = 3, family = "Arial") +
  scale_fill_manual(values = celltype_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Number of consensus modules", title = "Consensus modules per cell type") +
  base_theme + theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

bottom <- plot_grid(p_b, p_c, nrow = 1, rel_widths = c(2, 1), labels = c("B", "C"), label_size = 12, align = "h", axis = "tb")
fig <- plot_grid(p_a, bottom, ncol = 1, rel_heights = c(1, 1.5), labels = c("A", ""), label_size = 12)

out <- file.path(QC_DIR, "consensus_qc_combined.pdf")
ggsave(out, fig, width = 12, height = 11, units = "in", device = cairo_pdf)
log_msg("Wrote: ", out, "\n[11_consensus_qc_combined_figure] DONE.")
