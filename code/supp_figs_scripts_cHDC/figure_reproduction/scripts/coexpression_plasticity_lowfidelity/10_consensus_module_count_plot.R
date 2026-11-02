#!/usr/bin/env Rscript

stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages(library(ggplot2))

BASE_DIR <- "figure_reproduction/outputs/coexpression_plasticity_lowfidelity"
IN_DIR   <- file.path(BASE_DIR, "consensus_modules_prototype")
OUT_DIR  <- file.path(BASE_DIR, "consensus_module_qc")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CELL_TYPES <- c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells")
log_msg <- function(...) message("[10_consensus_module_count_plot] ", ...)

rows <- lapply(CELL_TYPES, function(ct) {
  f <- file.path(IN_DIR, paste0(ct, "_consensus_q0.5_module_membership.csv"))
  if (!file.exists(f)) stop("Missing ", f)
  m <- read.csv(f, stringsAsFactors = FALSE)
  stopifnot("in_module disagrees with module==0" = all((m$module != 0) == m$in_module))
  data.frame(cell_type = ct, n_modules = length(unique(m$module[m$in_module])),
             n_genes = nrow(m), n_unassigned = sum(!m$in_module))
})
df <- do.call(rbind, rows)
stopifnot("a cell type has no modules" = all(df$n_modules > 0))
write.csv(df, file.path(OUT_DIR, "consensus_module_count.csv"), row.names = FALSE)
log_msg(paste(df$cell_type, df$n_modules, sep = "=", collapse = ", "))

df$cell_type_label <- factor(gsub("_", " ", df$cell_type), levels = gsub("_", " ", CELL_TYPES))

celltype_colors <- c(
  "Activated fibroblasts" = "#117733", "Cardiac fibroblasts" = "#4477AA",
  "Cardiomyocytes" = "#882255", "Epithelial cells" = "#AA4499",
  "Foregut endodermal cells" = "#CC6677", "Mesenchymal cells" = "#DDCC77"
)

p <- ggplot(df, aes(x = cell_type_label, y = n_modules, fill = cell_type_label)) +
  geom_col(color = "black", linewidth = 0.3) +
  geom_text(aes(label = n_modules), vjust = -0.4, size = 4.3, family = "Arial") +
  scale_fill_manual(values = celltype_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Number of consensus modules",
       title = "Consensus modules per cell type",
       subtitle = "Unassigned genes (module 0) not counted. Epithelial cells: sub-cluster 24 excluded.") +
  theme_classic(base_family = "Arial") +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 11.7),
    plot.title = element_text(size = 15, face = "bold"), plot.subtitle = element_text(size = 10.8),
    plot.margin = margin(4, 4, 4, 4, unit = "mm")
  )

out <- file.path(OUT_DIR, "consensus_module_count.pdf")
ggsave(out, p, width = 7, height = 5, units = "in", device = cairo_pdf)
log_msg("Wrote: ", out, "\n[10_consensus_module_count_plot] DONE.")
