#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(grid)
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
analysis_dir <- file.path(repo_dir, "figure_reproduction")
fig2_dir <- file.path(analysis_dir, "outputs", "fig2_replication")

# Smaller overall than script 34's (panel 1.5in, wrap 18, height 4in)
# throughout -- fewer terms represented (8 here vs. Fig2E's 10).
panel_width_in <- 1.2
plot_height_in <- 3.2
# Narrowest wrap width that still caps every label at 3 lines (matching
# Fig2E's own 3-line cap for visual consistency between the two panels).
label_wrap_width <- 16

# Byte-for-byte copy of the base replication script's helper.
make_go_dotplot <- function(plot_df, title) {
  ggplot(plot_df, aes(x = GeneRatio_num, y = Description)) +
    geom_point(aes(size = Count, color = log10_padj)) +
    scale_color_gradient(low = "grey85", high = "#8B0000",
                         name = "-log10(adjusted p)") +
    scale_size(range = c(2, 6)) +
    labs(title = title, x = "Gene ratio", y = NULL) +
    theme_classic(base_family = "Arial") +
    theme(
      text = element_text(size = 6),
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 6),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      plot.title = element_text(size = 6, hjust = 0.5),
      legend.key.height = unit(3, "mm"),
      legend.key.width  = unit(3, "mm"),
      plot.margin = margin(2, 2, 2, 2, unit = "mm")
    )
}

# ------------------------------- LOAD -----------------------------------------
top8_csv <- file.path(fig2_dir, "fig2a_go_enrichment_retained_top8.csv")
if (!file.exists(top8_csv)) stop("Missing ", top8_csv, " -- run 04_replicate_fig2_panels.R first")
plot_df <- read.csv(top8_csv) %>%
  mutate(Description = str_wrap(gsub("_", " ", sub("^GOBP_", "", ID)),
                                width = label_wrap_width)) %>%
  arrange(desc(GeneRatio_num)) %>%
  mutate(Description = factor(Description, levels = rev(Description)))

# Title removed (matches script 34/Fig3A's identical earlier removal).
p <- make_go_dotplot(plot_df, title = NULL)

# 10% x-scale expansion, same fix as script 34 -- default 5% left too
# little absolute margin for the largest points' fixed real-unit radius at
# a narrow panel, clipping them against the panel edge.
p <- p + scale_x_continuous(expand = expansion(mult = 0.10))

# ------------------------------- GTABLE SURGERY: FIX PANEL WIDTH -----------------
# Same two-pass measure-then-draw pattern as script 34: measure inside a
# throwaway cairo_pdf device, then draw on the real one at the now-correctly-
# measured size.
measure_pdf <- tempfile(fileext = ".pdf")
cairo_pdf(measure_pdf, width = 9, height = plot_height_in)
g <- ggplotGrob(p)
panel_idx <- which(g$layout$name == "panel")
stopifnot("Expected exactly one panel grob in the gtable layout" = length(panel_idx) == 1)
g$widths[g$layout$l[panel_idx]] <- unit(panel_width_in, "in")
total_width_in <- sum(convertWidth(g$widths, "in", valueOnly = TRUE))
dev.off()
unlink(measure_pdf)
message(sprintf("[panel-width] panel fixed to %.2fin; total canvas %.2fx%.2fin",
                panel_width_in, total_width_in, plot_height_in))

out_pdf <- file.path(fig2_dir, "fig2a_narrow_wrapped.pdf")
cairo_pdf(out_pdf, width = total_width_in, height = plot_height_in)
grid.newpage()
grid.draw(g)
dev.off()
message("[done] wrote ", out_pdf)

# ------------------------------- REPORTABLE VALUES ------------------------------
write.csv(
  data.frame(
    key = c("fig2a_narrow_panel_width_in", "fig2a_narrow_total_width_in", "fig2a_narrow_n_terms"),
    value = c(panel_width_in, round(total_width_in, 2), nrow(plot_df)),
    provenance = "figure_reproduction/scripts/35_fig2a_narrow_wrapped_plot.R"
  ),
  file.path(analysis_dir, "outputs", "reportable_values_35.csv"),
  row.names = FALSE
)
