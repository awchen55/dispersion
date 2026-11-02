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
n_top <- 10
# Label wrap width: strip the "GOBP_" prefix and underscores first, then
# wrap. width=18 is the narrowest wrap width that still gives the longest
# term exactly 3 lines (verified in the source script against this exact
# term set).
label_wrap_width <- 18
# Fixed panel (data-area) width via gtable surgery (see below).
panel_width_in <- 1.5

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
reps_csv <- file.path(fig2_dir, "fig2e_godag_semantic_representatives.csv")
if (!file.exists(reps_csv)) stop("Missing ", reps_csv, " -- run 32_fig2e_go_dag_semantic_clustering.R first")
reps <- read.csv(reps_csv)

top10 <- reps %>%
  dplyr::arrange(p.adjust) %>%
  head(n_top)
message(sprintf("[fig2e:top10] %d/%d representatives are heart/cardiac-related",
                sum(top10$is_heart_related), nrow(top10)))

plot_df <- top10 %>%
  mutate(GeneRatio_num = as.numeric(str_split_fixed(GeneRatio, "/", 2)[, 1]) /
                          as.numeric(str_split_fixed(GeneRatio, "/", 2)[, 2]),
         log10_padj = -log10(p.adjust),
         Description = str_wrap(gsub("_", " ", sub("^GOBP_", "", Description)),
                                width = label_wrap_width)) %>%
  arrange(desc(GeneRatio_num)) %>%
  mutate(Description = factor(Description, levels = rev(Description)))

# Title removed (redundant with the manuscript's own figure caption).
p <- make_go_dotplot(plot_df, title = NULL)

# 10% x-scale expansion (doubled from the default 5%): at panel_width_in=1.5
# the default left too little absolute margin for the largest/rightmost
# points' fixed real-unit radius, clipping them against the panel's right edge.
p <- p + scale_x_continuous(expand = expansion(mult = 0.10))

# ------------------------------- GTABLE SURGERY: FIX PANEL WIDTH -----------------
# Convert to a gtable, locate the "panel" cell, force its width to
# panel_width_in -- every other column (y-axis label, legend) keeps whatever
# width it already auto-sized to.
#
# Measured inside a throwaway cairo_pdf device first (ggplotGrob()/grid width
# units are lazily resolved against whichever device is active at measurement
# time, and cairo's Arial metrics differ from R's default null device) then
# drawn on the real output device at the now-correctly-measured size.
measure_pdf <- tempfile(fileext = ".pdf")
cairo_pdf(measure_pdf, width = 9, height = 4)
g <- ggplotGrob(p)
panel_idx <- which(g$layout$name == "panel")
stopifnot("Expected exactly one panel grob in the gtable layout" = length(panel_idx) == 1)
g$widths[g$layout$l[panel_idx]] <- unit(panel_width_in, "in")
total_width_in <- sum(convertWidth(g$widths, "in", valueOnly = TRUE))
dev.off()
unlink(measure_pdf)
message(sprintf("[panel-width] panel fixed to %.2fin; total canvas width = %.2fin (was 9in)",
                panel_width_in, total_width_in))

out_pdf <- file.path(fig2_dir, "fig2e_godag_top10_representatives.pdf")
cairo_pdf(out_pdf, width = total_width_in, height = 4)
grid.newpage()
grid.draw(g)
dev.off()
message("[done] wrote ", out_pdf)

# ------------------------------- REPORTABLE VALUES ------------------------------
write.csv(
  data.frame(
    key = c("fig2e_top10_n_heart_related"),
    value = c(sum(top10$is_heart_related)),
    provenance = "figure_reproduction/scripts/34_fig2e_godag_top10_representatives_plot.R"
  ),
  file.path(analysis_dir, "outputs", "reportable_values_34.csv"),
  row.names = FALSE
)
