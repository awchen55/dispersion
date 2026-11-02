#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(png)
  library(grid)
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
analysis_dir <- file.path(repo_dir, "figure_reproduction")
fig2_dir <- file.path(analysis_dir, "outputs", "fig2_replication")
raster_dpi <- 300
target_width_in <- 11
target_height_in <- 8.5

# Each row: (file, canonical manuscript letter) pairs, in display order.
row_specs <- list(
  list(files = c("fig2a_narrow_wrapped.pdf", "Fig2b_VarPart_ResVar_relabeled.pdf",
                 "fig2d_scatter_retained.pdf"),
       letters = c("A", "B", "D")),
  list(files = c("Fig2a_ResVar_CorrelationHeatmap_relabeled.pdf",
                 "fig2e_godag_top10_representatives.pdf"),
       letters = c("C", "E"))
)

# Rasterize one PDF at raster_dpi and return a list(grob=rasterGrob,
# width_in=, height_in=) using its NATIVE aspect ratio (pixel dims / dpi).
rasterize_pdf <- function(pdf_name) {
  pdf_path <- file.path(fig2_dir, pdf_name)
  if (!file.exists(pdf_path)) stop("Missing ", pdf_path, " -- run the script that produces it first")
  png_path <- tempfile(fileext = ".png")
  status <- system2("gs", c("-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=png16m",
                            paste0("-r", raster_dpi),
                            paste0("-sOutputFile=", shQuote(png_path)), shQuote(pdf_path)))
  if (status != 0) stop("gs (Ghostscript) rasterization failed on ", pdf_path)
  if (!file.exists(png_path)) stop("Expected rasterized output not found: ", png_path)
  img <- readPNG(png_path)
  unlink(png_path)
  list(
    grob = rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE),
    width_in = dim(img)[2] / raster_dpi,
    height_in = dim(img)[1] / raster_dpi
  )
}

# Tags applied directly per-panel, NOT via plot_annotation(tag_levels="A")
# -- that auto-tags every leaf in the patchwork tree, including the
# plot_spacer() padding elements added below.
tag_panel <- function(grob, letter) {
  wrap_elements(full = grob) +
    labs(tag = letter) +
    theme(plot.tag = element_text(size = 14, face = "bold"),
          plot.tag.position = "topleft")
}

# Build one row: rasterize its panels, scale to the row's own tallest
# native panel (only ever scales UP), tag each with its canonical letter.
build_row <- function(files, letters) {
  panels <- lapply(files, rasterize_pdf)
  row_height <- max(vapply(panels, `[[`, numeric(1), "height_in"))
  widths <- vapply(panels, function(p) p$width_in * (row_height / p$height_in), numeric(1))
  tagged <- Map(tag_panel, lapply(panels, `[[`, "grob"), letters)
  list(row = Reduce(`+`, tagged), height_in = row_height, widths_in = widths)
}

built_rows <- lapply(row_specs, function(spec) build_row(spec$files, spec$letters))
row_heights <- vapply(built_rows, `[[`, numeric(1), "height_in")
row_totals <- vapply(built_rows, function(r) sum(r$widths_in), numeric(1))
for (i in seq_along(built_rows)) {
  message(sprintf("[composite] row %d: height forced to %.2fin, widths %s",
                  i, row_heights[i], paste(round(built_rows[[i]]$widths_in, 2), collapse = ", ")))
}

# A row narrower than the widest row gets stretched by patchwork's `/`
# vertical stack to span the FULL outer canvas width -- padded with an
# invisible plot_spacer() sized to the exact leftover width instead, so real
# panels keep their correct aspect ratio and the gap becomes blank space.
canvas_width_in <- max(row_totals)
padded_rows <- Map(function(built, total) {
  spacer_width <- canvas_width_in - total
  if (spacer_width > 1e-6) {
    built$row + plot_spacer() + plot_layout(widths = c(built$widths_in, spacer_width))
  } else {
    built$row + plot_layout(widths = built$widths_in)
  }
}, built_rows, row_totals)

fig2_composite <- Reduce(`/`, padded_rows) +
  plot_layout(heights = row_heights)

canvas_height_in <- sum(row_heights)

# Uniform scale-down to fit within 11x8.5in -- min() of the two ratios so
# neither dimension overflows; no distortion.
fit_scale <- min(target_width_in / canvas_width_in, target_height_in / canvas_height_in)
out_width_in <- canvas_width_in * fit_scale
out_height_in <- canvas_height_in * fit_scale

out_pdf <- file.path(fig2_dir, "Figure2_composite.pdf")
ggsave(out_pdf, fig2_composite, device = cairo_pdf,
       width = out_width_in, height = out_height_in, units = "in", limitsize = FALSE)
message(sprintf("[done] wrote %s (%.2fx%.2fin, scaled %.3fx to fit %sx%sin)",
                out_pdf, out_width_in, out_height_in, fit_scale, target_width_in, target_height_in))
