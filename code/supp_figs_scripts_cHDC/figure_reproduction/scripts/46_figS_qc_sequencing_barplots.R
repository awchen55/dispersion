#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
metrics_csv <- file.path(repo_dir, "cellranger_cluster-mode_trial", "sequencing_metrics_per_lane.csv")
out_dir <- file.path(repo_dir, "figure_reproduction", "outputs", "figS_qc_barplots")

# Published panels: metric column -> y limits (NULL = ggplot default).
panels <- list(
  Estimated.Number.of.Cells = c(0, 40000),
  Mean.Reads.per.Cell = NULL,
  Median.Genes.per.Cell = NULL,
  Reads.Mapped.to.Genome = c(0, 1),
  Reads.Mapped.Confidently.to.Genome = c(0, 1),
  Fraction.Reads.in.Cells = c(0, 1),
  Total.Genes.Detected = c(0, 40000)
)
n_lanes_expected <- 7

if (!file.exists(metrics_csv)) stop("Missing input: ", metrics_csv)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Same parsing as the source: comma-separated counts -> numeric (cols 19:20),
# percentages -> fractions (cols 6:18), after row.names = 1.
seq_metrics <- read.csv(metrics_csv, row.names = 1)
if (nrow(seq_metrics) != n_lanes_expected) stop("Expected ", n_lanes_expected, " lanes, found ", nrow(seq_metrics))
if (!identical(colnames(seq_metrics)[c(6, 18, 19, 20)],
               c("Valid.Barcodes", "Fraction.Reads.in.Cells", "Total.Genes.Detected", "Median.UMI.Counts.per.Cell"))) {
  stop("Unexpected column layout in ", metrics_csv, "; the source's positional parsing would be wrong")
}
seq_metrics[, 19:20] <- lapply(seq_metrics[, 19:20], function(x) as.numeric(gsub(",", "", x)))
seq_metrics[, 6:18] <- lapply(seq_metrics[, 6:18], function(x) as.numeric(gsub("%", "", x)) / 100)
# Count columns 2-4 are plotted too but were not converted in the source chunk.
count_cols <- c("Estimated.Number.of.Cells", "Mean.Reads.per.Cell", "Median.Genes.per.Cell")
seq_metrics[, count_cols] <- lapply(seq_metrics[, count_cols], function(x) as.numeric(gsub(",", "", x)))

for (metric in names(panels)) {
  if (!metric %in% colnames(seq_metrics)) stop("Missing metric column ", metric)
  vals <- seq_metrics[[metric]]
  if (anyNA(vals)) stop("NA after parsing ", metric)
  lims <- panels[[metric]]
  if (!is.null(lims) && (min(vals) < lims[1] || max(vals) > lims[2])) {
    stop(sprintf("%s values [%g, %g] fall outside the fixed limits [%g, %g]", metric, min(vals), max(vals), lims[1], lims[2]))
  }
  label <- gsub("\\.", " ", metric)
  plot <- ggplot(seq_metrics, aes(x = lane, y = .data[[metric]])) +
    geom_bar(stat = "identity") +
    labs(x = NULL, y = label, title = label) +
    (if (is.null(lims)) NULL else scale_y_continuous(limits = lims)) +
    theme_classic(base_family = "Arial") +
    theme(
      text = element_text(size = 6),
      axis.title = element_text(size = 6),
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
      axis.line = element_line(linewidth = 0.4),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
      plot.margin = margin(2, 2, 2, 2, unit = "mm"))
  out_pdf <- file.path(out_dir, paste0(gsub("\\.", "_", metric), "_barplot.pdf"))
  ggsave(out_pdf, plot, device = cairo_pdf, width = 2.5, height = 2, units = "in")
  message(sprintf("[46] %s: %d lanes, range %g-%g -> %s", metric, length(vals), min(vals), max(vals), basename(out_pdf)))
}
