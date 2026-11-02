#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
fig_repro_dir <- file.path(repo_dir, "figure_reproduction")
out_dir <- file.path(fig_repro_dir, "outputs", "dd_trend_correction")
effects_csv <- file.path(out_dir, "dd_trend_corrected_effects.csv")

if (!file.exists(effects_csv)) {
  stop("Missing ", effects_csv, " -- run 14_dd_trend_corrected_effects.R first")
}
effects_df <- read.csv(effects_csv, stringsAsFactors = FALSE)
stopifnot(
  "input table must contain the expected columns" =
    all(c("contrast", "gene", "dd_call_current", "observed_logFC",
          "expected_artifact_gaussian") %in% colnames(effects_df)),
  "input table must contain no non-finite observed_logFC / expected_artifact_gaussian among DD calls" =
    all(is.finite(effects_df$observed_logFC[effects_df$dd_call_current])) &&
    all(is.finite(effects_df$expected_artifact_gaussian[effects_df$dd_call_current])))

dd_df <- effects_df[effects_df$dd_call_current, ]
message(sprintf("[load] %d current DD calls (of %d total gene-contrast pairs) across %d contrasts",
                nrow(dd_df), nrow(effects_df), length(unique(dd_df$contrast))))

# --------------------------------- MAIN --------------------------------------

long_df <- rbind(
  data.frame(contrast = dd_df$contrast, quantity = "Observed |logFC|",
             value = abs(dd_df$observed_logFC), stringsAsFactors = FALSE),
  data.frame(contrast = dd_df$contrast, quantity = "Expected |artifact| (gaussian)",
             value = abs(dd_df$expected_artifact_gaussian), stringsAsFactors = FALSE))
long_df$quantity <- factor(long_df$quantity,
                           levels = c("Observed |logFC|", "Expected |artifact| (gaussian)"))

p <- ggplot(long_df, aes(x = value, fill = quantity)) +
  geom_histogram(position = "identity", alpha = 0.55, bins = 40, linewidth = 0) +
  scale_fill_manual(name = NULL, values = c("Observed |logFC|" = "#4575b4",
                                            "Expected |artifact| (gaussian)" = "#d73027")) +
  facet_wrap(~ contrast, ncol = 5, scales = "free_y") +
  labs(x = "Absolute value", y = "Count (current DD calls)",
       title = "Observed DD effect size vs. mean-expression-trend-predicted artifact, per contrast",
       subtitle = "Only genes currently called DD (adj.P.Val<=0.05)") +
  theme_classic(base_family = "Arial") +
  theme(
    text = element_text(size = 10), axis.text = element_text(size = 10),
    axis.title = element_text(size = 12), plot.title = element_text(size = 14),
    plot.subtitle = element_text(size = 10), strip.text = element_text(size = 6),
    legend.text = element_text(size = 10), legend.title = element_text(size = 10),
    legend.key.size = unit(10, "pt"), legend.position = "top",
    axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
    plot.margin = margin(2, 2, 2, 2, unit = "mm"))

out_file <- file.path(out_dir, "dd_effect_vs_artifact_histograms.pdf")
ggsave(out_file, p, width = 14, height = 9, units = "in", device = cairo_pdf, limitsize = FALSE)
message("[done] Saved: ", out_file)
