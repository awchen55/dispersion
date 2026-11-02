library(variancePartition)
library(limma)
library(ggplot2)

source("scripts/00_config.R")

OUT_DIR <- "outputs"
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)

# Erik's sweep covers two more (finer) thresholds than the shared
# Brendan-vs-Erik comparison in script 02 -- defined locally rather than
# added to THRESHOLDS in 00_config.R, since THRESHOLDS is shared with scripts
# that require a matching Brendan fit at every threshold (script 02/03), which
# does not exist at 100/150.
ERIK_THRESHOLDS <- c(THRESHOLDS, 100, 150)

n_sig_genes <- function(fit, contrast) {
  tt <- topTable(fit, coef = contrast, number = Inf, adjust.method = "BH", sort.by = "none")
  sum(!is.na(tt$adj.P.Val) & tt$adj.P.Val < 0.05)
}

rows <- list()

for (thresh in ERIK_THRESHOLDS) {
  for (i in seq_len(nrow(SHARED_CONTRASTS))) {
    contrast   <- SHARED_CONTRASTS$contrast[i]
    erik_batch <- SHARED_CONTRASTS$erik_batch[i]
    erik_fit   <- readRDS(file.path(erik_dd_fit_dir(erik_batch), paste0("dd_fit_cell_count_", thresh, ".rds")))

    stopifnot(contrast %in% colnames(erik_fit$coefficients))

    rows[[length(rows) + 1]] <- data.frame(
      cell_count_threshold = thresh,
      contrast             = contrast,
      erik_batch           = toupper(erik_batch),
      n_sig_dd_genes       = n_sig_genes(erik_fit, contrast),
      stringsAsFactors = FALSE
    )
  }
}

df <- do.call(rbind, rows)
df$contrast_label <- gsub("_", " ", sub("^compare_", "", df$contrast))
df <- df[order(-df$cell_count_threshold, df$contrast), ]

out_csv <- file.path(OUT_DIR, "dd_gene_count_erik_by_threshold.csv")
write.csv(df, out_csv, row.names = FALSE)
message("Wrote: ", out_csv)
print(df)

df$cell_count_group <- factor(
  as.character(df$cell_count_threshold),
  levels = as.character(sort(ERIK_THRESHOLDS, decreasing = FALSE))
)

p <- ggplot(df, aes(x = cell_count_group, y = n_sig_dd_genes)) +
  geom_boxplot(outlier.shape = NA, fill = "#4477AA") +
  geom_jitter(width = 0.08, size = 1.5, alpha = 0.7) +
  labs(
    title = "Significant DD genes in Erik (B2+B6) data by cell count",
    subtitle = "each point = 1 of the 4 shared contrasts",
    x = "Cells per sample", y = "Significant DD genes (per contrast)"
  ) +
  theme_bw()

out_fig <- file.path(OUT_DIR, "figures", "dd_gene_count_erik_by_threshold.pdf")
ggsave(out_fig, p, width = 7, height = 5.5)
message("Wrote: ", out_fig)
