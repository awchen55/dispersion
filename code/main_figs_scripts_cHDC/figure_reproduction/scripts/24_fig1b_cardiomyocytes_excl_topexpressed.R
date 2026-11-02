#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(viridisLite)
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
analysis_dir <- file.path(repo_dir, "figure_reproduction")
filt_mat_dir <- file.path(analysis_dir, "outputs", "filtered_matrices")
out_dir <- file.path(analysis_dir, "outputs", "fig1_replication")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cell_type <- "Cardiomyocytes"
sample <- "NA19099_Cardiomyocytes"  # source Rmd's representative Cardiomyocytes sample
published_rho <- 0.0115  # Fig1b_card_dispvsmean_spearman.pdf, full gene universe

# ------------------------------- LOAD -----------------------------------------
mean_path <- file.path(filt_mat_dir, paste0(cell_type, "_mean_retained.csv"))
disp_path <- file.path(filt_mat_dir, paste0(cell_type, "_dispersion_retained.csv"))
if (!file.exists(mean_path) || !file.exists(disp_path)) {
  stop("Missing filtered matrix for ", cell_type,
       " -- run 01_exclude_topexpressed_genes.R first: ", mean_path, " / ", disp_path)
}

mean_mat <- read.csv(mean_path, row.names = 1, check.names = FALSE)
disp_mat <- read.csv(disp_path, row.names = 1, check.names = FALSE)
if (!sample %in% colnames(mean_mat) || !sample %in% colnames(disp_mat)) {
  stop("Expected column ", sample, " not found in retained matrices")
}
if (!identical(rownames(mean_mat), rownames(disp_mat))) {
  stop("Retained mean/dispersion matrices have non-identical row names/order")
}

df <- data.frame(
  gene = rownames(mean_mat),
  mean = mean_mat[[sample]],
  res_var = disp_mat[[sample]]
)

n_retained <- nrow(df)
message(sprintf("[%s] %d retained genes (-6-excluded, ribo-dropped)", cell_type, n_retained))

# ------------------------------- STATS -----------------------------------------
cor_test <- cor.test(df$mean, df$res_var, method = "spearman")
rho_val <- round(cor_test$estimate, 4)
message(sprintf("[%s] retained-set Spearman rho = %.4f (published, full universe: %.4f)",
                cell_type, rho_val, published_rho))

# ------------------------------- PLOT (byte-for-byte copy of source Rmd spec) --
plot <- ggplot(df, aes(x = mean, y = res_var)) +
  geom_point(
    color = viridis(1, option = "C"),
    alpha = 0.2,
    size = 0.4
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("Spearman ρ = ", rho_val),
    hjust = 1.1,
    vjust = 1.5,
    size = 2  # ~6 pt at final 90 mm size
  ) +
  theme_classic(base_family = "Arial") +
  labs(title = cell_type) +
  theme(
    axis.text.x = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(size = 7, angle = 90),
    title = element_text(size = 7)
  ) +
  labs(y = "Dispersion", x = "Mean")

out_pdf <- file.path(out_dir, "Fig1b_card_dispvsmean_spearman_excl_topexpressed.pdf")
ggsave(
  filename = out_pdf,
  plot = plot,
  device = cairo_pdf,   # embeds fonts properly, matches source
  width = 3,            # 90 mm, matches source ggsave call
  height = 2,
  units = "in"
)

# ------------------------------- REPORTABLE VALUES ------------------------------
write.csv(
  data.frame(
    key = c("fig1b_cardiomyocytes_retained_spearman_rho",
            "fig1b_cardiomyocytes_retained_n_genes"),
    value = c(rho_val, n_retained),
    provenance = "figure_reproduction/scripts/24_fig1b_cardiomyocytes_excl_topexpressed.R"
  ),
  file.path(analysis_dir, "outputs", "reportable_values_24.csv"),
  row.names = FALSE
)

message("[done] wrote ", out_pdf)
