#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
analysis_dir <- file.path(repo_dir, "figure_reproduction")
fig6_dir <- file.path(analysis_dir, "outputs", "fig6_dqtl")
dir.create(fig6_dir, recursive = TRUE, showWarnings = FALSE)

RUN_ID <- "20260925_143116"
RUN_DIR <- file.path(repo_dir, "analysis/dqtl_mapping/outputs/dqtl_results", RUN_ID)
if (!dir.exists(RUN_DIR)) stop("Missing ", RUN_DIR)
CELL_TYPES <- c("Activated_fibroblasts", "Cardiac_fibroblasts", "Cardiomyocytes",
                "Foregut_endodermal_cells")

# Same values as 02_qq_plots.R.
THIN_P_CUTOFF <- 0.001     # keep every test below this; subsample the rest
THIN_N_BULK <- 3000        # random tests kept per cell type above the cutoff
ENVELOPE_ALPHA <- 0.05     # 95% pointwise envelope
N_ENVELOPE_POINTS <- 500

# ------------------------------- QQ PANELS -----------------------------------
set.seed(42)
qq_panel <- function(ct) {
  f <- file.path(RUN_DIR, ct, paste0("dQTL_", ct, "_AllChroms.txt"))
  if (!file.exists(f)) stop("Missing ", f)
  dt <- fread(f, sep = "\t", select = "pvalue")
  n <- nrow(dt)
  if (n == 0) stop("No tests in ", f)
  if (anyNA(dt$pvalue) || any(dt$pvalue < 0 | dt$pvalue > 1)) stop("Invalid p-values in ", f)
  message(sprintf("[QQ] %s: %d cis tests", ct, n))

  setorder(dt, pvalue)
  dt[, expected := -log10(seq_len(n) / (n + 1))]
  dt[, observed := -log10(pmax(pvalue, 1e-300))]
  bulk <- dt[pvalue >= THIN_P_CUTOFF]
  plot_dt <- rbind(dt[pvalue < THIN_P_CUTOFF], bulk[sample(.N, min(.N, THIN_N_BULK))])

  k <- unique(round(10 ^ seq(0, log10(n), length.out = N_ENVELOPE_POINTS)))
  env <- data.table(
    expected = -log10(k / (n + 1)),
    lower = -log10(qbeta(1 - ENVELOPE_ALPHA / 2, k, n - k + 1)),
    upper = -log10(qbeta(ENVELOPE_ALPHA / 2, k, n - k + 1))
  )

  p <- ggplot(plot_dt, aes(x = expected, y = observed)) +
    geom_ribbon(data = env, aes(x = expected, ymin = lower, ymax = upper),
                inherit.aes = FALSE, fill = "grey70", alpha = 0.5) +
    geom_point(size = 0.3, alpha = 0.5, color = "steelblue") +
    geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.5) +
    labs(title = gsub("_", " ", ct), x = "Expected -log10(p)", y = "Observed -log10(p)") +
    theme_classic(base_family = "Arial") +
    theme(text = element_text(size = 6), plot.title = element_text(hjust = 0.5))
  list(plot = p, n = n)
}
qq <- lapply(CELL_TYPES, qq_panel)
names(qq) <- CELL_TYPES

# ------------------------------- DOSAGE PANEL (HADHB) ------------------------
# Per-donor data CSV written by analysis/dqtl_mapping_q02/scripts/04_dosage_boxplots.R,
# plotted with that script's own panel helper.
setwd(repo_dir)  # utils.R's helpers use repo-root-relative paths
source(file.path(repo_dir, "analysis/dqtl_genotype_dosage_plots/scripts/utils.R"))

dosage_gene <- "HADHB"
dosage_ct <- "Cardiomyocytes"
dosage_snp <- "chr2:26200545:G_A"  # geno_snp_id from the standalone PDF's title
dosage_csv <- file.path(repo_dir, "analysis/dqtl_mapping_q02/outputs/genotype_dosage_plots/data",
                        paste0(dosage_gene, "_", dosage_ct, "_dosage_data_", RUN_ID, ".csv"))
if (!file.exists(dosage_csv)) stop("Missing ", dosage_csv)
dosage_df <- read.csv(dosage_csv, colClasses = c(dosage = "character"))
dosage_df$dosage <- factor(dosage_df$dosage, levels = c("0", "1", "2"))
if (anyNA(dosage_df$dosage)) stop("Unexpected dosage value(s) in ", dosage_csv)

snp_parts <- strsplit(dosage_snp, "[:_]")[[1]]  # chr, pos, ref, alt
dosage_ref <- snp_parts[3]
dosage_alt <- snp_parts[4]
# Major/minor from this cohort's own allele frequency (dosage = ALT copies).
alt_freq <- mean(as.numeric(as.character(dosage_df$dosage))) / 2
major_allele <- if (alt_freq > 0.5) dosage_alt else dosage_ref
minor_allele <- if (alt_freq > 0.5) dosage_ref else dosage_alt
message(sprintf("[dosage] %s %s: n=%d donors, ALT (%s) freq=%.3f -> major %s, minor %s",
                dosage_gene, dosage_snp, nrow(dosage_df), dosage_alt, alt_freq, major_allele, minor_allele))

p_dos_mean <- dosage_boxplot_panel(dosage_df, "mean_memento", "Mean expression", dosage_gene,
                                   dosage_ref, dosage_alt)
p_dos_disp <- dosage_boxplot_panel(dosage_df, "dispersion", "Dispersion", NULL,
                                   dosage_ref, dosage_alt)

# ------------------------------- FONT SIZES ----------------------------------
# Same composite font sizes as the Fig3 composite (script 28).
axis_text_bump <- theme(axis.text = element_text(size = 7))
axis_title_bump <- theme(axis.title = element_text(size = 8))
subplot_title_bump <- theme(plot.title = element_text(size = 8))

qq_plots <- lapply(qq, function(x) x$plot + axis_text_bump + axis_title_bump + subplot_title_bump)
p_dos_mean <- p_dos_mean + axis_text_bump + axis_title_bump + subplot_title_bump
p_dos_disp <- p_dos_disp + axis_text_bump + axis_title_bump

# ------------------------------- ASSEMBLY ------------------------------------
# Each side is wrapped with wrap_elements() so the outer auto-tagging gives it
# one letter (A = QQ grid, B = dosage) instead of tagging every leaf plot.
qq_grid <- wrap_plots(qq_plots, ncol = 2)
qq_unit <- wrap_elements(full = qq_grid)

dos_title <- paste0(gsub("_", " ", dosage_ct), "\n", snp_parts[1], ":", snp_parts[2], ", ",
                    major_allele, " (major) / ", minor_allele, " (minor)")
p_dos_combined <- (p_dos_mean / p_dos_disp) +
  plot_annotation(title = dos_title,
                  theme = theme(plot.title = element_text(size = 8, family = "Arial", hjust = 0.5)))
dos_unit <- wrap_elements(full = p_dos_combined)

fig6_composite <- qq_unit + dos_unit +
  plot_layout(widths = c(2, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

out_pdf <- file.path(fig6_dir, "Figure6_dqtl_composite.pdf")
ggsave(out_pdf, fig6_composite, device = cairo_pdf, width = 7.2, height = 5.4, units = "in")
message("[done] wrote ", out_pdf)

# Standalone dosage panel (reviewer response; previously written by script 28).
out_dosage_pdf <- file.path(fig6_dir, "HADHB_Cardiomyocytes_dosage_standalone.pdf")
ggsave(out_dosage_pdf, p_dos_combined, device = cairo_pdf, width = 2.4, height = 5.0, units = "in")
ggsave(sub("\\.pdf$", ".png", out_dosage_pdf), p_dos_combined, device = ragg::agg_png,
       width = 2.4, height = 5.0, units = "in", dpi = 600, bg = "white")
message("[done] wrote ", out_dosage_pdf)

# ------------------------------- REPORTABLE VALUES ---------------------------
reportable <- data.frame(
  key = c(paste0("fig6_n_cis_tests_", CELL_TYPES), "fig6_hadhb_n_donors", "fig6_hadhb_alt_allele_freq"),
  value = c(vapply(qq, function(x) x$n, numeric(1)), nrow(dosage_df), alt_freq),
  provenance = "figure_reproduction/scripts/42_fig6_dqtl_composite.R"
)
write.csv(reportable, file.path(analysis_dir, "outputs", "reportable_values_42.csv"), row.names = FALSE)
