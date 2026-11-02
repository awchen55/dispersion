#!/usr/bin/env Rscript

stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(readxl)
})

# ------------------------------- CONSTANTS -----------------------------------
FR <- "figure_reproduction"
OUT_DIR <- file.path(FR, "outputs", "subcluster_retained_genesets")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
LOW_TXT  <- file.path(FR, "outputs", "ranked", "consistently_low_retained.txt")
HIGH_TXT <- file.path(FR, "outputs", "ranked", "consistently_high_retained.txt")
RETAINED_MATRIX_CSV <- file.path(FR, "outputs", "retained_gene_matrix.csv")

SUBCLUSTER_DIR <- "analysis/subcluster_dispersion_robustness/outputs/mincells1000_minind2"
RESOLUTIONS <- c("0.5", "1.0", "2.0")
# Sub-clusters composed of the excluded Epithelial sub-cluster-24 cells
# (cross-tabulated from subcluster_raw_counts.h5ad obs, 2026-10-08).
DROP_SUBCLUSTERS <- c("0.5" = "Epithelial_cells::12", "1.0" = "Epithelial_cells::13",
                      "2.0" = "Epithelial_cells::24")

CELL_TYPES <- c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells")
RIBO_PATTERN <- "^RP"
MITO_PATTERN <- "^MT-"
MEAN_EXCL_THRESHOLD <- -6
N_QUINTILES <- 5L

# External feature tables (identical to functional_enrichment_percelltype/04-05).
TSS_TXT <- "/project2/gilad/wenhe/cox_baboon_diet/combined_analysis/gene_annotation/fantom/short_descrip_split_gene_only.txt"
ABC_ENHANCER_TXT <- "/project2/gilad/wenhe/cox_baboon_diet/combined_analysis/gene_annotation/ABC_prediction_score005.txt"
PPI_TSV <- "cellranger_cluster-mode_trial/analysis/datasets/external_data/protein_interaction_data/interaction_consensus.tsv.zip"
PPI_TISSUE_TSV <- "cellranger_cluster-mode_trial/analysis/datasets/external_data/protein_interaction_data/has_data_Tissue.tsv"
TWN_S4_XLSX <- "cellranger_cluster-mode_trial/analysis/datasets/external_data/saha_coexp_net_TWN_2017/Supplemental_Data_S4.xlsx"
MANE_SUMMARY_TXT <- "data/raw/mane_v1_5/MANE.GRCh38.v1.5.summary.txt.gz"
MANE_GFF_GZ <- "data/raw/mane_v1_5/MANE.GRCh38.v1.5.ensembl_genomic.gff.gz"
EPD_BED <- "cellranger_cluster-mode_trial/analysis/datasets/external_data/epd/human_epdnew_3sa0y.bed"

for (f in c(LOW_TXT, HIGH_TXT, RETAINED_MATRIX_CSV, TSS_TXT, ABC_ENHANCER_TXT, PPI_TSV,
            PPI_TISSUE_TSV, TWN_S4_XLSX, MANE_SUMMARY_TXT, MANE_GFF_GZ, EPD_BED,
            file.path(SUBCLUSTER_DIR, sprintf("res_var_by_subcluster_res_%s.csv", RESOLUTIONS)),
            file.path(SUBCLUSTER_DIR, sprintf("mean_expression_by_subcluster_res_%s.csv", RESOLUTIONS)))) {
  if (!file.exists(f)) stop("Missing required input: ", f)
}

log_msg <- function(...) message("[44_subcluster_retained_genesets] ", ...)

# ------------------------------- GENE SETS -----------------------------------
base_low  <- readLines(LOW_TXT)
base_high <- readLines(HIGH_TXT)
retained_ind <- as.matrix(read.csv(RETAINED_MATRIX_CSV, row.names = 1, check.names = FALSE))
union_universe <- rownames(retained_ind)[rowSums(retained_ind[, CELL_TYPES] == 1) == length(CELL_TYPES)]
stopifnot(all(base_low %in% union_universe), all(base_high %in% union_universe))
log_msg("Base sets: ", length(base_low), " low, ", length(base_high), " high; universe ",
        length(union_universe), " genes")

# Same loader as subcluster 03's load_memento_csv(): genes x groups, complete cases only.
load_memento_csv <- function(path) {
  df <- read.csv(path, row.names = 1, check.names = FALSE)
  if ("gene" %in% colnames(df)) { rownames(df) <- df$gene; df$gene <- NULL }
  as.matrix(df[stats::complete.cases(df), , drop = FALSE])
}

quintile_rank <- function(v) {
  r <- rank(v, ties.method = "average", na.last = "keep")
  as.integer(cut(r, breaks = quantile(r, probs = seq(0, 1, length.out = N_QUINTILES + 1), na.rm = TRUE),
                 include.lowest = TRUE, labels = FALSE))
}

retained_sets <- list()
size_rows <- list()
for (res in RESOLUTIONS) {
  mat <- load_memento_csv(file.path(SUBCLUSTER_DIR, sprintf("res_var_by_subcluster_res_%s.csv", res)))
  mean_mat <- load_memento_csv(file.path(SUBCLUSTER_DIR, sprintf("mean_expression_by_subcluster_res_%s.csv", res)))
  common <- intersect(rownames(mat), rownames(mean_mat))
  common <- common[!(grepl(RIBO_PATTERN, common) | grepl(MITO_PATTERN, common))]
  mat <- mat[common, , drop = FALSE]
  mean_mat <- mean_mat[common, , drop = FALSE]

  cols <- colnames(mat)
  ind <- regmatches(cols, regexpr("NA[0-9]+", cols))
  sc <- regmatches(cols, regexpr("[A-Za-z_]+::[0-9]+", cols))
  if (length(ind) != length(cols) || length(sc) != length(cols)) stop("Unparseable group columns at res ", res)
  if (!all(sub("::.*$", "", sc) %in% CELL_TYPES)) stop("Unexpected parent cell type at res ", res)
  if (!DROP_SUBCLUSTERS[[res]] %in% sc) stop(DROP_SUBCLUSTERS[[res]], " not found at res ", res)

  # -6 any-individual mask per sub-cluster, then mean across individuals.
  subclusters <- setdiff(unique(sc), DROP_SUBCLUSTERS[[res]])
  sub_agg <- sapply(subclusters, function(s) {
    s_cols <- cols[sc == s]
    v <- rowMeans(mat[, s_cols, drop = FALSE])
    masked <- apply(mean_mat[, s_cols, drop = FALSE], 1, function(x) any(x > MEAN_EXCL_THRESHOLD))
    v[masked] <- NA_real_
    v
  })
  sub_agg <- sub_agg[intersect(rownames(sub_agg), union_universe), , drop = FALSE]
  quint <- apply(sub_agg, 2, quintile_rank)
  rownames(quint) <- rownames(sub_agg)

  all_q <- function(genes, q) {
    g <- intersect(genes, rownames(quint))
    keep <- apply(quint[g, , drop = FALSE], 1, function(x) all(!is.na(x) & x == q))
    g[keep]
  }
  low_r  <- all_q(base_low, 1L)
  high_r <- all_q(base_high, N_QUINTILES)
  retained_sets[[res]] <- list(low = low_r, high = high_r)
  writeLines(low_r, file.path(OUT_DIR, sprintf("retained_low_res_%s.txt", res)))
  writeLines(high_r, file.path(OUT_DIR, sprintf("retained_high_res_%s.txt", res)))
  size_rows[[res]] <- data.frame(
    resolution = res, n_subclusters = length(subclusters), n_genes_ranked = nrow(quint),
    n_low_base = length(base_low), n_low_measured = length(intersect(base_low, rownames(quint))),
    n_low_retained = length(low_r),
    n_high_base = length(base_high), n_high_measured = length(intersect(base_high, rownames(quint))),
    n_high_retained = length(high_r))
  log_msg(sprintf("res %s: %d sub-clusters (dropped %s), %d genes ranked; retained %d/%d low, %d/%d high",
                  res, length(subclusters), DROP_SUBCLUSTERS[[res]], nrow(quint),
                  length(low_r), length(base_low), length(high_r), length(base_high)))
  if (length(low_r) == 0 || length(high_r) == 0) stop("Empty retained set at res ", res)
}
size_df <- do.call(rbind, size_rows); rownames(size_df) <- NULL
write.csv(size_df, file.path(OUT_DIR, "retained_set_sizes.csv"), row.names = FALSE)
print(size_df)

# ------------------------------- FEATURE LOOKUPS -----------------------------
# Copied from functional_enrichment_percelltype/05 (same files, same parsing).
tss_usage <- read.csv(TSS_TXT, header = FALSE, stringsAsFactors = FALSE)
tss_table <- as.data.frame(table(tss_usage$V1), stringsAsFactors = FALSE)
tss_lookup <- setNames(tss_table$Freq, tss_table$Var1)

abc_pred <- readr::read_table(ABC_ENHANCER_TXT, show_col_types = FALSE)
stopifnot(all(c("TargetGene", "avg_length_across_celltype") %in% colnames(abc_pred)),
          !anyDuplicated(abc_pred$TargetGene))
enh_lookup <- setNames(abc_pred$avg_length_across_celltype, abc_pred$TargetGene)

ppi <- readr::read_tsv(PPI_TSV, show_col_types = FALSE)
exp_evidence <- readr::read_tsv(PPI_TISSUE_TSV, show_col_types = FALSE)
ppi$gene_1 <- exp_evidence$Gene[match(ppi$ensembl_gene_id_1, exp_evidence$Ensembl)]
ppi_table <- as.data.frame(table(ppi$gene_1), stringsAsFactors = FALSE)
ppi_lookup <- setNames(ppi_table$Freq, ppi_table$Var1)

s4_sheets <- readxl::excel_sheets(TWN_S4_XLSX)
stopifnot("Supplemental_Data_S4.xlsx sheet 10 layout changed" = grepl("^Group3", s4_sheets[10]))
shared_te_te <- readxl::read_excel(TWN_S4_XLSX, sheet = s4_sheets[10], col_names = TRUE)
twn_lookup <- setNames(as.integer(shared_te_te$rank), shared_te_te$hub)

strip_version <- function(x) sub("\\..*$", "", x)
mane_summary <- read.delim(MANE_SUMMARY_TXT, stringsAsFactors = FALSE, comment.char = "")
colnames(mane_summary)[1] <- sub("^#", "", colnames(mane_summary)[1])
mane_select <- mane_summary[mane_summary$MANE_status == "MANE Select", , drop = FALSE]
mane_select$ensembl_transcript_id <- strip_version(mane_select$Ensembl_nuc)
mane_gff <- read.delim(MANE_GFF_GZ, header = FALSE, comment.char = "#", stringsAsFactors = FALSE,
                       col.names = c("seqid", "source", "feature", "start", "end",
                                     "score", "strand", "frame", "attributes"))
exon_gff <- mane_gff[mane_gff$feature == "exon", , drop = FALSE]
exon_gff$ensembl_transcript_id <- strip_version(sub("^transcript_id=", "",
  regmatches(exon_gff$attributes, regexpr("transcript_id=[^;]+", exon_gff$attributes))))
exon_counts <- table(exon_gff$ensembl_transcript_id)
mane_select$num_introns <- as.integer(exon_counts[mane_select$ensembl_transcript_id]) - 1L
stopifnot("some MANE Select transcripts have no exon rows" = !any(is.na(mane_select$num_introns)))
intron_lookup <- setNames(mane_select$num_introns, mane_select$symbol)

FEATURE_LOOKUPS <- list(tss_count = tss_lookup, enhancer_length = enh_lookup,
                        ppi_connectivity = ppi_lookup, twn_hub_rank = twn_lookup,
                        intron_count = intron_lookup)
FEATURE_LABELS <- c(tss_count = "TSS count", enhancer_length = "Enhancer length (bp)",
                    ppi_connectivity = "PPI connections", twn_hub_rank = "TWN shared-hub rank",
                    intron_count = "Intron count")

epd_tata <- read.table(EPD_BED, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
tata_genes <- intersect(unique(sub("_[0-9]+$", "", epd_tata$V4)), union_universe)

# ------------------------------- TESTS ---------------------------------------
hl_rows <- list(); tata_rows <- list()
for (res in RESOLUTIONS) {
  low_r <- retained_sets[[res]]$low; high_r <- retained_sets[[res]]$high
  for (feat in names(FEATURE_LOOKUPS)) {
    x <- unname(FEATURE_LOOKUPS[[feat]][low_r]);  x <- x[!is.na(x)]
    y <- unname(FEATURE_LOOKUPS[[feat]][high_r]); y <- y[!is.na(y)]
    if (length(x) == 0 || length(y) == 0) stop("Empty low/high values for ", feat, " at res ", res)
    wt <- wilcox.test(x, y, conf.int = TRUE)
    hl_rows[[paste(res, feat)]] <- data.frame(
      resolution = res, feature = feat, n_low = length(x), n_high = length(y),
      median_low = median(x), median_high = median(y), hl_shift = unname(wt$estimate),
      ci_lower = wt$conf.int[1], ci_upper = wt$conf.int[2], p_value = wt$p.value)
  }
  for (set_name in c("low", "high")) {
    target <- if (set_name == "low") low_r else high_r
    a <- length(intersect(tata_genes, target)); b <- length(setdiff(tata_genes, target))
    c_ <- length(setdiff(target, tata_genes))
    d <- length(setdiff(union_universe, union(tata_genes, target)))
    stopifnot(a + b + c_ + d == length(union_universe))
    ft <- fisher.test(matrix(c(a, b, c_, d), nrow = 2, byrow = TRUE))
    tata_rows[[paste(res, set_name)]] <- data.frame(
      resolution = res, gene_set = set_name, n_target = length(target),
      n_tata_in_target = a, odds_ratio = unname(ft$estimate),
      ci_lower = ft$conf.int[1], ci_upper = ft$conf.int[2], p_value = ft$p.value)
  }
}
hl_df <- do.call(rbind, hl_rows); rownames(hl_df) <- NULL
tata_df <- do.call(rbind, tata_rows); rownames(tata_df) <- NULL
stopifnot(nrow(hl_df) == length(RESOLUTIONS) * length(FEATURE_LOOKUPS))
write.csv(hl_df, file.path(OUT_DIR, "hl_shift_by_resolution.csv"), row.names = FALSE)
write.csv(tata_df, file.path(OUT_DIR, "tata_by_resolution.csv"), row.names = FALSE)
print(hl_df); print(tata_df)

# ------------------------------- PLOT (05b style) ----------------------------
FONT_SCALE <- 1.5
res_levels <- c("0.5" = "0.5", "1.0" = "1", "2.0" = "2")
hl_df$res_label <- factor(res_levels[hl_df$resolution], levels = res_levels)
tata_df$res_label <- factor(res_levels[tata_df$resolution], levels = res_levels)
tata_df$log2_or <- log2(tata_df$odds_ratio)
tata_df$log2_ci_lower <- log2(tata_df$ci_lower)
tata_df$log2_ci_upper <- log2(tata_df$ci_upper)
tata_df$gene_set <- factor(tata_df$gene_set, levels = c("low", "high"))

# Plotmath p-value labels (parsed by geom_text): "P < 10^-15" below 1e-15
# (bucket used in the other enrichment figures), otherwise the exact value as
# "P = 5.01 x 10^-11" (main-figure style).
format_pvalue_label <- function(p) {
  exponent <- as.integer(ifelse(p > 0, floor(log10(p)), 0))  # p = 0 falls in the < 1e-15 bucket
  mantissa <- sprintf("%.2f", p / 10^exponent)
  ifelse(p < 1e-15, "P < 10^{-15}",
         sprintf('P == "%s" %%*%% 10^{%d}', mantissa, exponent))
}
SET_COLORS <- c(low = "steelblue", high = "coral")

# One point per resolution (per gene set if `by_set`); p-value label placed
# just beyond the CI end farther from 0.
make_panel <- function(df, y_col, ymin_col, ymax_col, panel_title, y_lab = "", show_x_text = TRUE,
                       by_set = FALSE) {
  away_up <- df[[y_col]] >= 0
  df$label_y <- ifelse(away_up, df[[ymax_col]], df[[ymin_col]])
  df$label_vjust <- ifelse(away_up, -0.6, 1.6)
  df$p_label <- format_pvalue_label(df$p_value)
  df$grp <- if (by_set) df$gene_set else factor("all")
  dodge <- position_dodge(width = if (by_set) 0.45 else 0)
  p <- ggplot(df, aes(x = res_label, y = .data[[y_col]], group = grp)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey40")
  if (by_set) {
    p <- p +
      geom_errorbar(aes(ymin = .data[[ymin_col]], ymax = .data[[ymax_col]], color = gene_set),
                    width = 0.08, linewidth = 0.3, alpha = 0.6, position = dodge) +
      geom_line(aes(color = gene_set), linewidth = 0.5, position = dodge) +
      geom_point(aes(color = gene_set), size = 1.4, position = dodge) +
      scale_color_manual(values = SET_COLORS, labels = c(low = "Low dispersion", high = "High dispersion"),
                         name = NULL) +
      theme_classic(base_family = "Arial") +
      theme(legend.position = "inside", legend.position.inside = c(0.02, 0.5),
            legend.justification = c(0, 0.5), legend.text = element_text(size = 6),
            legend.key.size = unit(7, "pt"), legend.background = element_blank())
  } else {
    p <- p +
      geom_errorbar(aes(ymin = .data[[ymin_col]], ymax = .data[[ymax_col]]), width = 0.08,
                    linewidth = 0.3, alpha = 0.5) +
      geom_line(linewidth = 0.5) +
      geom_point(size = 1.4) +
      theme_classic(base_family = "Arial")
  }
  p <- p +
    geom_text(aes(y = label_y, label = p_label, vjust = label_vjust), size = 1.9,
              family = "Arial", position = dodge, parse = TRUE) +
    scale_y_continuous(expand = expansion(mult = 0.18)) +
    labs(x = NULL, y = y_lab, title = panel_title) +
    theme(
      text = element_text(size = 7 * FONT_SCALE), axis.title = element_text(size = 7 * FONT_SCALE),
      axis.text = element_text(size = 6), axis.line = element_line(linewidth = 0.4),
      plot.title = element_text(size = 7 * FONT_SCALE, face = "bold"),
      plot.margin = margin(2, 2, 2, 2, unit = "mm"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4)
    )
  if (!show_x_text) p <- p + theme(axis.text.x = element_blank())
  p
}

TOP_ROW <- c("tss_count", "enhancer_length", "ppi_connectivity")
feature_panels <- lapply(names(FEATURE_LABELS), function(f) {
  make_panel(hl_df[hl_df$feature == f, ], "hl_shift", "ci_lower", "ci_upper",
             FEATURE_LABELS[[f]], show_x_text = !(f %in% TOP_ROW))
})
names(feature_panels) <- names(FEATURE_LABELS)
tata_panel <- make_panel(tata_df, "log2_or", "log2_ci_lower", "log2_ci_upper",
                         "TATA box enrichment", y_lab = expression(log[2]("Odds ratio")),
                         by_set = TRUE)
y_label_grob <- wrap_elements(grid::textGrob("Hodges-Lehmann shift (low - high)", rot = 90,
                              gp = grid::gpar(fontsize = 7 * FONT_SCALE, fontfamily = "Arial")))

combined <- wrap_plots(list(L = y_label_grob,
                            A = feature_panels[["tss_count"]], B = feature_panels[["enhancer_length"]],
                            C = feature_panels[["ppi_connectivity"]], D = feature_panels[["twn_hub_rank"]],
                            E = feature_panels[["intron_count"]], G = tata_panel),
                       design = "\nLABC\nLDEG\n") +
  plot_layout(widths = c(0.05, 1, 1, 1)) +
  plot_annotation(
    caption = "Clustering resolution",
    theme = theme(plot.caption = element_text(size = 7 * FONT_SCALE, family = "Arial", hjust = 0.5)))

out_pdf <- file.path(OUT_DIR, "hl_shift_by_resolution.pdf")
ggsave(out_pdf, combined, width = 8, height = 6, units = "in", device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", out_pdf), combined, width = 8, height = 6, units = "in",
       device = ragg::agg_png, dpi = 300)

write.csv(rbind(
  data.frame(key = sprintf("sc_retained_n_%s_res_%s", rep(c("low", "high"), each = 3), RESOLUTIONS),
             value = c(size_df$n_low_retained, size_df$n_high_retained)),
  data.frame(key = sprintf("sc_retained_hl_%s_res_%s", hl_df$feature, hl_df$resolution),
             value = signif(hl_df$hl_shift, 4)),
  data.frame(key = sprintf("sc_retained_hl_p_%s_res_%s", hl_df$feature, hl_df$resolution),
             value = signif(hl_df$p_value, 4))) |>
  transform(provenance = "figure_reproduction/scripts/44_subcluster_retained_genesets_hl_by_resolution.R"),
  file.path(FR, "outputs", "reportable_values_44.csv"), row.names = FALSE)

log_msg("DONE. Wrote ", out_pdf, " and .png")
