#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)      # read_table, as in the source
  library(ggplot2)
  library(patchwork)
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
fr_dir <- file.path(repo_dir, "figure_reproduction")
cluster_dir <- file.path(repo_dir, "cellranger_cluster-mode_trial")
res_var_csv <- file.path(fr_dir, "data", "res_var_excl_cluster24.csv")
retained_matrix_csv <- file.path(fr_dir, "outputs", "retained_gene_matrix.csv")
pub_res_var_csv <- file.path(cluster_dir, "analysis", "datasets", "memento",
                             "new_res_var_by_cell_type_ind_1000_cts_new_predicted_id.csv")
# ANALYSIS_OK[external-path]: read-only path to Erik's frozen MatrixEQTL results (same as script 03 Part A)
eqtl_file_tmpl <- file.path("/project/gilad/emcintire/chemo_project_2023/pipeline/map_eqtls",
                            "archive", "no_cvrt_condox_only", "eqtl_results", "Control",
                            "%s", "steadyState_%s_Control_AllChroms.txt")
lfsr_csv <- file.path(cluster_dir, "analysis", "datasets", "erik_egenes", "mash_qtl_lfsr_05")
beta_csv <- file.path(cluster_dir, "analysis", "datasets", "erik_egenes", "premash_betas")
out_dir <- file.path(fr_dir, "outputs", "figS8_S9_egene_sliding_window")
reportable_csv <- file.path(fr_dir, "outputs", "reportable_values_45.csv")

individuals <- c("NA19099", "NA19093", "NA18522")
ribo_pattern <- "^RP"
seed <- 123            # source: set.seed(123) per cell type
n_permutations <- 15   # source value
step_size <- 100
s8 <- list(name = "matrixeqtl", window = 1500, padj = 0.05,
           cell_types = c("Cardiomyocytes", "Cardiac_fibroblasts", "Activated_fibroblasts",
                          "Foregut_endodermal_cells"),
           ylim = c(0, 0.02), ncol = 2, width = 4.4, height = 3)
s9 <- list(name = "mash", window = 3000, lfsr = 0.05,
           cell_types = c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                          "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells"),
           ncol = 3, width = 5.4, height = 3)

log_msg <- function(...) message(sprintf("[45] %s", sprintf(...)))

for (f in c(res_var_csv, retained_matrix_csv, pub_res_var_csv, lfsr_csv, beta_csv)) {
  if (!file.exists(f)) stop("Missing input: ", f)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------- eGene sets -----------------------------------
matrixeqtl_egenes <- function(cell_type) {
  f <- sprintf(eqtl_file_tmpl, cell_type, cell_type)
  if (!file.exists(f)) stop("Missing MatrixEQTL file: ", f)
  eq <- suppressMessages(read_table(f))
  if (!all(c("gene", "bh_adj_pvalue") %in% colnames(eq))) stop("Unexpected columns in ", f)
  unique(eq$gene[eq$bh_adj_pvalue <= s8$padj])
}
lfsr <- read.csv(lfsr_csv, row.names = 1)
premash_beta <- read.csv(beta_csv, row.names = 1)
mash_egenes <- function(cell_type) {
  col <- paste0(cell_type, "_Control_lfsr")
  if (!col %in% colnames(lfsr)) stop("Missing lfsr column ", col)
  rownames(lfsr)[lfsr[[col]] < s9$lfsr]
}
mash_beta <- function(cell_type) {
  col <- paste0(cell_type, "_Control_premash_beta")
  if (!col %in% colnames(premash_beta)) stop("Missing beta column ", col)
  setNames(premash_beta[[col]], rownames(premash_beta))
}

# ------------------------------- Window function -------------------------------
# ranks: named vector of dispersion ranks over the gene universe (names in
# matrix order, as rownames(rank_df) in the source). egenes: eGenes within the
# universe. beta: optional named vector (mash premash betas) for the median
# |beta| column. Draws the null with the same RNG calls as the source.
sliding_window <- function(ranks, egenes, window_size, beta = NULL) {
  starts <- seq(1, length(ranks) - window_size + 1, by = step_size)
  ends <- starts + window_size - 1
  members <- lapply(seq_along(starts), function(w) names(ranks)[ranks >= starts[w] & ranks <= ends[w]])
  n_in <- lengths(members)
  prop_of <- function(set) vapply(members, function(g) if (length(g) > 0) sum(g %in% set) / length(g) else NA_real_, numeric(1))
  res <- data.frame(start = starts, end = ends, genes_in_window = n_in,
                    egenes_in_window = vapply(members, function(g) sum(g %in% egenes), integer(1)),
                    proportion_egenes = prop_of(egenes))
  if (!is.null(beta)) {
    eg_beta <- beta[names(beta) %in% egenes]
    res$median_beta_in_window <- vapply(members, function(g) median(abs(eg_beta[names(eg_beta) %in% g])), numeric(1))
  }
  res$midpoint <- (res$start + res$end) / 2
  set.seed(seed)
  null_matrix <- t(vapply(seq_len(n_permutations), function(i) {
    prop_of(sample(names(ranks), length(egenes), replace = FALSE))
  }, numeric(length(starts))))
  res$null_mean <- colMeans(null_matrix, na.rm = TRUE)
  res$null_lower <- apply(null_matrix, 2, quantile, probs = 0.025, na.rm = TRUE)
  res$null_upper <- apply(null_matrix, 2, quantile, probs = 0.975, na.rm = TRUE)
  res
}

median_ranks <- function(res_var_df, cell_type, genes) {
  cols <- paste0(individuals, "_", cell_type)
  if (!all(cols %in% colnames(res_var_df))) stop("Missing columns for ", cell_type)
  med <- apply(res_var_df[genes, cols, drop = FALSE], 1, median)
  rank(med, ties.method = "average")
}

# ---------------- Check: reproduce the source loop on the published matrix -----
# Verbatim transcription of the source loop body (Supp S8 / S9), apart from
# taking res_var_df / egenes / beta as arguments.
source_loop <- function(res_var_df, cell_type, egenes, window_size, ct_beta = NULL, beta_col = NULL) {
  cols <- paste0(c("NA19099_", "NA19093_", "NA18522_"), cell_type)
  ct_res_var <- res_var_df %>% dplyr::select(all_of(cols))
  median_res_var <- data.frame(med_vals = apply(ct_res_var, 1, function(x) median(x)))
  egenes <- egenes[egenes %in% rownames(ct_res_var)]
  rank_df <- apply(median_res_var, 2, function(x) rank(x, ties.method = "average"))
  results <- data.frame()
  max_rank <- nrow(rank_df)
  for (start in seq(1, max_rank - window_size + 1, by = step_size)) {
    end <- start + window_size - 1
    in_window <- apply(rank_df, 1, function(ranks) all(ranks >= start & ranks <= end))
    genes_in_all <- rownames(rank_df)[in_window]
    num_egenes <- sum(genes_in_all %in% egenes)
    proportion <- if (length(genes_in_all) > 0) num_egenes / length(genes_in_all) else NA
    row <- data.frame(start = start, end = end, proportion_egenes = proportion,
                      genes_in_window = length(genes_in_all), egenes_in_window = num_egenes)
    if (!is.null(ct_beta)) row$median_beta_in_window <- median(abs(ct_beta[[beta_col]][ct_beta$gene %in% genes_in_all]))
    results <- rbind(results, row)
  }
  set.seed(123)
  n_windows <- nrow(results)
  null_matrix <- matrix(NA, nrow = n_permutations, ncol = n_windows)
  all_genes <- rownames(rank_df)
  for (i in 1:n_permutations) {
    random_egenes <- sample(all_genes, length(egenes), replace = FALSE)
    tmp_results <- c()
    for (w in seq_len(n_windows)) {
      in_window <- apply(rank_df, 1, function(ranks) all(ranks >= results$start[w] & ranks <= results$end[w]))
      genes_in_all <- rownames(rank_df)[in_window]
      prop <- if (length(genes_in_all) > 0) sum(genes_in_all %in% random_egenes) / length(genes_in_all) else NA
      tmp_results <- c(tmp_results, prop)
    }
    null_matrix[i, ] <- tmp_results
  }
  results$null_mean <- colMeans(null_matrix, na.rm = TRUE)
  results$null_lower <- apply(null_matrix, 2, quantile, probs = 0.025, na.rm = TRUE)
  results$null_upper <- apply(null_matrix, 2, quantile, probs = 0.975, na.rm = TRUE)
  results
}

pub_res_var <- read.csv(pub_res_var_csv, row.names = 1)
# ANALYSIS_OK[filtering]: ^RP removal as in the source setup chunk (lines 47-52); check-only input
pub_res_var <- pub_res_var[!grepl(ribo_pattern, rownames(pub_res_var)), ]
check_cols <- c("start", "end", "genes_in_window", "egenes_in_window", "proportion_egenes",
                "null_mean", "null_lower", "null_upper")
for (spec in list(s8, s9)) {
  for (ct in spec$cell_types) {
    egenes <- if (spec$name == "matrixeqtl") matrixeqtl_egenes(ct) else mash_egenes(ct)
    ct_beta <- NULL; beta <- NULL
    if (spec$name == "mash") {
      beta <- mash_beta(ct)
      ct_beta <- data.frame(gene = names(beta), b = unname(beta))
      ct_beta <- ct_beta[ct_beta$gene %in% egenes, ]
      colnames(ct_beta)[2] <- paste0(ct, "_Control_premash_beta")
    }
    ref <- source_loop(pub_res_var, ct, egenes, spec$window, ct_beta, paste0(ct, "_Control_premash_beta"))
    genes <- rownames(pub_res_var)
    mine <- sliding_window(median_ranks(pub_res_var, ct, genes), egenes[egenes %in% genes], spec$window, beta)
    cols <- c(check_cols, if (spec$name == "mash") "median_beta_in_window")
    if (!isTRUE(all.equal(mine[, cols], ref[, cols], tolerance = 0, check.attributes = FALSE))) {
      stop(sprintf("Port check failed: %s %s does not reproduce the source loop", spec$name, ct))
    }
  }
  log_msg("Check passed: %s windows reproduce the source loop exactly (published matrix, %d cell types)",
          spec$name, length(spec$cell_types))
}

# ---------------- Run on the excl-cluster-24 matrix, per-CT retained genes ------
res_var <- read.csv(res_var_csv, row.names = 1, check.names = FALSE)
# ANALYSIS_OK[filtering]: ^RP removal (repo convention, as script 01); count logged below
res_var <- res_var[!grepl(ribo_pattern, rownames(res_var)), ]
retained_ind <- as.matrix(read.csv(retained_matrix_csv, row.names = 1, check.names = FALSE))
if (!all(rownames(retained_ind) %in% rownames(res_var))) stop("retained_gene_matrix has genes absent from the matrix")
log_msg("%d non-ribo genes in the excl-cluster-24 dispersion matrix", nrow(res_var))

plot_panel <- function(results, ct, spec) {
  if (spec$name == "matrixeqtl") {
    plot <- ggplot(results, aes(x = midpoint, y = proportion_egenes)) +
      geom_line(color = "red", size = 0.6) +
      geom_ribbon(aes(ymin = null_lower, ymax = null_upper), fill = "gray80", alpha = 0.5) +
      geom_line(aes(y = null_mean), linetype = "dashed", color = "black") +
      ylim(spec$ylim[1], spec$ylim[2]) +
      labs(x = NULL, y = NULL, title = paste0(ct))
    y_cts <- c("Cardiomyocytes", "Activated_fibroblasts")
    x_cts <- c("Activated_fibroblasts", "Foregut_endodermal_cells")
  } else {
    plot <- ggplot(results, aes(x = midpoint)) +
      geom_line(aes(y = proportion_egenes), color = "red", linewidth = 0.6) +
      geom_ribbon(aes(ymin = null_lower, ymax = null_upper), fill = "gray80", alpha = 0.5) +
      geom_line(aes(y = null_mean), linetype = "dashed", color = "black") +
      scale_y_continuous(name = "eGene proportion") +
      labs(x = NULL, y = NULL, title = ct)
    y_cts <- c("Cardiomyocytes", "Epithelial_cells")
    x_cts <- c("Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells")
  }
  plot <- plot +
    theme(title = element_text(size = 6)) +
    theme_classic(base_family = "Arial") +
    theme(
      text = element_text(size = 6),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 6),
      plot.title = element_text(size = 6, hjust = 0.5),
      axis.line = element_line(linewidth = 0.4),
      plot.margin = margin(2, 2, 2, 2, unit = "mm"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4))
  if (ct %in% y_cts) {
    plot <- plot + theme(axis.text.y = element_text(size = 6),
                         axis.title.y = element_text(size = 7, angle = 90)) +
      labs(y = "Proportion of eGenes")
  }
  if (ct %in% x_cts) {
    plot <- plot + theme(axis.text.x = element_text(size = 6),
                         axis.title.x = element_text(size = 7)) +
      labs(x = "Rank Window Midpoint")
  }
  plot
}

reportable <- data.frame()
for (spec in list(s8, s9)) {
  all_results <- list(); plots <- list()
  for (ct in spec$cell_types) {
    genes <- rownames(retained_ind)[retained_ind[, ct] == 1]
    egenes_all <- if (spec$name == "matrixeqtl") matrixeqtl_egenes(ct) else mash_egenes(ct)
    egenes <- egenes_all[egenes_all %in% genes]
    beta <- if (spec$name == "mash") mash_beta(ct) else NULL
    res <- sliding_window(median_ranks(res_var, ct, genes), egenes, spec$window, beta)
    log_msg("%s %s: %d retained genes ranked, %d of %d eGenes in that set, %d windows",
            spec$name, ct, length(genes), length(egenes), length(egenes_all), nrow(res))
    if (spec$name == "matrixeqtl") {
      top <- max(c(res$proportion_egenes, res$null_upper), na.rm = TRUE)
      if (top > spec$ylim[2]) stop(sprintf("%s: max proportion %.4f exceeds the source's fixed ylim %.2f", ct, top, spec$ylim[2]))
    }
    all_results[[ct]] <- cbind(cell_type = ct, res)
    plots[[ct]] <- plot_panel(res, ct, spec)
    reportable <- rbind(reportable, data.frame(
      figure = if (spec$name == "matrixeqtl") "S8" else "S9", egene_source = spec$name, cell_type = ct,
      n_genes_ranked = length(genes), n_egenes = length(egenes), n_windows = nrow(res),
      n_windows_above_null_upper = sum(res$proportion_egenes > res$null_upper, na.rm = TRUE),
      n_windows_below_null_lower = sum(res$proportion_egenes < res$null_lower, na.rm = TRUE)))
  }
  out <- do.call(rbind, all_results)
  rownames(out) <- NULL
  write.csv(out, file.path(out_dir, paste0(spec$name, "_sliding_window_results.csv")), row.names = FALSE)
  combined_plot <- wrap_plots(plots) + plot_layout(ncol = spec$ncol, axis_titles = "collect")
  fig_name <- if (spec$name == "matrixeqtl") "FigureS8_matrixeqtl_sliding_window.pdf" else "FigureS9_mash_sliding_window.pdf"
  ggsave(file.path(out_dir, fig_name), combined_plot, device = cairo_pdf,
         width = spec$width, height = spec$height, units = "in")
  log_msg("Wrote %s", fig_name)
}
write.csv(reportable, reportable_csv, row.names = FALSE)
log_msg("Wrote %s", reportable_csv)
