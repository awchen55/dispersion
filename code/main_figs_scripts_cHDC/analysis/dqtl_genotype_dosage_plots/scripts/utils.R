# Rownames format: GENE_CHR:POS:REF_ALT (gene symbols contain no
# underscores) -- same convention as dqtl_plots.R::summarize_mash and
# dqtl_egene_overlap's gene_of(). Returns NULL if `gene` isn't in lfsr_mat.
resolve_gene_snp <- function(lfsr_mat, gene) {
  rn_genes <- sub("_[0-9]+:.*$", "", rownames(lfsr_mat))
  stopifnot("mash rownames are not 1 row per gene" = !anyDuplicated(rn_genes))
  hit <- which(rn_genes == gene)
  if (length(hit) == 0) return(NULL)
  snp_full <- rownames(lfsr_mat)[hit]
  snp_part <- sub(paste0("^", gene, "_"), "", snp_full)
  pieces <- regmatches(snp_part, regexec("^([0-9]+):([0-9]+):([A-Z]+)_([A-Z]+)$", snp_part))[[1]]
  stopifnot("Could not parse CHR:POS:REF_ALT from mash rowname" = length(pieces) == 5)
  list(chr = pieces[2], pos = pieces[3], ref = pieces[4], alt = pieces[5],
       geno_snp_id = paste0("chr", pieces[2], ":", pieces[3], ":", pieces[4], "_", pieces[5]))
}

# Genotype dosage (0/1/2 ALT-allele copies) for one SNP, via a single awk
# pass over the relevant chromosome's inputs/genotypes/{chr}_snp.txt --
# avoids loading the whole (~500MB-1GB) file into R for one row.
lookup_genotype_dosage <- function(geno_dir, chr, geno_snp_id) {
  geno_f <- file.path(geno_dir, paste0(chr, "_snp.txt"))
  stopifnot(file.exists(geno_f))
  geno_header <- strsplit(readLines(geno_f, n = 1), "\t")[[1]]

  awk_cmd <- sprintf("awk -F'\t' -v id=%s '$1==id{print;exit}' %s",
                      shQuote(geno_snp_id), shQuote(geno_f))
  geno_line <- system(awk_cmd, intern = TRUE)
  if (length(geno_line) == 0) {
    stop("Genotype row not found for ", geno_snp_id, " in ", geno_f,
         " -- inputs/genotypes/ may have changed since this mash was built")
  }
  geno_values <- strsplit(geno_line, "\t")[[1]]
  stopifnot("Genotype row length does not match header" = length(geno_values) == length(geno_header))
  names(geno_values) <- geno_header
  geno_values <- geno_values[-1]  # drop the leading "id" column
  dosage <- suppressWarnings(as.numeric(geno_values))
  names(dosage) <- names(geno_values)
  stopifnot("Non-finite genotype dosage value(s) for this SNP" = all(is.finite(dosage)))
  dosage
}

# Cohort for one cell type: restrict to exactly the individuals dqtl_mapping
# tested (source_run_id's metadata), not every donor in Erik's raw matrices.
get_cohort <- function(dqtl_mapping_dir, source_run_id, cell_type) {
  meta_header <- read.table(file.path(dqtl_mapping_dir, "inputs/metadata", source_run_id,
                                       paste0(cell_type, "_metadata_sex_pcs.tsv")),
                             sep = "\t", header = TRUE, nrows = 0, check.names = FALSE)
  sample_ids <- setdiff(colnames(meta_header), "id")
  suffix <- paste0("_", cell_type)
  stopifnot("Sample id does not end in the expected _{cell_type} suffix" =
              all(endsWith(sample_ids, suffix)))
  donor_batch <- substr(sample_ids, 1, nchar(sample_ids) - nchar(suffix))
  donor       <- sub("_B[0-9]+$", "", donor_batch)
  mat_cols    <- paste0("sg^", donor_batch, "^", cell_type)
  list(donor_batch = donor_batch, donor = donor, mat_cols = mat_cols)
}

# Match donor_batch keys against a raw-count matrix's Control-column keys,
# falling back to an individual-only match when a specific batch is absent
# (e.g. NA19114_B5 has no raw-count entry but NA19114_B1 does -- same donor,
# same rationale as 01_prepare_dispersion_inputs.R's identical fallback).
match_with_individual_fallback <- function(donor_batch, control_ind_batch) {
  idx <- match(donor_batch, control_ind_batch)
  fallback <- is.na(idx)
  if (any(fallback)) {
    donor_only   <- sub("_B[0-9]+$", "", donor_batch)
    control_only <- sub("_B[0-9]+$", "", control_ind_batch)
    idx[fallback] <- match(donor_only[fallback], control_only)
  }
  idx
}

dosage_colors <- c("0" = "#deebf7", "1" = "#9ecae1", "2" = "#3182bd")

# x-axis tick labels as genotype calls (REF/REF, REF/ALT, ALT/ALT) instead of
# raw dosage counts, same left-to-right ALT-allele-increasing order as the
# "0"/"1"/"2" factor levels.
make_dosage_labels <- function(ref, alt) {
  c("0" = paste0(ref, "/", ref), "1" = paste0(ref, "/", alt), "2" = paste0(alt, "/", alt))
}

dosage_base_theme <- function() {
  theme_classic(base_family = "Arial") +
    theme(
      text = element_text(size = 6), axis.text = element_text(size = 6),
      axis.title = element_text(size = 6), plot.title = element_text(size = 6),
      plot.subtitle = element_text(size = 6), legend.position = "none",
      axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
      plot.margin = margin(2, 2, 2, 2, unit = "mm")
    )
}

# One boxplot + jitter panel: `value_col` vs. dosage, for one (gene, cell
# type). `df` must have a "dosage" factor column (levels "0"/"1"/"2") and a
# numeric column named `value_col`.
dosage_boxplot_panel <- function(df, value_col, y_lab, title, ref, alt) {
  ggplot(df, aes(x = dosage, y = .data[[value_col]], fill = dosage)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.4) +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.6) +
    scale_fill_manual(values = dosage_colors, drop = FALSE) +
    scale_x_discrete(drop = FALSE, labels = make_dosage_labels(ref, alt)) +
    labs(x = NULL, y = y_lab, title = title) +
    dosage_base_theme()
}

# matrixEQTL beta/SE for one (gene, SNP) pair in one cell type, via a single
# awk pass over the per-chromosome matrixEQTL output (much smaller than the
# AllChroms file -- ~30MB vs ~700MB -- since the SNP's chromosome is already
# known). SE derived from beta/t-stat (matrixEQTL doesn't write SE directly);
# sign(t-stat) always matches sign(beta) since t = beta/SE and SE > 0, so
# beta/t_stat is always positive.
#
# Returns NULL (not an error) if the row is absent: 04_map_dqtls.R applies a
# per-cell-type MAF > 0.05 filter (each cell type has its own cohort) BEFORE
# running matrixEQTL, so a gene's single representative mash SNP can pass
# that filter in some cell types and not others even though mash's own
# PosteriorMean/PosteriorSD/lfsr matrices are complete for every gene x
# condition cell (05_reformat_for_mash.R / fastqtl2mash does not re-derive
# its association stats from this same on-disk, already-MAF-filtered
# matrixEQTL text output -- confirmed empirically 2026-08-21: CAMK2D's SNP
# has no row in Cardiac_fibroblasts's matrixEQTL output despite a complete
# mash entry for that gene x cell type).
lookup_matrixeqtl_effect <- function(dqtl_mapping_dir, run_id, cell_type, chr, gene, geno_snp_id) {
  f <- file.path(dqtl_mapping_dir, "outputs/dqtl_results", run_id, cell_type,
                 paste0("dQTL_", cell_type, "_chr", chr, ".txt"))
  stopifnot(file.exists(f))
  awk_cmd <- sprintf("awk -F'\t' -v snp=%s -v gene=%s '$1==snp && $2==gene {print; exit}' %s",
                      shQuote(geno_snp_id), shQuote(gene), shQuote(f))
  line <- system(awk_cmd, intern = TRUE)
  if (length(line) == 0) return(NULL)
  parts <- strsplit(line, "\t")[[1]]
  beta   <- as.numeric(parts[3])
  t_stat <- as.numeric(parts[4])
  stopifnot("Non-finite matrixEQTL beta/t-stat" = is.finite(beta), is.finite(t_stat))
  list(beta = beta, se = beta / t_stat)
}

# Forest plot: one row per condition (fixed order via `cell_type_order`,
# first entry at top), point = effect size, horizontal bar = 95% Wald CI
# (est +/- 1.96*se), dashed reference line at 0 -- mirrors Fig 4 of the mash
# paper (Urbut et al. 2019, Nat Genet). `df` needs columns cell_type, est, se,
# is_dgene (logical, used to color the point/range). `df` may have fewer than
# length(cell_type_order) rows (e.g. matrixEQTL has no row for a condition
# where the SNP failed that cohort's own MAF filter) -- drop = FALSE keeps
# every condition's row position on the axis regardless, so that missing
# condition renders as a visibly empty row rather than silently compressing
# the layout, and every gene's figure keeps the same 4-row axis.
#
# Optional `has_real_data` logical column (used by the mash version): when a
# condition's matrixEQTL row is missing, 05_reformat_for_mash.R placeholder-
# fills mash's own input for that (gene, SNP, condition) cell with beta=0,
# se=1e16 (near-zero information) before mash ever sees it -- so that
# condition's mash posterior is almost entirely PRIOR-driven (mash borrowing
# strength from the covariance structure learned across other genes/
# conditions), not evidence from this condition's own data. Rendered as an
# open point (vs. filled for a condition with a real matrixEQTL observation)
# so this is visible rather than looking identical to a genuinely-observed
# estimate. Absent column => treated as all real (matrixEQTL's own plot).
effect_forest_plot <- function(df, cell_type_order, x_lab, title, subtitle) {
  df$cell_type <- factor(df$cell_type, levels = cell_type_order)
  df$ci_lo <- df$est - 1.96 * df$se
  df$ci_hi <- df$est + 1.96 * df$se
  if (!"has_real_data" %in% names(df)) df$has_real_data <- TRUE
  ggplot(df, aes(y = cell_type, x = est, xmin = ci_lo, xmax = ci_hi,
                 color = is_dgene, shape = has_real_data)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey50") +
    geom_pointrange(orientation = "y", size = 0.3, linewidth = 0.4) +
    scale_y_discrete(limits = rev(cell_type_order), drop = FALSE,
                      labels = function(x) gsub("_", " ", x)) +
    scale_color_manual(values = c("TRUE" = "#3182bd", "FALSE" = "grey50"), guide = "none") +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), guide = "none") +
    labs(x = x_lab, y = NULL, title = title, subtitle = subtitle) +
    dosage_base_theme() +
    theme(plot.subtitle = element_text(size = 6, family = "Arial"))
}
