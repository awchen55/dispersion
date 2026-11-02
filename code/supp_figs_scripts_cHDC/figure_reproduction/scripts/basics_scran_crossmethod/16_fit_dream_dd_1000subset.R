stopifnot("must run from project root (CLAUDE.md not found)" = file.exists("CLAUDE.md"))

suppressPackageStartupMessages({
  library(limma)
  library(variancePartition)
  library(BiocParallel)
})

# ---- named constants (match 03_rerun_dd_rank_analyses.R) -------------------
CELL_TYPES <- c("Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
                "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells")
INDIVIDUALS <- c("NA19099", "NA19093", "NA18522")
N_EXPECTED_COLS <- length(INDIVIDUALS) * length(CELL_TYPES)  # 18
RIBO_PATTERN <- "^RP"
N_DREAM_CORES <- 4L
MIN_GENES <- 5000L  # canary: ~11-12k genes expected; far fewer means the all-18-columns intersection broke
# ANALYSIS_OK[model-formula]: dream design formula (cell_type fixed effect, no intercept), identical to 03's published DD model; no numeric weights.
DREAM_FORMULA <- ~ 0 + cell_type

IN_DIR  <- "figure_reproduction/outputs/basics_scran_crossmethod/memento_1000subset"
OUT_RDS <- file.path(IN_DIR, "dream_fit_1000subset.rds")

mats <- lapply(CELL_TYPES, function(ct) {
  f <- file.path(IN_DIR, paste0("dispersion_by_ind_", ct, ".csv"))
  stopifnot("1000-subset Memento CSV not found" = file.exists(f))
  m <- read.csv(f, row.names = 1, check.names = FALSE)
  stopifnot("expected 3 individual columns" = ncol(m) == 3,
            "unexpected column names" = setequal(colnames(m), paste0(INDIVIDUALS, "_", ct)),
            "duplicate gene names" = !anyDuplicated(rownames(m)))
  message(sprintf("  %s: %d genes (%d with non-finite value)", ct, nrow(m), sum(!is.finite(rowSums(as.matrix(m))))))
  as.matrix(m)
})
names(mats) <- CELL_TYPES

# ANALYSIS_OK[filtering]: genes must have finite values in all 6 cell types for every pairwise contrast to be estimable; kept/dropped counts logged below and floor-checked by MIN_GENES.
shared <- Reduce(intersect, lapply(mats, function(m) rownames(m)[is.finite(rowSums(m))]))
message(sprintf("[load] union of genes = %d; genes finite in all 18 columns = %d",
                length(unique(unlist(lapply(mats, rownames)))), length(shared)))
# ANALYSIS_OK[filtering]: drop ^RP ribosomal genes (repo convention, same as 03); count logged on next line.
shared <- shared[!grepl(RIBO_PATTERN, shared)]
message(sprintf("[load] after dropping ^RP: %d genes", length(shared)))
stopifnot("too few genes after intersection" = length(shared) >= MIN_GENES)

res_var_mat <- do.call(cbind, lapply(mats, function(m) m[shared, , drop = FALSE]))
stopifnot("expected 18 columns" = ncol(res_var_mat) == N_EXPECTED_COLS,
          "non-finite values remain" = all(is.finite(res_var_mat)))

sample_cols <- colnames(res_var_mat)
metadata <- data.frame(
  cell_type = factor(sub("^[^_]*_", "", sample_cols), levels = CELL_TYPES),
  cell_line = factor(sub("_.*$", "", sample_cols), levels = INDIVIDUALS),
  row.names = sample_cols)
stopifnot("unparsed sample names" = !anyNA(metadata$cell_type) && !anyNA(metadata$cell_line))

pairs <- combn(CELL_TYPES, 2, simplify = FALSE)
contrast_matrix <- setNames(
  vapply(pairs, function(x) paste0("cell_type", x[2], " - cell_type", x[1]), character(1)),
  vapply(pairs, function(x) paste0("compare_", x[2], "_", x[1]), character(1)))

mat <- normalizeQuantiles(res_var_mat)
contrast_dream <- makeContrastsDream(DREAM_FORMULA, metadata, contrasts = contrast_matrix)
param <- SnowParam(N_DREAM_CORES, type = "SOCK", progressbar = FALSE, exportglobals = FALSE)
message(sprintf("dream fit: %d genes x %d samples, %d contrasts", nrow(mat), ncol(mat), length(contrast_matrix)))
fit <- eBayes(dream(mat, DREAM_FORMULA, metadata, contrast_dream, BPPARAM = param))
stopifnot("expected 15 compare_ contrasts" = sum(grepl("^compare_", colnames(fit$coefficients))) == 15L)
saveRDS(fit, OUT_RDS)
message("Saved: ", OUT_RDS)
