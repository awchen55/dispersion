library(variancePartition)
library(limma)
library(BiocParallel)

# Run from the repo root (analysis/dd_sweep_reviewer_response/), as run.sh does.
source("scripts/00_config.R")

OUT_DIR <- "outputs/brendan_full_reference_no_random_effect"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

load_matrix <- function(prefix) {
  mats <- lapply(BRENDAN_CELL_TYPES, function(ct) {
    read.csv(
      file.path(BRENDAN_FULL_MEMENTO_DIR, paste0(prefix, "_by_ind_", ct, ".csv")),
      row.names = 1, check.names = FALSE
    )
  })
  common_genes <- Reduce(intersect, lapply(mats, rownames))
  mats <- lapply(mats, function(m) m[common_genes, , drop = FALSE])
  do.call(cbind, mats)
}

res_var_mat <- load_matrix("dispersion")
message("Samples: ", ncol(res_var_mat), " | Genes before filtering: ", nrow(res_var_mat))
stopifnot(ncol(res_var_mat) == 3 * length(BRENDAN_CELL_TYPES))  # 3 individuals x 6 cell types = 18

n_before <- nrow(res_var_mat)
res_var_mat <- res_var_mat[apply(res_var_mat, 1, function(x) all(is.finite(x))), , drop = FALSE]
n_dropped <- n_before - nrow(res_var_mat)
message("Removed ", n_dropped, " gene(s) with non-finite dispersion values (",
        round(100 * n_dropped / n_before, 2), "%)")
stopifnot(n_dropped < 5)  # full-depth data; expect ~0 non-finite genes, unlike the low-cell-count sweep points

res_var_norm <- normalizeQuantiles(as.matrix(res_var_mat))

split_parts <- strsplit(colnames(res_var_norm), "_")
individuals <- sapply(split_parts, `[`, 1)
cell_types  <- sapply(split_parts, function(x) paste(x[-1], collapse = "_"))

metadata <- data.frame(
  cell_type = factor(cell_types),
  cell_line = factor(individuals),
  row.names = colnames(res_var_norm)
)
stopifnot(nlevels(metadata$cell_line) == 3, nlevels(metadata$cell_type) == length(BRENDAN_CELL_TYPES))

celltypes <- levels(metadata$cell_type)
pairs <- combn(celltypes, 2, simplify = FALSE)
contrast_matrix <- setNames(
  sapply(pairs, function(x) paste0("cell_type", x[2], " - ", "cell_type", x[1])),
  sapply(pairs, function(x) paste0("compare_", x[2], "_", x[1]))
)

# The one deliberate change from the manuscript pipeline's model: no
# (1|cell_line) random effect, matching every fit in the downsampled sweep.
formula <- ~ 0 + cell_type
L <- makeContrastsDream(formula = formula, data = metadata, contrasts = contrast_matrix)

param <- SnowParam(4, type = "SOCK", progressbar = TRUE, exportglobals = FALSE)
fit <- dream(res_var_norm, formula, metadata, L, BPPARAM = param)
fit <- eBayes(fit)

out_file <- file.path(OUT_DIR, "dd_fit_brendan_full_no_random_effect.rds")
saveRDS(fit, out_file)
message("Saved: ", out_file)

# Sanity check: confirm all 4 shared contrasts used elsewhere in this repo
# are present in this fit's coefficients.
stopifnot(all(SHARED_CONTRASTS$contrast %in% colnames(fit$coefficients)))
message("All ", nrow(SHARED_CONTRASTS), " shared contrasts confirmed present.")
