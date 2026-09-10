library(Seurat)
library(SeuratObject)
library(Matrix)

seurat_path         <- snakemake@input[["seurat"]]
out_path            <- snakemake@output[["pb"]]
min_cells_per_group <- snakemake@params[["min_cells_per_group"]]

generate_pseudobulk <- function(object, labels, assay = "RNA", slot = "counts") {
  factorlist <- lapply(labels, function(l) unique(object@meta.data[, l]))
  names(factorlist) <- labels
  meta <- expand.grid(factorlist, stringsAsFactors = FALSE)
  rownames(meta) <- apply(meta, 1, paste0, collapse = ".")

  counts_mat <- slot(object[[assay]], slot)
  n <- nrow(meta)
  out <- matrix(0L, nrow = nrow(counts_mat), ncol = n,
                dimnames = list(rownames(counts_mat), rownames(meta)))

  ncells <- integer(n)
  for (i in seq_len(n)) {
    sel <- seq_len(ncol(counts_mat))
    for (j in labels) sel <- sel[object@meta.data[[j]][sel] == meta[i, j]]
    ncells[i] <- length(sel)
    if (length(sel) == 1L) out[, i] <- counts_mat[, sel]
    else if (length(sel) > 1L) out[, i] <- Matrix::rowSums(counts_mat[, sel, drop = FALSE])
  }
  meta$ncells <- ncells
  list(counts = out, meta = meta)
}

add_sex <- function(bulk) {
  bulk$meta$sex <- "F"
  male_ids <- c("NA19160", "NA19210", "NA18913", "NA18519", "C3649", "C8861", "C4955")
  bulk$meta$sex[bulk$meta$individual %in% male_ids] <- "M"
  bulk
}

seurat_obj <- readRDS(seurat_path)
seurat_obj[["RNA"]] <- as(seurat_obj[["RNA"]], "Assay")

pb <- generate_pseudobulk(seurat_obj, labels = c("donor_id", "species"))
# Same threshold as memento_diploid.py's per-donor cell-count filter: a
# donor x species pseudobulk sample built from too few cells is a noisy,
# low-depth input to the DE model, not just an empty one. Applying the same
# cutoff here keeps the DE and DD pipelines consistent about which
# donor/cell-type combinations are excluded (see q_threshold_sensitivity.Rmd
# Section 8b for the original diagnosis of this issue in memento).
keep <- pb$meta$ncells >= min_cells_per_group
pb$counts <- pb$counts[, keep, drop = FALSE]
pb$meta   <- pb$meta[keep, ]
pb <- add_sex(pb)

saveRDS(pb, out_path)
