#!/usr/bin/env Rscript

# Column names ({IND}_{CT}) for one cell type across the three individuals.
ct_columns <- function(cell_type) paste0(individuals, "_", cell_type)

# Jaccard index of two character vectors (NA-safe: empty union -> NA).
jaccard <- function(a, b) {
  if (length(union(a, b)) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}
