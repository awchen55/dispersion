library(tidyverse)

BASE       <- "analysis/dqtl_mapping_q02/outputs/cell_count_threshold_sweep_q02"
THRESHOLDS <- c(150)
CELL_TYPES <- c("Activated_fibroblasts", "Cardiac_fibroblasts", "Cardiomyocytes",
                "Epithelial_cells", "Foregut_endodermal_cells", "Mesenchymal_cells")

# Rename columns from "sg^NA18510_B1^Activated_fibroblasts" to "NA18510_B1_Activated_fibroblasts"
# Batch is retained so it can be used as a covariate in downstream models.
# (Identical to reformat_sweep_dataframes.R's rename_cols() -- not re-derived.)
rename_cols <- function(df) {
  rownames(df) <- df$gene
  df$gene <- NULL
  parts       <- strsplit(colnames(df), "\\.")
  ind_batch   <- sapply(parts, `[`, 2)   # e.g. "NA18510_B1"
  cell_types  <- sapply(parts, `[`, 3)
  colnames(df) <- paste(ind_batch, cell_types, sep = "_")
  df
}

for (thresh in THRESHOLDS) {
  dir_path <- file.path(BASE, paste0("cell_count_", thresh))
  stopifnot("cell_count_{thresh} directory not found -- run 01_memento_q02_threshold150.py first" =
              dir.exists(dir_path))
  for (ct in CELL_TYPES) {
    for (prefix in c("dispersion", "mean")) {
      in_file  <- file.path(dir_path, paste0(prefix, "_by_ind_", ct, ".csv"))
      out_file <- file.path(dir_path, paste0("new_", prefix, "_by_ind_", ct, ".csv"))
      stopifnot(file.exists(in_file))
      df <- read.csv(in_file, row.names = 1)
      df <- rename_cols(df)
      write.csv(df, out_file, row.names = TRUE)
      message("Wrote: ", out_file)
    }
  }
}
