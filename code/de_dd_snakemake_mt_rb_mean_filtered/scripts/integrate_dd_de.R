library(dplyr)
library(tidyr)
library(mashr)

tet_de_paths    <- snakemake@input[["tet_de_csvs"]]
dip_de_logfc    <- snakemake@input[["dip_de_logfc"]]
dip_de_se       <- snakemake@input[["dip_de_se"]]
dip_de_pval     <- snakemake@input[["dip_de_pval"]]
dip_de_adj_pval <- snakemake@input[["dip_de_adj_pval"]]
tet_dd_logfc    <- snakemake@input[["tet_dd_logfc"]]
tet_dd_se       <- snakemake@input[["tet_dd_se"]]
tet_dd_mash     <- snakemake@input[["tet_dd_mash"]]
dip_dd_logfc    <- snakemake@input[["dip_dd_logfc"]]
dip_dd_se       <- snakemake@input[["dip_dd_se"]]
dip_dd_mash     <- snakemake@input[["dip_dd_mash"]]
out_master      <- snakemake@output[["master"]]

ct_dip <- snakemake@params[["celltypes_diploid"]]
ct_tet <- snakemake@params[["celltypes_tetraploid"]]

# ----- helper -----
read_csv_gene_col <- function(path) {
  df <- read.csv(path, check.names = FALSE)
  if (colnames(df)[1] %in% c("X", "")) colnames(df)[1] <- "gene"
  df
}

# ----- tetraploid DE (Wilcoxon) -----
master_df <- NULL
for (i in seq_along(ct_tet)) {
  ct  <- ct_tet[i]
  ct2 <- gsub("_", ".", ct)
  res <- read.csv(tet_de_paths[i])
  colnames(res)[colnames(res) == "X"] <- "gene"

  res <- res %>%
    mutate(adj.p.val = ifelse(obs.paired < 5, NA,
                              p.adjust(pvalue, method = "BH")))
  res <- res[, c("gene", "pvalue", "logFC", "adj.p.val")]
  colnames(res) <- c("gene",
                     paste0("de.pvalue.allotetraploid.", ct2),
                     paste0("de.logFC.allotetraploid.", ct2),
                     paste0("de.adj.pvalue.allotetraploid.", ct2))

  master_df <- if (is.null(master_df)) res else full_join(master_df, res, by = "gene")
}

# ----- diploid DE -----
for (path in c(dip_de_logfc, dip_de_se, dip_de_pval, dip_de_adj_pval)) {
  master_df <- full_join(master_df, read_csv_gene_col(path), by = "gene")
}

# ----- tetraploid DD: logFC + SE from merged matrices -----
for (path in c(tet_dd_logfc, tet_dd_se)) {
  master_df <- full_join(master_df, read_csv_gene_col(path), by = "gene")
}

# ----- tetraploid DD: LFSR from MASH -----
mash_tet  <- readRDS(tet_dd_mash)
lfsrs_tet <- get_lfsr(mash_tet) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene")
colnames(lfsrs_tet) <- gsub("dd\\.logFC\\.allotetraploid\\.",
                             "dd.lfsr.allotetraploid.", colnames(lfsrs_tet))
master_df <- full_join(master_df, lfsrs_tet, by = "gene")

# ----- diploid DD: logFC + SE from filtered matrices -----
for (path in c(dip_dd_logfc, dip_dd_se)) {
  master_df <- full_join(master_df, read_csv_gene_col(path), by = "gene")
}

# ----- diploid DD: LFSR from MASH -----
mash_dip  <- readRDS(dip_dd_mash)
lfsrs_dip <- get_lfsr(mash_dip) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene")
colnames(lfsrs_dip) <- gsub("dd\\.logFC\\.diploid\\.",
                             "dd.lfsr.diploid.", colnames(lfsrs_dip))
master_df <- full_join(master_df, lfsrs_dip, by = "gene")

# Keep only genes actually tested for DD in both diploid and allotetraploid.
# logFC is imputed to 0 for untested genes, so != 0 identifies tested genes.
lfc_dip_cols  <- grep("^dd\\.logFC\\.diploid\\.",        colnames(master_df), value = TRUE)
lfc_allo_cols <- grep("^dd\\.logFC\\.allotetraploid\\.", colnames(master_df), value = TRUE)
tested_dip  <- rowSums(master_df[, lfc_dip_cols]  != 0, na.rm = TRUE) > 0
tested_allo <- rowSums(master_df[, lfc_allo_cols] != 0, na.rm = TRUE) > 0
master_df   <- master_df[tested_dip & tested_allo, ]

# ----- Derived per-cell-type fields (recapitulation, cis/trans) -----
# Restricted to cell types tested in BOTH diploid and allotetraploid.
# ct_dip and ct_tet are now the same 20 cell types (config.yaml), so this
# intersection is just that shared list -- kept as an intersect() rather
# than using ct_dip directly so this stays correct if the two lists are
# ever allowed to diverge again.
ct_common <- gsub("_", ".", intersect(ct_dip, ct_tet))

# Recapitulated DD: significant (lfsr < 0.05) in both diploid and
# allotetraploid, with the same sign of effect size (logFC). 1 = recapitulated.
for (ct in ct_common) {
  lfsr_dip   <- paste0("dd.lfsr.diploid.", ct)
  lfsr_allo  <- paste0("dd.lfsr.allotetraploid.", ct)
  logfc_dip  <- paste0("dd.logFC.diploid.", ct)
  logfc_allo <- paste0("dd.logFC.allotetraploid.", ct)
  new_col    <- paste0("dd.recapitulated.same.sign.", ct)

  master_df[[new_col]] <- as.integer(
    !is.na(master_df[[lfsr_dip]])  & master_df[[lfsr_dip]]  < 0.05 &
    !is.na(master_df[[lfsr_allo]]) & master_df[[lfsr_allo]] < 0.05 &
    sign(master_df[[logfc_dip]]) == sign(master_df[[logfc_allo]])
  )
}

master_df <- master_df %>%
  mutate(dd.recapitulated.same.sign.n.celltypes =
           rowSums(select(., starts_with("dd.recapitulated.same.sign.")), na.rm = TRUE))

# Cis proportion per cell type: the allotetraploid (cis-only) effect size as a
# fraction of the total diploid effect size, where the remainder
# (diploid - allotetraploid) is attributed to trans effects.
for (ct in ct_common) {
  allo_col <- paste0("dd.logFC.allotetraploid.", ct)
  dip_col  <- paste0("dd.logFC.diploid.", ct)
  prop_col <- paste0("dd.cis.proportion.", ct)

  master_df[[prop_col]] <- abs(master_df[[allo_col]]) /
    (abs(master_df[[allo_col]]) + abs(master_df[[dip_col]] - master_df[[allo_col]]))
}

# Mean cis proportion, averaged only over cell types where the gene is DD in
# diploid (lfsr < 0.05); NaN (no DD cell types) is handled below.
master_df <- master_df %>%
  rowwise() %>%
  mutate(
    dd.in.diploid.mean.cis.proportion = mean(
      ifelse(
        c_across(all_of(paste0("dd.lfsr.diploid.", ct_common))) < 0.05,
        c_across(all_of(paste0("dd.cis.proportion.", ct_common))),
        NA_real_
      ),
      na.rm = TRUE
    )
  ) %>%
  ungroup()

# Cis/trans category per cell type, only among diploid-DD genes:
# cis proportion > 0.70 -> cis, < 0.30 -> trans, else neither; not DD if the
# gene isn't diploid-significant in that cell type.
cis_cutoff   <- 0.70
trans_cutoff <- 0.30
for (ct in ct_common) {
  lfsr_col    <- paste0("dd.lfsr.diploid.", ct)
  cisprop_col <- paste0("dd.cis.proportion.", ct)
  out_col     <- paste0("dd.cis.category.", ct)

  lfsr    <- master_df[[lfsr_col]]
  cisprop <- master_df[[cisprop_col]]

  cat_vec <- rep("not DD", nrow(master_df))
  is_dd   <- !is.na(lfsr) & lfsr < 0.05 & !is.na(cisprop)

  cat_vec[is_dd & cisprop > cis_cutoff]                             <- "cis"
  cat_vec[is_dd & cisprop < trans_cutoff]                           <- "trans"
  cat_vec[is_dd & cisprop <= cis_cutoff & cisprop >= trans_cutoff]  <- "neither"

  master_df[[out_col]] <- factor(cat_vec, levels = c("cis", "neither", "trans", "not DD"))
}

# Overall cis/trans category across cell types: cis/trans if only that
# category appears (alongside neither/not DD), mixed if both cis and trans
# appear, else neither or not DD.
cat_cols <- paste0("dd.cis.category.", ct_common)
master_df <- master_df %>%
  rowwise() %>%
  mutate(
    dd.cis.category.overall = {
      vals <- as.character(c_across(all_of(cat_cols)))
      vals <- vals[!is.na(vals)]

      has_cis     <- "cis"     %in% vals
      has_trans   <- "trans"   %in% vals
      has_neither <- "neither" %in% vals

      if (!has_cis && !has_trans && !has_neither) "not DD"
      else if (has_cis && has_trans)              "mixed"
      else if (has_cis)                           "cis"
      else if (has_trans)                         "trans"
      else if (has_neither)                       "neither"
      else                                        "not DD"
    }
  ) %>%
  ungroup()

# Overall cis/trans category based on the mean cis proportion across
# diploid-DD cell types, rather than a per-celltype vote.
master_df <- master_df %>%
  mutate(
    dd.cis.category.overall.mean.definition = case_when(
      is.nan(dd.in.diploid.mean.cis.proportion) ~ "not DD",
      dd.in.diploid.mean.cis.proportion > 0.7    ~ "cis",
      dd.in.diploid.mean.cis.proportion < 0.3    ~ "trans",
      TRUE                                       ~ "neither"
    )
  )

write.csv(master_df, out_master, row.names = FALSE)
message(sprintf("Master table: %d genes × %d columns", nrow(master_df), ncol(master_df)))
