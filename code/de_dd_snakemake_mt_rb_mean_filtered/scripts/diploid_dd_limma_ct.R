library(limma)
library(variancePartition)
library(dreamlet)

var_csv  <- snakemake@input[["var_csv"]]
out_path <- snakemake@output[["fit"]]
voomspan <- snakemake@params[["voomspan"]]
ddf      <- snakemake@params[["ddf"]]

memento_var <- read.csv(var_csv, row.names = 1, check.names = FALSE)

# parse metadata from column names: sg^species^donor_id^pool
cols       <- colnames(memento_var)
split_cols <- strsplit(cols, "\\^")
meta_df    <- do.call(rbind, lapply(split_cols, function(x) {
  data.frame(species    = x[2],
             individual = x[3],
             stringsAsFactors = FALSE)
}))
rownames(meta_df) <- cols

male_ids <- c("NA19160", "NA19210", "NA18913", "NA18519", "C3649", "C8861", "C4955")
meta_df$sex <- ifelse(meta_df$individual %in% male_ids, "M", "F")

norm_var <- limma::normalizeQuantiles(memento_var)

form <- ~ 0 + species
L    <- makeContrastsDream(formula = form, data = meta_df,
                           contrasts = c("specieshuman-specieschimp"))
fit  <- dream(norm_var, form, meta_df, L)
fit  <- eBayes(fit)

saveRDS(fit, out_path)
