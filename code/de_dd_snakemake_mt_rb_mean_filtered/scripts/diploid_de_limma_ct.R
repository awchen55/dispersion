library(edgeR)
library(variancePartition)
library(dreamlet)

pb_path  <- snakemake@input[["pb"]]
out_path <- snakemake@output[["fit"]]
voomspan  <- snakemake@params[["voomspan"]]
min_count <- snakemake@params[["min_count"]]
ddf       <- snakemake@params[["ddf"]]

bulk <- readRDS(pb_path)

keep <- rowSums(bulk$counts) > 0
bulk$counts <- bulk$counts[keep, , drop = FALSE]

d0   <- DGEList(bulk$counts)
keep <- filterByExpr(d0$counts, group = bulk$meta$species,
                     lib.size = d0$samples$lib.size, min.count = min_count)
if (any(!keep)) d0 <- d0[keep, ]
d0 <- calcNormFactors(d0, method = "TMMwsp")

form <- ~ 0 + species
vobj <- voomWithDreamWeights(d0, form, bulk$meta,
                              normalize.method = "cyclicloess",
                              plot = FALSE, save.plot = FALSE,
                              span = voomspan)
L   <- getContrast(vobj, form, bulk$meta, c("specieshuman", "specieschimp"))
fit <- dream(vobj, form, bulk$meta, L, ddf = ddf)
fit <- eBayes(fit)

saveRDS(fit, out_path)
