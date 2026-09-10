library(Seurat)
library(matrixTests)

seurat_path <- snakemake@input[["seurat"]]
out_path    <- snakemake@output[["csv"]]

tet <- readRDS(seurat_path)

human_cells <- rownames(tet@meta.data[tet@meta.data$species == "human", ])
chimp_cells <- rownames(tet@meta.data[tet@meta.data$species == "chimp", ])

rc_human <- GetAssayData(tet, assay = "RNA", slot = "counts")[, human_cells, drop = FALSE]
rc_chimp <- GetAssayData(tet, assay = "RNA", slot = "counts")[, chimp_cells, drop = FALSE]

DE <- row_wilcoxon_paired(
  as.matrix(rc_chimp), as.matrix(rc_human),
  alternative = "two.sided", null = 0, exact = NA, correct = TRUE
)

hmean      <- rowSums(rc_human)
cmean      <- rowSums(rc_chimp)
DE$logFC   <- log2((hmean + 1) / (cmean + 1))
DE$umi     <- hmean + cmean

write.csv(DE, out_path)
