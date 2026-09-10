import scanpy as sc
import memento
from scipy.sparse import issparse, csr_matrix

h5ad_path        = snakemake.input["h5ad"]
out_csv          = snakemake.output["csv"]
out_mean         = snakemake.output["mean"]
out_var          = snakemake.output["var"]
q                = float(snakemake.params["q"])
min_perc_group   = snakemake.params["min_perc_group"]
mean_expr_cutoff = float(snakemake.params["mean_expr_cutoff"])
num_cpus         = snakemake.params["num_cpus"]
num_boot         = snakemake.params["num_boot"]

adata = sc.read(h5ad_path)
if not (issparse(adata.X) and adata.X.format == "csr"):
    adata.X = csr_matrix(adata.X)

adata.var_names_make_unique()
adata.obs["capture_rate"] = q
memento.setup_memento(adata, q_column="capture_rate")
# Group by species x celltype only (celltype is constant per file, so this is
# effectively per-species pooled moments -- there is no donor/replicate level
# for the allotetraploid data, unlike the diploid samples).
memento.create_groups(adata, label_columns=["species", "celltype"])
memento.compute_1d_moments(adata, min_perc_group=min_perc_group)

mean, var, _ = memento.get_1d_moments(adata)

if "gene" in var.columns:
    var = var.set_index("gene")
if "gene" in mean.columns:
    mean = mean.set_index("gene")

# Same MT/ribosomal + mean-expression filter as the diploid pipeline
# (memento_diploid.py): drop mitochondrial (MT-*) and ribosomal protein
# (RPL*/RPS*, excluding the RPS6K* kinase family) genes, and drop genes with
# log mean expression above the cutoff in EITHER species group, before DD
# testing.
is_mt        = mean.index.str.match(r"^MT-")
is_ribo      = mean.index.str.match(r"^RP[LS]") & ~mean.index.str.startswith("RPS6K")
is_high_expr = (mean > mean_expr_cutoff).any(axis=1)

retained_genes = mean.index[~(is_mt | is_ribo | is_high_expr)]
var  = var.loc[var.index.intersection(retained_genes)]
mean = mean.loc[retained_genes]

var.to_csv(out_var)
mean.to_csv(out_mean)

# Restrict the DD test itself to the same retained gene set so the final
# tetraploid DD results are consistent with the mean/var filter above.
adata = adata[:, adata.var_names.isin(retained_genes)].copy()

# encode species: chimp=0, human=1
adata.obs["species"] = adata.obs["species"].apply(lambda x: 0 if x == "chimp" else 1)

result = memento.binary_test_1d(
    adata=adata,
    capture_rate=q,
    treatment_col="species",
    num_cpus=num_cpus,
    num_boot=num_boot,
)

result["de_pval_adj"] = memento.util._fdrcorrect(result["de_pval"])
result["dv_pval_adj"] = memento.util._fdrcorrect(result["dv_pval"])

result.to_csv(out_csv, index=False)
