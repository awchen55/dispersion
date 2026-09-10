import scanpy as sc
import memento
from scipy.sparse import issparse, csr_matrix

h5ad_path        = snakemake.input["h5ad"]
out_var          = snakemake.output["var"]
out_mean         = snakemake.output["mean"]
q                   = float(snakemake.params["q"])
min_perc_group      = snakemake.params["min_perc_group"]
mean_expr_cutoff    = float(snakemake.params["mean_expr_cutoff"])
min_cells_per_group = int(snakemake.params["min_cells_per_group"])

adata = sc.read(h5ad_path)
if not (issparse(adata.X) and adata.X.format == "csr"):
    adata.X = csr_matrix(adata.X)

# Drop species x donor_id groups with too few cells for this cell type: with
# very few cells, most genes have 0-1 nonzero counts in that group, and
# memento's moment-based dispersion estimator collapses to near-identical,
# artificially very-negative "residual variance" values across hundreds of
# unrelated genes at once (a numerical artifact of the sparse-count regime,
# not a biological signal) -- see the donor cell-count investigation in
# q_threshold_sensitivity.Rmd (Section 8), which traced this exact pattern in
# Glioblast/ENS_glia/Retinal_neurons back to individual chimp donors with as
# few as 4-53 cells in that cell type, vs. hundreds-to-thousands for others.
group_counts = adata.obs.groupby(["species", "donor_id"], observed=True).size()
small_groups = group_counts[group_counts < min_cells_per_group]
if len(small_groups) > 0:
    drop_idx = adata.obs.set_index(["species", "donor_id"]).index.isin(small_groups.index)
    adata = adata[~drop_idx].copy()

adata.obs["capture_rate"] = q
memento.setup_memento(adata, q_column="capture_rate")
# Group by species × donor_id only — no pool.
memento.create_groups(adata, label_columns=["species", "donor_id"])
memento.compute_1d_moments(adata, min_perc_group=min_perc_group)

mean, var, _ = memento.get_1d_moments(adata)

if "gene" in var.columns:
    var = var.set_index("gene")
if "gene" in mean.columns:
    mean = mean.set_index("gene")

# Exclude mitochondrial and ribosomal protein genes: mitochondrial genes
# (MT-*) and cytoplasmic ribosomal protein genes (RPL*/RPS*) show
# translation- and metabolism-driven expression noise that is not of
# interest for interspecies DD, so they are dropped up front. The RPS6K*
# kinase family (RPS6KA1-6, RPS6KB1-2, RPS6KC1, RPS6KL1, and their antisense
# transcripts) is excluded from the ribosomal pattern since those genes are
# S6 kinases, not ribosomal proteins, despite the RPS-prefixed name.
is_mt   = mean.index.str.match(r"^MT-")
is_ribo = mean.index.str.match(r"^RP[LS]") & ~mean.index.str.startswith("RPS6K")

# On top of the MT/ribosomal drop, also apply the per-sample log-mean-expression
# cutoff: drop a gene if its log mean expression exceeds the cutoff in ANY
# species x donor group (not just on average), since memento's mean-variance
# correction is unreliable at the high-expression tail (see
# dispersion_diagnostics.Rmd, Sections 6-7).
is_high_expr = (mean > mean_expr_cutoff).any(axis=1)

retained_genes = mean.index[~(is_mt | is_ribo | is_high_expr)]
var  = var.loc[var.index.intersection(retained_genes)]
mean = mean.loc[retained_genes]

var.to_csv(out_var)
mean.to_csv(out_mean)
