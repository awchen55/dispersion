"""
Rerun the production memento DD test (binary_test_1d) for one allotetraploid
cell type after downsampling to a fixed number of physical (hybrid) cells.

Mirrors de_dd_snakemake_mt_rb_mean_filtered/scripts/memento_tetraploid_dd_ct.py
exactly (same filters/params), except the input AnnData is first subset to a
random sample of physical cells before the memento pipeline runs.

Each physical hybrid cell contributes two rows in the source h5ad -- one
"human_<barcode>" row (human-allele counts) and one "chimp_<barcode>" row
(chimp-allele counts), sharing the same barcode suffix. Downsampling must
therefore select whole physical cells (both rows), not rows independently,
or the paired species design memento's binary_test_1d relies on breaks.
"""
import argparse

import numpy as np
import scanpy as sc
import memento
from scipy.sparse import issparse, csr_matrix

parser = argparse.ArgumentParser()
parser.add_argument("--h5ad", required=True)
parser.add_argument("--celltype", required=True)
parser.add_argument("--n-cells", type=int, required=True,
                     help="number of physical hybrid cells to keep (per species this yields the same N, since each physical cell has one human + one chimp row)")
parser.add_argument("--seed", type=int, required=True)
parser.add_argument("--q", type=float, default=0.2)
parser.add_argument("--min-perc-group", type=float, default=0.95)
parser.add_argument("--mean-expr-cutoff", type=float, default=-6.0)
parser.add_argument("--num-cpus", type=int, default=8)
parser.add_argument("--num-boot", type=int, default=5000)
parser.add_argument("--out-csv", required=True)
args = parser.parse_args()

adata = sc.read(args.h5ad)
if not (issparse(adata.X) and adata.X.format == "csr"):
    adata.X = csr_matrix(adata.X)

adata.var_names_make_unique()

# --- downsample to N physical cells, keeping both species rows per cell ---
obs_names = adata.obs_names.to_numpy()
species = adata.obs["species"].astype(str).to_numpy()
suffix = np.array([n.split("_", 1)[1].strip() for n in obs_names])

human_suffix = set(suffix[species == "human"])
chimp_suffix = set(suffix[species == "chimp"])
paired_suffix = np.array(sorted(human_suffix & chimp_suffix))
assert len(paired_suffix) == len(human_suffix) == len(chimp_suffix), (
    f"expected fully paired human/chimp barcodes for {args.celltype}, "
    f"got {len(human_suffix)} human / {len(chimp_suffix)} chimp / "
    f"{len(paired_suffix)} paired"
)
assert args.n_cells <= len(paired_suffix), (
    f"requested {args.n_cells} cells but only {len(paired_suffix)} physical "
    f"cells available for {args.celltype}"
)

rng = np.random.default_rng(args.seed)
keep_suffix = set(rng.choice(paired_suffix, size=args.n_cells, replace=False))
keep_mask = np.isin(suffix, list(keep_suffix))
adata = adata[keep_mask].copy()

print(
    f"[{args.celltype}] downsampled to {args.n_cells} physical cells "
    f"({adata.n_obs} rows, seed={args.seed})"
)

adata.obs["capture_rate"] = args.q
memento.setup_memento(adata, q_column="capture_rate")
memento.create_groups(adata, label_columns=["species", "celltype"])
memento.compute_1d_moments(adata, min_perc_group=args.min_perc_group)

mean, var, _ = memento.get_1d_moments(adata)

if "gene" in var.columns:
    var = var.set_index("gene")
if "gene" in mean.columns:
    mean = mean.set_index("gene")

is_mt = mean.index.str.match(r"^MT-")
is_ribo = mean.index.str.match(r"^RP[LS]") & ~mean.index.str.startswith("RPS6K")
is_high_expr = (mean > args.mean_expr_cutoff).any(axis=1)

retained_genes = mean.index[~(is_mt | is_ribo | is_high_expr)]
var = var.loc[var.index.intersection(retained_genes)]
mean = mean.loc[retained_genes]

adata = adata[:, adata.var_names.isin(retained_genes)].copy()

adata.obs["species"] = adata.obs["species"].apply(lambda x: 0 if x == "chimp" else 1)

result = memento.binary_test_1d(
    adata=adata,
    capture_rate=args.q,
    treatment_col="species",
    num_cpus=args.num_cpus,
    num_boot=args.num_boot,
)

result["de_pval_adj"] = memento.util._fdrcorrect(result["de_pval"])
result["dv_pval_adj"] = memento.util._fdrcorrect(result["dv_pval"])

result.to_csv(args.out_csv, index=False)
