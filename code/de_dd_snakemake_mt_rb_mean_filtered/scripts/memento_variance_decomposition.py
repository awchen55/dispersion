"""
Estimate, per gene x cell type x species x donor (diploid) or gene x cell
type x species (tetraploid, no donor level), memento's raw 2nd-moment
variance ("total.variability") alongside its mean-corrected dispersion
(residual variance) estimate, plus -- per group -- how much of the
across-gene variance is explained by the mean-variance trend itself.

Why this script exists: memento.get_1d_moments() only returns log(mean) and
log(residual variance) -- the raw 2nd-moment variance it computes internally
(adata.uns['memento']['1d_moments'][group][1]) is discarded and never
written to any of this project's existing *_var_*.csv files (those files are
actually log(residual variance), i.e. what this project calls "dispersion"
throughout dispersion_diagnostics.Rmd / simulation_dispersion_diagnostics.Rmd
-- see memento/estimator.py's _fit_mv_regressor/_residual_variance and
memento/main.py's compute_1d_moments/get_1d_moments). This script re-runs
memento's own moment computation (identical params to the production
mt_rb_mean_filtered pipeline) and pulls the raw variance out of
adata.uns['memento'] directly, before it gets discarded.

NOTE on why there is no per-gene "percent of variance explained by
dispersion" column: memento's dispersion (residual variance) is defined as
raw_variance / exp(fitted_trend(log(mean))) -- a ratio RELATIVE to a fitted
mean-variance trend curve, deliberately renormalized to hover near 1
regardless of a gene's absolute expression level. It is not a variance
component in the same (additive) units as the raw 2nd moment, so
dispersion/variance does not reduce to a meaningful 0-1 fraction (verified:
it comes out anywhere from ~1e5 to ~1e9 across genes here, tracking almost
exactly 1/exp(fitted_trend), i.e. essentially the inverse mean-variance
trend curve, not "percent of variability"). What IS well-defined at the
group level is the R^2 of that same mean-variance trend fit itself --
report that instead, in a separate per-group summary file.

Outputs:
  --out-csv-genes  (long/tidy, one row per gene x sample):
    gene, celltype, ploidy, species, donor
    mean               -- log(mean), matches existing project convention
    dispersion.linear  -- residual variance (memento's mean-corrected dispersion),
                          LINEAR scale (existing project files store log() of this)
    total.variability  -- raw 2nd-moment variance, LINEAR scale

  --out-csv-summary  (one row per sample / group):
    celltype, ploidy, species, donor, n_genes
    mean.variance.trend.r_squared  -- R^2 of log(variance) ~ poly(log(mean), degree=2),
                                       i.e. how much of the across-gene variance
                                       differences are explained by mean expression
                                       level alone (matches memento's own
                                       _fit_mv_regressor exactly); the complement
                                       (1 - R^2) is the share left for dispersion
                                       to capture.

Usage:
  python memento_variance_decomposition.py \
      --h5ad <path> --celltype <name> --ploidy diploid|tetraploid \
      --q 0.2 --min-perc-group 0.95 --mean-expr-cutoff -6 \
      --min-cells-per-group 100 \
      --out-csv-genes <path> --out-csv-summary <path>
"""
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import memento
from scipy.sparse import issparse, csr_matrix

parser = argparse.ArgumentParser()
parser.add_argument("--h5ad", required=True)
parser.add_argument("--celltype", required=True)
parser.add_argument("--ploidy", required=True, choices=["diploid", "tetraploid"])
parser.add_argument("--q", type=float, default=0.2)
parser.add_argument("--min-perc-group", type=float, default=0.95)
parser.add_argument("--mean-expr-cutoff", type=float, default=-6.0)
parser.add_argument("--min-cells-per-group", type=int, default=100,
                     help="Diploid only -- species x donor_id groups smaller than this are dropped")
parser.add_argument("--out-csv-genes", required=True)
parser.add_argument("--out-csv-summary", required=True)
args = parser.parse_args()

adata = sc.read(args.h5ad)
if not (issparse(adata.X) and adata.X.format == "csr"):
    adata.X = csr_matrix(adata.X)
adata.var_names_make_unique()

if args.ploidy == "diploid":
    label_columns = ["species", "donor_id"]
    # Same low-cell-count group exclusion as memento_diploid.py, to keep this
    # consistent with the production dispersion/mean estimates.
    group_counts = adata.obs.groupby(["species", "donor_id"], observed=True).size()
    small_groups = group_counts[group_counts < args.min_cells_per_group]
    if len(small_groups) > 0:
        drop_idx = adata.obs.set_index(["species", "donor_id"]).index.isin(small_groups.index)
        adata = adata[~drop_idx].copy()
else:
    label_columns = ["species", "celltype"]

adata.obs["capture_rate"] = args.q
memento.setup_memento(adata, q_column="capture_rate")
memento.create_groups(adata, label_columns=label_columns)
memento.compute_1d_moments(adata, min_perc_group=args.min_perc_group)

genes = adata.var.index.to_numpy()
groups = [g for g in adata.uns["memento"]["groups"] if g in adata.uns["memento"]["1d_moments"]]

# Build log-mean matrix across groups first, to replicate the production
# mito/ribo + mean-expression-cutoff gene filter (memento_diploid.py /
# memento_tetraploid_dd_ct.py) before assembling the long output.
log_mean_df = pd.DataFrame({g: np.log(adata.uns["memento"]["1d_moments"][g][0]) for g in groups}, index=genes)

is_mt = pd.Series(genes, index=genes).str.match(r"^MT-")
is_ribo = pd.Series(genes, index=genes).str.match(r"^RP[LS]") & ~pd.Series(genes, index=genes).str.startswith("RPS6K")
is_high_expr = (log_mean_df > args.mean_expr_cutoff).any(axis=1)
retained_genes = genes[~(is_mt.to_numpy() | is_ribo.to_numpy() | is_high_expr.to_numpy())]
retained_mask = np.isin(genes, retained_genes)

gene_rows = []
summary_rows = []
for g in groups:
    parts = g.split("^")[1:]  # drop leading "sg"
    species = parts[0]
    donor = parts[1] if args.ploidy == "diploid" else None

    mean_arr, var_arr, res_var_arr = adata.uns["memento"]["1d_moments"][g][:3]
    mean_arr, var_arr, res_var_arr = mean_arr[retained_mask], var_arr[retained_mask], res_var_arr[retained_mask]

    gene_rows.append(pd.DataFrame({
        "gene": retained_genes,
        "celltype": args.celltype,
        "ploidy": args.ploidy,
        "species": species,
        "donor": donor,
        "mean": np.log(mean_arr),
        "dispersion.linear": res_var_arr,
        "total.variability": var_arr,
    }))

    # R^2 of memento's own mean-variance trend fit (log(var) ~ poly(log(mean), 2)),
    # i.e. _fit_mv_regressor in memento/estimator.py, reproduced here exactly.
    cond = (mean_arr > 0) & (var_arr > 0)
    m, v = np.log(mean_arr[cond]), np.log(var_arr[cond])
    if cond.sum() >= 3:
        poly = np.polyfit(m, v, 2)
        v_pred = np.poly1d(poly)(m)
        ss_res = np.sum((v - v_pred) ** 2)
        ss_tot = np.sum((v - v.mean()) ** 2)
        r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    else:
        r_squared = np.nan

    summary_rows.append({
        "celltype": args.celltype,
        "ploidy": args.ploidy,
        "species": species,
        "donor": donor,
        "n_genes": int(cond.sum()),
        "mean.variance.trend.r_squared": r_squared,
    })

gene_df = pd.concat(gene_rows, ignore_index=True)
gene_df.to_csv(args.out_csv_genes, index=False)

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(args.out_csv_summary, index=False)

print(f"Wrote {len(gene_df)} gene rows ({len(retained_genes)} genes x {len(groups)} groups) to {args.out_csv_genes}")
print(f"Wrote {len(summary_df)} summary rows to {args.out_csv_summary}")
