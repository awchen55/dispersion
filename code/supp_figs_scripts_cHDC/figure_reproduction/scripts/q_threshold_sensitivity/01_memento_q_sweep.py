#!/usr/bin/env python

import gc
import os
import sys
import pandas as pd
import scanpy as sc
import memento
from scipy.sparse import csr_matrix

REPO_DIR = "/project/gilad/brendan/dispersion/pilot/cHDC_data"
# PATCHED for figure_reproduction (see PORT_CHANGELOG_qthreshold.md): reads
# the Epithelial_cells sub-cluster-24-excluded joint source h5ad (built by
# 00_build_source_h5ad_excl_epithelial_cluster24.py) instead of the original
# 6-cell-type raw-counts h5ad. Only Epithelial_cells' cells differ; every
# other cell type/individual is byte-identical to the original.
H5AD_PATH = os.path.join(
    REPO_DIR, "figure_reproduction", "data",
    "q_threshold_source_excl_epithelial_cluster24.h5ad")
OUT_ROOT = os.path.join(
    REPO_DIR, "figure_reproduction", "outputs", "q_threshold_sensitivity", "memento")

RES_NAME = "new.predicted.id"
CELL_TYPES = {
    "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
    "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells",
}
INDIVIDUALS = {"NA19099", "NA19093", "NA18522"}
MIN_PERC_GROUP = 0.95  # matches memento_pip_forloop.py
EXPECTED_NCOLS = len(INDIVIDUALS) * len(CELL_TYPES)  # 18

if "Q_VALUE" not in os.environ:
    sys.exit("Q_VALUE environment variable not set -- aborting (no default q)")
q_str = os.environ["Q_VALUE"]
Q_VALUE = float(q_str)
if not (0 < Q_VALUE <= 1):
    sys.exit(f"Q_VALUE={Q_VALUE} is not a valid capture rate in (0, 1]")

out_dir = os.path.join(OUT_ROOT, f"q{q_str}")
os.makedirs(out_dir, exist_ok=True)

print(f"[q={q_str}] loading {H5AD_PATH}")
adata = sc.read(H5AD_PATH)
# ANALYSIS_OK[layer-choice]: .X holds raw UMI counts, per this h5ad's own filename
# (100k_..._raw_counts_only.h5ad) and memento_pip_forloop.py's identical use of
# adata.X unmodified; checked immediately below, not assumed silently.
if type(adata.X) != csr_matrix:
    sys.exit("adata.X is not a csr_matrix -- raw-count assumption violated")

adata_subset = adata[adata.obs[RES_NAME].isin(CELL_TYPES)].copy()
obs_cts = set(adata_subset.obs[RES_NAME].unique())
if obs_cts != CELL_TYPES:
    sys.exit(f"Subset cell types {obs_cts} != expected {CELL_TYPES}")
print(f"[q={q_str}] subset: {adata_subset.n_obs} cells x {adata_subset.n_vars} genes")

# The source h5ad stores a full duplicate of the counts matrix under .raw (confirmed
# via h5ls: /X and /raw/X are both ~982M-nonzero CSR matrices) that memento_pip_forloop.py
# never reads. Dropping it -- and the now-unneeded full adata -- before the memory-heavy
# Memento moment computation avoids holding two redundant copies of the expression matrix
# at once (root cause of this script's first-run OOM at 100G; see .living/learnings.md).
adata_subset.raw = None
del adata
gc.collect()

adata_subset.obs["capture_rate"] = Q_VALUE

memento.setup_memento(adata_subset, q_column="capture_rate")
memento.create_groups(adata_subset, label_columns=["vireo.individual", RES_NAME])
memento.compute_1d_moments(adata_subset, min_perc_group=MIN_PERC_GROUP)

mean_df, var_df, counts_df = memento.get_1d_moments(adata_subset)


def reformat(df, label):
    """sg^{individual}^{cell_type} columns -> {individual}_{cell_type}; drop
    any gene row with an NA in any group (matches na.omit's row-wise
    behaviour, the repo's existing precedent for this reformat step)."""
    n_before = df.shape[0]
    # ANALYSIS_OK[feature-filter]: drops genes with an NA moment estimate in any
    # of the 18 individual x cell-type groups (Memento's internal per-group
    # filtering, which varies with capture_rate -- the exact effect under test
    # in this analysis); dropped-row count logged immediately below.
    df = df.dropna(how="any")
    print(f"[{label}] dropna: {n_before - df.shape[0]} of {n_before} genes dropped "
          f"(NA in >=1 of 18 groups)")
    if "gene" in df.columns:
        df = df.set_index("gene")
    parts = [c.split("^") for c in df.columns]
    bad_cols = [c for c, p in zip(df.columns, parts) if len(p) != 3]
    if bad_cols:
        raise ValueError(f"[{label}] expected 'sg^individual^cell_type' columns, "
                          f"got unparseable: {bad_cols}")
    individuals = [p[1] for p in parts]
    cell_types = [p[2] for p in parts]
    bad_ind = set(individuals) - INDIVIDUALS
    bad_ct = set(cell_types) - CELL_TYPES
    if bad_ind or bad_ct:
        raise ValueError(f"[{label}] unexpected individual/cell-type tokens: "
                          f"individuals={bad_ind}, cell_types={bad_ct}")
    df.columns = [f"{i}_{c}" for i, c in zip(individuals, cell_types)]
    return df


mean_fmt = reformat(mean_df, "mean")
var_fmt = reformat(var_df, "var")

if mean_fmt.shape[1] != EXPECTED_NCOLS or var_fmt.shape[1] != EXPECTED_NCOLS:
    raise ValueError(f"Expected {EXPECTED_NCOLS} columns; got "
                      f"mean={mean_fmt.shape[1]}, var={var_fmt.shape[1]}")

mean_fmt.to_csv(os.path.join(out_dir, "mean_expression.csv"))
var_fmt.to_csv(os.path.join(out_dir, "res_var.csv"))
print(f"[q={q_str}] wrote {mean_fmt.shape[0]} genes x {mean_fmt.shape[1]} cols to {out_dir}")
