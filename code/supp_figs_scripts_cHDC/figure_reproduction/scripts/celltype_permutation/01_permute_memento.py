#!/usr/bin/env python

import gc
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import memento
from scipy.sparse import csr_matrix

REPO_DIR = "/project/gilad/brendan/dispersion/pilot/cHDC_data"
# ANALYSIS_OK[file-selection]: one pinned file (rule flags the "old" in "threshold"); fingerprinted by total cell count below.
H5AD_PATH = os.path.join(
    REPO_DIR, "figure_reproduction", "data",
    "q_threshold_source_excl_epithelial_cluster24.h5ad")
OUT_DIR = os.path.join(
    REPO_DIR, "figure_reproduction", "outputs", "celltype_permutation", "memento")

LABEL_COL = "new.predicted.id"
IND_COL = "vireo.individual"
PERM_COL = "permuted_celltype"
CELL_TYPES = {
    "Cardiomyocytes", "Mesenchymal_cells", "Cardiac_fibroblasts",
    "Epithelial_cells", "Activated_fibroblasts", "Foregut_endodermal_cells",
}
INDIVIDUALS = {"NA19099", "NA19093", "NA18522"}
CAPTURE_RATE = 0.25    # matches memento_pip_forloop.py
MIN_PERC_GROUP = 0.95  # matches memento_pip_forloop.py
# ANALYSIS_OK[seed]: fixed seed for the single within-individual permutation; the original was unseeded (documented in PORT_CHANGELOG.md).
PERM_SEED = 123
EXPECTED_NCOLS = len(INDIVIDUALS) * len(CELL_TYPES)  # 18
EXPECTED_TOTAL_CELLS = 132784  # n_obs of H5AD_PATH as built by qthresh_build_source_h5ad

os.makedirs(OUT_DIR, exist_ok=True)

print(f"loading {H5AD_PATH}")
adata = sc.read(H5AD_PATH)
if adata.n_obs != EXPECTED_TOTAL_CELLS:
    sys.exit(f"{H5AD_PATH} has {adata.n_obs} cells, expected {EXPECTED_TOTAL_CELLS} -- wrong/rebuilt file")
# ANALYSIS_OK[layer-choice]: .X holds raw UMI counts (same h5ad and same check as q_threshold_sensitivity/01).
if type(adata.X) != csr_matrix:
    sys.exit("adata.X is not a csr_matrix -- raw-count assumption violated")

# ANALYSIS_OK[filtering]: the 6 manuscript cell types, identical to memento_pip_forloop.py; cell count logged below.
adata_subset = adata[adata.obs[LABEL_COL].isin(CELL_TYPES)].copy()
obs_cts = set(adata_subset.obs[LABEL_COL].unique())
if obs_cts != CELL_TYPES:
    sys.exit(f"Subset cell types {obs_cts} != expected {CELL_TYPES}")
obs_inds = set(adata_subset.obs[IND_COL].unique())
if obs_inds != INDIVIDUALS:
    sys.exit(f"Subset individuals {obs_inds} != expected {INDIVIDUALS}")
print(f"subset: {adata_subset.n_obs} of {adata.n_obs} cells x {adata_subset.n_vars} genes")

# Drop the duplicate .raw matrix before Memento (q_threshold_sensitivity/01's OOM fix).
adata_subset.raw = None
del adata
gc.collect()

# ------------------------- within-individual permutation ----------------------
rng = np.random.default_rng(PERM_SEED)
true_labels = adata_subset.obs[LABEL_COL].astype(str).to_numpy()
inds = adata_subset.obs[IND_COL].astype(str).to_numpy()
perm_labels = true_labels.copy()
for ind in sorted(INDIVIDUALS):
    idx = np.flatnonzero(inds == ind)
    perm_labels[idx] = rng.permutation(true_labels[idx])
adata_subset.obs[PERM_COL] = pd.Categorical(perm_labels)

true_counts = pd.crosstab(inds, true_labels)
perm_counts = pd.crosstab(inds, perm_labels)
if not true_counts.equals(perm_counts):
    sys.exit("Per-individual cell-type counts changed under permutation -- bug")
frac_same = float(np.mean(true_labels == perm_labels))
print(f"permutation (seed {PERM_SEED}): {frac_same:.4f} of cells keep their true label")
counts_out = perm_counts.stack().rename("n_cells").reset_index()
counts_out.columns = ["individual", "cell_type", "n_cells"]
counts_out["frac_cells_keeping_true_label"] = frac_same
counts_out.to_csv(os.path.join(OUT_DIR, "permutation_cell_counts.csv"), index=False)

# Composition of each permuted group by true label (should mirror each individual's mix).
comp = pd.crosstab([inds, perm_labels], true_labels, normalize="index")
print("true-label composition of each permuted group (rows sum to 1):")
print(comp.round(3).to_string())

# ------------------------------- Memento --------------------------------------
adata_subset.obs["capture_rate"] = CAPTURE_RATE
memento.setup_memento(adata_subset, q_column="capture_rate")
memento.create_groups(adata_subset, label_columns=[IND_COL, PERM_COL])
memento.compute_1d_moments(adata_subset, min_perc_group=MIN_PERC_GROUP)
mean_df, var_df, counts_df = memento.get_1d_moments(adata_subset)


def reformat(df, label):
    """sg^{individual}^{cell_type} -> {individual}_{cell_type}; drop genes with
    an NA in any group. Copied from q_threshold_sensitivity/01."""
    n_before = df.shape[0]
    # ANALYSIS_OK[feature-filter]: drops genes with an NA moment in any of the 18 groups (Memento's per-group filter, same as the published reformat); count logged.
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

mean_fmt.to_csv(os.path.join(OUT_DIR, "mean_expression_permuted.csv"))
var_fmt.to_csv(os.path.join(OUT_DIR, "res_var_permuted.csv"))
print(f"wrote {var_fmt.shape[0]} genes x {var_fmt.shape[1]} cols to {OUT_DIR}")
