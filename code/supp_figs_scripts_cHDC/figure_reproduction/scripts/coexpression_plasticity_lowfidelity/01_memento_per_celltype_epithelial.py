#!/usr/bin/env python

import gc
import os

import numpy as np
import scanpy as sc
import memento
from scipy.sparse import csr_matrix, issparse

REPO_DIR = "/project/gilad/brendan/dispersion/pilot/cHDC_data"
H5AD_PATH = os.path.join(
    REPO_DIR, "analysis/epithelial_subcluster24_excl_sensitivity/outputs/",
    "epithelial_excl_cluster24_raw_counts.h5ad")
OUT_DIR = os.path.join(
    REPO_DIR, "figure_reproduction/outputs/coexpression_plasticity_lowfidelity",
    "memento_per_celltype")

RES_NAME = "new.predicted.id"
INDIVIDUAL_COL = "vireo.individual"
CELL_TYPE = "Epithelial_cells"
INDIVIDUALS = {"NA19099", "NA19093", "NA18522"}

CAPTURE_RATE = 0.25
MIN_PERC_GROUP = 0.9
EXPECTED_NCOLS = len(INDIVIDUALS)


def log_msg(msg):
    print(f"[01_memento_per_celltype_epithelial] {msg}", flush=True)


def assert_raw_counts(adata):
    x = adata.X
    if not issparse(x):
        raise ValueError("adata.X is dense; expected a sparse raw-counts matrix.")
    if not isinstance(x, csr_matrix):
        log_msg(f"adata.X is {type(x).__name__}; converting to CSR.")
        x = csr_matrix(x)
        adata.X = x
    data = x.data
    if data.size == 0:
        raise ValueError("adata.X has no non-zero entries -- not raw counts.")
    if data.min() < 0:
        raise ValueError("adata.X contains negative values -- not raw counts.")
    max_int_err = float(np.max(np.abs(data - np.round(data))))
    if max_int_err > 1e-8:
        raise ValueError(
            f"adata.X is NOT integer-valued (max deviation {max_int_err:.3g}) -- "
            "not raw counts. Aborting.")
    log_msg(
        f"Raw-count assertion passed: CSR, non-negative, integer "
        f"(max deviation {max_int_err:.3g}); {adata.n_obs} cells x {adata.n_vars} genes.")


def reformat(df, label):
    n_before = df.shape[0]
    df = df.dropna(how="any")
    log_msg(f"  [{label}] dropna: {n_before - df.shape[0]} of {n_before} genes "
            f"dropped (NA in >=1 of {EXPECTED_NCOLS} groups)")
    if "gene" in df.columns:
        df = df.set_index("gene")
    parts = [c.split("^") for c in df.columns]
    bad_cols = [c for c, p in zip(df.columns, parts) if len(p) != 2]
    if bad_cols:
        raise ValueError(
            f"[{label}] expected 'sg^individual' columns (single label_column), "
            f"got unparseable: {bad_cols}")
    individuals = [p[1] for p in parts]
    bad_ind = set(individuals) - INDIVIDUALS
    if bad_ind:
        raise ValueError(f"[{label}] unexpected individual token(s): {bad_ind}")
    df.columns = individuals
    return df


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    log_msg(f"Cell type: {CELL_TYPE} (MIN_PERC_GROUP={MIN_PERC_GROUP}, "
            f"CAPTURE_RATE={CAPTURE_RATE})")

    log_msg(f"Loading {H5AD_PATH}")
    adata = sc.read(H5AD_PATH)
    log_msg(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes.")

    if INDIVIDUAL_COL not in adata.obs.columns:
        raise KeyError(f"Missing obs column: {INDIVIDUAL_COL}")

    # This h5ad is already Epithelial_cells-only (sub-cluster 24 excluded, by
    # epithelial_subcluster24_excl_sensitivity's script 01). Assert that
    # explicitly instead of re-subsetting by cell type, so a future change to
    # the upstream h5ad's contents fails loudly rather than silently
    # including the wrong cells.
    if RES_NAME in adata.obs.columns:
        ct_vals = set(adata.obs[RES_NAME].astype(str).unique())
        if ct_vals != {CELL_TYPE}:
            raise ValueError(
                f"Expected this h5ad to contain ONLY {CELL_TYPE!r}, found: {ct_vals}. "
                "This script assumes no further cell-type subsetting is needed.")
    else:
        log_msg(f"  NOTE: '{RES_NAME}' column not present in this h5ad's obs -- "
                f"cannot double-check cell-type purity; proceeding on the assumption "
                f"(from epithelial_subcluster24_excl_sensitivity's own script 01) that "
                f"every cell here is {CELL_TYPE}.")

    adata_ct = adata.copy()
    del adata
    gc.collect()

    obs_ind = set(adata_ct.obs[INDIVIDUAL_COL].unique())
    unexpected_ind = obs_ind - INDIVIDUALS
    if unexpected_ind:
        raise ValueError(f"Unexpected individual(s) present: {unexpected_ind}")
    missing_ind = INDIVIDUALS - obs_ind
    if missing_ind:
        log_msg(f"  WARNING: individual(s) with zero cells: {missing_ind}")
    ind_counts = adata_ct.obs[INDIVIDUAL_COL].value_counts()
    for ind in sorted(INDIVIDUALS):
        log_msg(f"  {ind}: {int(ind_counts.get(ind, 0))} cells")

    assert_raw_counts(adata_ct)

    adata_ct.obs["capture_rate"] = CAPTURE_RATE
    memento.setup_memento(adata_ct, q_column="capture_rate")
    memento.create_groups(adata_ct, label_columns=[INDIVIDUAL_COL])
    memento.compute_1d_moments(adata_ct, min_perc_group=MIN_PERC_GROUP)

    mean_df, var_df, counts_df = memento.get_1d_moments(adata_ct)
    log_msg(f"  Memento moments (pre-reformat): mean {mean_df.shape}, var {var_df.shape} "
            f"(genes x groups).")

    mean_fmt = reformat(mean_df, "mean")
    var_fmt = reformat(var_df, "var")

    if mean_fmt.shape[1] != EXPECTED_NCOLS or var_fmt.shape[1] != EXPECTED_NCOLS:
        raise ValueError(
            f"Expected {EXPECTED_NCOLS} individual columns; got "
            f"mean={mean_fmt.shape[1]}, var={var_fmt.shape[1]}")

    mean_path = os.path.join(OUT_DIR, f"{CELL_TYPE}_mean_expression.csv")
    var_path = os.path.join(OUT_DIR, f"{CELL_TYPE}_res_var.csv")
    mean_fmt.to_csv(mean_path)
    var_fmt.to_csv(var_path)
    log_msg(f"[{CELL_TYPE}] wrote {mean_fmt.shape[0]} genes x {mean_fmt.shape[1]} "
            f"individuals to {mean_path} / {var_path}")


if __name__ == "__main__":
    main()
