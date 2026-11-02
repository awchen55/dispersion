#!/usr/bin/env python

import os
import numpy as np
import scanpy as sc
import memento
from scipy.sparse import csr_matrix, issparse

# ---- Named constants --------------------------------------------------------
REPO_DIR = "/project/gilad/brendan/dispersion/pilot/cHDC_data"
ANALYSIS_DIR = os.path.join(
    REPO_DIR, "analysis/epithelial_subcluster24_excl_sensitivity"
)
OUT_DIR = os.path.join(ANALYSIS_DIR, "outputs")

INPUT_H5AD = os.path.join(OUT_DIR, "epithelial_excl_cluster24_raw_counts.h5ad")
MEAN_OUT_CSV = os.path.join(OUT_DIR, "mean_expression_epithelial_excl_cluster24.csv")
VAR_OUT_CSV = os.path.join(OUT_DIR, "res_var_epithelial_excl_cluster24.csv")

INDIVIDUAL_COL = "vireo.individual"
CELL_TYPE_SUFFIX = "Epithelial_cells"
INDIVIDUALS = ["NA19099", "NA19093", "NA18522"]

# Canonical manuscript Memento parameters (memento_pip_forloop.py).
CAPTURE_RATE = 0.25
MIN_PERC_GROUP = 0.95


def log_msg(msg):
    print(f"[02_memento_epi_excl] {msg}", flush=True)


def assert_raw_counts(adata):
    """Fail loudly unless adata.X is an integer-valued, non-negative, sparse matrix."""
    x = adata.X
    if not issparse(x):
        raise ValueError("adata.X is dense; expected a sparse raw-counts matrix.")
    if not isinstance(x, csr_matrix):
        log_msg(f"adata.X is {type(x).__name__}; converting to CSR.")
        x = csr_matrix(x)
        # ANALYSIS_OK[layer-choice]: writes the CSR-converted raw counts back to
        # .X (same raw-count layer just read above); no layer semantics change.
        adata.X = x
    data = x.data
    if data.size == 0:
        raise ValueError("adata.X has no non-zero entries -- not raw counts.")
    if data.min() < 0:
        raise ValueError("adata.X contains negative values -- not raw counts.")
    max_int_err = float(np.max(np.abs(data - np.round(data))))
    if max_int_err > 1e-8:
        raise ValueError(
            f"adata.X is NOT integer-valued (max deviation {max_int_err:.3g}). "
            "Aborting to avoid feeding non-raw counts to Memento."
        )
    log_msg(
        f"Raw-count assertion passed: CSR, non-negative, integer "
        f"(max deviation {max_int_err:.3g}); {adata.n_obs} cells x {adata.n_vars} genes."
    )


def gene_indexed(df):
    """Normalize a Memento moment dataframe to a gene-indexed frame (defensive:
    Memento's get_1d_moments output has been observed with either the gene id
    already as the index, or as a separate 'gene' column -- see
    subcluster_dispersion_robustness/scripts/03_robustness_correlations.R's
    identical defensive handling)."""
    if "gene" in df.columns:
        df = df.set_index("gene")
    return df


def main():
    if not os.path.exists(INPUT_H5AD):
        raise FileNotFoundError(
            f"Missing input h5ad: {INPUT_H5AD}. Run "
            "01_build_excl_cluster24_epithelial_h5ad.py first."
        )
    os.makedirs(OUT_DIR, exist_ok=True)

    log_msg(f"Reading {INPUT_H5AD}")
    adata = sc.read(INPUT_H5AD)
    log_msg(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes.")
    assert_raw_counts(adata)

    if INDIVIDUAL_COL not in adata.obs.columns:
        raise KeyError(f"Missing obs column: {INDIVIDUAL_COL}")
    present = set(adata.obs[INDIVIDUAL_COL].astype(str))
    if not set(INDIVIDUALS).issubset(present):
        raise ValueError(f"Expected individuals {INDIVIDUALS}; found {sorted(present)}")

    n_genes_in = adata.n_vars

    adata.obs["capture_rate"] = CAPTURE_RATE
    memento.setup_memento(adata, q_column="capture_rate")
    memento.create_groups(adata, label_columns=[INDIVIDUAL_COL])
    memento.compute_1d_moments(adata, min_perc_group=MIN_PERC_GROUP)

    mean_df, var_df, counts_df = memento.get_1d_moments(adata)
    mean_df = gene_indexed(mean_df)
    var_df = gene_indexed(var_df)
    log_msg(
        f"Memento moments: mean {mean_df.shape}, var {var_df.shape} (genes x groups); "
        f"{n_genes_in} genes in -> {mean_df.shape[0]} genes with a computed estimate "
        f"({n_genes_in - mean_df.shape[0]} dropped by min_perc_group={MIN_PERC_GROUP})."
    )

    # Memento's group columns are Individual-only labels here (single label
    # column). Rename to the published matrix's "{IND}_{CellType}" convention.
    def rename_cols(df):
        new_cols = {}
        for c in df.columns:
            c_str = str(c)
            matches = [ind for ind in INDIVIDUALS if ind in c_str]
            if len(matches) != 1:
                raise ValueError(
                    f"Could not uniquely map Memento group column '{c}' to one of "
                    f"{INDIVIDUALS}. Aborting rather than guessing."
                )
            new_cols[c] = f"{matches[0]}_{CELL_TYPE_SUFFIX}"
        return df.rename(columns=new_cols)

    mean_df = rename_cols(mean_df)
    var_df = rename_cols(var_df)

    expected_cols = sorted(f"{ind}_{CELL_TYPE_SUFFIX}" for ind in INDIVIDUALS)
    if sorted(mean_df.columns) != expected_cols or sorted(var_df.columns) != expected_cols:
        raise ValueError(
            f"Expected exactly columns {expected_cols} after renaming; got "
            f"mean={sorted(mean_df.columns)}, var={sorted(var_df.columns)}. "
            "An individual was likely dropped entirely by min_perc_group filtering."
        )

    mean_df.to_csv(MEAN_OUT_CSV, index=True)
    var_df.to_csv(VAR_OUT_CSV, index=True)
    log_msg(f"Wrote {MEAN_OUT_CSV}")
    log_msg(f"Wrote {VAR_OUT_CSV}")
    log_msg("DONE. Next: scripts/03_build_modified_matrices.R")


if __name__ == "__main__":
    main()
