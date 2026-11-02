#!/usr/bin/env python

import os
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix, issparse

# ---- Named constants --------------------------------------------------------
REPO_DIR = "/project/gilad/brendan/dispersion/pilot/cHDC_data"
PARENT_ANALYSIS_DIR = os.path.join(
    REPO_DIR, "analysis/subcluster_dispersion_robustness"
)
ANALYSIS_DIR = os.path.join(
    REPO_DIR, "analysis/epithelial_subcluster24_excl_sensitivity"
)
OUT_DIR = os.path.join(ANALYSIS_DIR, "outputs")

INPUT_H5AD = os.path.join(PARENT_ANALYSIS_DIR, "outputs/subcluster_raw_counts.h5ad")
OUTPUT_H5AD = os.path.join(OUT_DIR, "epithelial_excl_cluster24_raw_counts.h5ad")

CELL_TYPE_COL = "new.predicted.id"
INDIVIDUAL_COL = "vireo.individual"
SUBCLUSTER_COL = "subcluster_res_2"  # resolution 2.0, where cluster 24 was flagged

TARGET_CELL_TYPE = "Epithelial_cells"
EXCLUDED_SUBCLUSTER = "Epithelial_cells::24"

INDIVIDUALS = ["NA19099", "NA19093", "NA18522"]

# Manuscript's per-individual inclusion threshold (matches
# subcluster_dispersion_robustness / the canonical cell-type-level pipeline).
MIN_CELLS_PER_INDIVIDUAL = 1000


def log_msg(msg):
    print(f"[01_build_excl_h5ad] {msg}", flush=True)


def assert_raw_counts(adata):
    """Fail loudly unless adata.X is an integer-valued, non-negative, sparse matrix."""
    x = adata.X
    if not issparse(x):
        raise ValueError("adata.X is dense; expected a sparse raw-counts matrix.")
    if not isinstance(x, csr_matrix):
        log_msg(f"adata.X is {type(x).__name__}; converting to CSR.")
        x = csr_matrix(x)
        # ANALYSIS_OK[layer-choice]: writes the CSR-converted raw counts back to
        # .X (same raw-count layer just read/subset above); no layer semantics
        # change, matches subcluster_dispersion_robustness's own assertion.
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


def main():
    if not os.path.exists(INPUT_H5AD):
        raise FileNotFoundError(
            f"Missing input h5ad: {INPUT_H5AD}. "
            "subcluster_dispersion_robustness/scripts/01_subcluster_and_export.R "
            "must have already been run (it was, 2026-07-30, job 52832619)."
        )
    os.makedirs(OUT_DIR, exist_ok=True)

    log_msg(f"Reading {INPUT_H5AD}")
    adata = sc.read(INPUT_H5AD)
    log_msg(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes.")
    assert_raw_counts(adata)

    for col in (CELL_TYPE_COL, INDIVIDUAL_COL, SUBCLUSTER_COL):
        if col not in adata.obs.columns:
            raise KeyError(
                f"Missing obs column '{col}'. Available: {list(adata.obs.columns)}"
            )

    present_individuals = set(adata.obs[INDIVIDUAL_COL].astype(str))
    if not set(INDIVIDUALS).issubset(present_individuals):
        raise ValueError(
            f"Expected individuals {INDIVIDUALS} not all present. "
            f"Found: {sorted(present_individuals)}"
        )

    epi_mask = adata.obs[CELL_TYPE_COL].astype(str) == TARGET_CELL_TYPE
    n_epi_total = int(epi_mask.sum())
    if n_epi_total == 0:
        raise ValueError(f"No cells found with {CELL_TYPE_COL} == '{TARGET_CELL_TYPE}'.")

    excl_mask = adata.obs[SUBCLUSTER_COL].astype(str) == EXCLUDED_SUBCLUSTER
    n_excl_total = int(excl_mask.sum())
    if n_excl_total == 0:
        raise ValueError(
            f"No cells found with {SUBCLUSTER_COL} == '{EXCLUDED_SUBCLUSTER}' -- "
            "the sub-cluster to exclude does not exist in this h5ad. Aborting "
            "rather than silently keeping all Epithelial cells."
        )

    keep_mask = epi_mask & (~excl_mask)
    n_keep_total = int(keep_mask.sum())
    log_msg(
        f"Epithelial_cells: {n_epi_total} total cells; "
        f"{EXCLUDED_SUBCLUSTER} accounts for {n_excl_total} of them; "
        f"{n_keep_total} cells retained after exclusion."
    )

    # Log + assert per-individual counts before and after exclusion against the
    # manuscript's >=1000-cells/individual inclusion threshold. Fail loudly if
    # any individual would drop below it -- silently proceeding with an
    # under-powered individual would invalidate the downstream Memento estimate.
    obs_epi = adata.obs.loc[epi_mask]
    obs_kept = adata.obs.loc[keep_mask]
    below_threshold = []
    for ind in INDIVIDUALS:
        n_before = int((obs_epi[INDIVIDUAL_COL].astype(str) == ind).sum())
        n_after = int((obs_kept[INDIVIDUAL_COL].astype(str) == ind).sum())
        log_msg(f"  {ind}: {n_before} Epithelial cells -> {n_after} after exclusion.")
        if n_after < MIN_CELLS_PER_INDIVIDUAL:
            below_threshold.append((ind, n_after))
    if below_threshold:
        raise ValueError(
            f"Individuals below the {MIN_CELLS_PER_INDIVIDUAL}-cell threshold after "
            f"excluding {EXCLUDED_SUBCLUSTER}: {below_threshold}. Aborting -- the "
            "downstream Epithelial Memento re-estimate would not be on the same "
            "footing as the published cell-type-level estimates."
        )

    adata_out = adata[keep_mask].copy()
    assert_raw_counts(adata_out)  # re-assert post-subset, matching repo convention

    # The input h5ad's SeuratDisk-produced .raw slot carries a var index
    # literally named "_index" (a Seurat->AnnData conversion artifact), which
    # collides with anndata's own reserved "_index" HDF5 dataset name on
    # write (ValueError: '_index' is a reserved name for dataframe columns).
    # We only need .X (raw counts) + .obs for Memento, so drop .raw rather
    # than work around the name collision.
    adata_out.raw = None

    if os.path.exists(OUTPUT_H5AD):
        os.remove(OUTPUT_H5AD)
    adata_out.write(OUTPUT_H5AD)
    log_msg(f"Wrote {OUTPUT_H5AD} ({adata_out.n_obs} cells x {adata_out.n_vars} genes).")
    log_msg("DONE. Next: scripts/02_memento_epithelial_excl_cluster24.py")


if __name__ == "__main__":
    main()
