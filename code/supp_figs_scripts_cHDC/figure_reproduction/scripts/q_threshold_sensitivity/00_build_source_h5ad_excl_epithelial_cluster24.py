#!/usr/bin/env python

import os
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix, issparse

REPO_DIR = "/project/gilad/brendan/dispersion/pilot/cHDC_data"

# Same source subcluster_dispersion_robustness/epithelial_subcluster24_excl_sensitivity
# already used to identify cluster 24's cells -- reused here only to read
# obs (cell barcodes + subcluster label), not to re-derive anything.
SUBCLUSTER_H5AD = os.path.join(
    REPO_DIR, "analysis", "subcluster_dispersion_robustness", "outputs",
    "subcluster_raw_counts.h5ad")
SUBCLUSTER_COL = "subcluster_res_2"
EXCLUDED_SUBCLUSTER = "Epithelial_cells::24"

# q_threshold_sensitivity's own joint raw-counts source (untouched apart from
# this cell-level filter).
TARGET_H5AD = os.path.join(
    REPO_DIR, "cellranger_cluster-mode_trial", "analysis", "datasets",
    "100k_permuted_annots_seurat_for_adata_conv_raw_counts_only.h5ad")

CELL_TYPE_COL = "new.predicted.id"
TARGET_CELL_TYPE = "Epithelial_cells"
INDIVIDUAL_COL = "vireo.individual"
INDIVIDUALS = ["NA19099", "NA19093", "NA18522"]
MIN_CELLS_PER_INDIVIDUAL = 1000  # manuscript's per-individual inclusion threshold

OUT_H5AD = os.path.join(
    REPO_DIR, "figure_reproduction", "data",
    "q_threshold_source_excl_epithelial_cluster24.h5ad")


def log_msg(msg):
    print(f"[00_build_source_h5ad_excl_cluster24] {msg}", flush=True)


def assert_raw_counts(adata):
    x = adata.X
    if not issparse(x):
        raise ValueError("adata.X is dense; expected a sparse raw-counts matrix.")
    if not isinstance(x, csr_matrix):
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
    log_msg(f"Raw-count assertion passed: {adata.n_obs} cells x {adata.n_vars} genes "
            f"(max int deviation {max_int_err:.3g}).")


def main():
    os.makedirs(os.path.dirname(OUT_H5AD), exist_ok=True)

    log_msg(f"Reading sub-cluster source: {SUBCLUSTER_H5AD}")
    sub_adata = sc.read(SUBCLUSTER_H5AD)
    if SUBCLUSTER_COL not in sub_adata.obs.columns:
        raise KeyError(f"Missing obs column '{SUBCLUSTER_COL}' in {SUBCLUSTER_H5AD}")
    excl_mask = sub_adata.obs[SUBCLUSTER_COL].astype(str) == EXCLUDED_SUBCLUSTER
    excluded_barcodes = set(sub_adata.obs_names[excl_mask])
    if len(excluded_barcodes) == 0:
        raise ValueError(
            f"No cells found with {SUBCLUSTER_COL} == '{EXCLUDED_SUBCLUSTER}' in "
            f"{SUBCLUSTER_H5AD} -- aborting rather than silently excluding nothing.")
    log_msg(f"Identified {len(excluded_barcodes)} cell(s) to exclude "
            f"({EXCLUDED_SUBCLUSTER}).")
    del sub_adata

    log_msg(f"Reading target (joint, 6-cell-type) source: {TARGET_H5AD}")
    adata = sc.read(TARGET_H5AD)
    log_msg(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes.")
    assert_raw_counts(adata)

    # Barcode-identity assertion -- fail loudly on any mismatch between the
    # two export pipelines rather than silently excluding a partial/wrong set.
    found_in_target = excluded_barcodes & set(adata.obs_names)
    missing = excluded_barcodes - found_in_target
    if len(missing) > 0:
        raise ValueError(
            f"{len(missing)} of {len(excluded_barcodes)} excluded barcodes from "
            f"{SUBCLUSTER_H5AD} were NOT found in {TARGET_H5AD} -- the two "
            "h5ad files' cell-barcode namespaces do not match as assumed. "
            "Aborting rather than proceeding with a partial/incorrect exclusion. "
            f"Example missing barcode(s): {sorted(missing)[:5]}")

    excl_mask_target = adata.obs_names.isin(excluded_barcodes)
    # Cross-check: every excluded barcode must actually be labeled
    # Epithelial_cells in the target object too.
    wrong_ct = adata.obs.loc[excl_mask_target, CELL_TYPE_COL].astype(str) != TARGET_CELL_TYPE
    if wrong_ct.any():
        raise ValueError(
            f"{int(wrong_ct.sum())} excluded barcode(s) are NOT labeled "
            f"'{TARGET_CELL_TYPE}' in the target h5ad -- barcode/label mismatch "
            "between the two sources. Aborting.")

    n_epi_total = int((adata.obs[CELL_TYPE_COL].astype(str) == TARGET_CELL_TYPE).sum())
    log_msg(f"Target Epithelial_cells total: {n_epi_total}; "
            f"excluding {int(excl_mask_target.sum())} (sub-cluster 24) -> "
            f"{n_epi_total - int(excl_mask_target.sum())} retained.")

    # Per-individual count check within Epithelial_cells, same threshold as
    # epithelial_subcluster24_excl_sensitivity's own script 01 (defense in
    # depth -- that analysis already validated this exact exclusion clears
    # the threshold; re-checked here since this is a different source object).
    epi_mask = adata.obs[CELL_TYPE_COL].astype(str) == TARGET_CELL_TYPE
    keep_epi_mask = epi_mask & (~excl_mask_target)
    below_threshold = []
    for ind in INDIVIDUALS:
        n_after = int((adata.obs.loc[keep_epi_mask, INDIVIDUAL_COL].astype(str) == ind).sum())
        log_msg(f"  Epithelial_cells / {ind}: {n_after} cells after exclusion.")
        if n_after < MIN_CELLS_PER_INDIVIDUAL:
            below_threshold.append((ind, n_after))
    if below_threshold:
        raise ValueError(
            f"Individual(s) below the {MIN_CELLS_PER_INDIVIDUAL}-cell Epithelial_cells "
            f"threshold after exclusion: {below_threshold}. Aborting.")

    keep_mask = ~excl_mask_target
    adata_out = adata[keep_mask].copy()
    del adata
    assert_raw_counts(adata_out)

    # Same SeuratDisk .raw reserved-name gotcha documented in
    # epithelial_subcluster24_excl_sensitivity's script 01 -- drop .raw before
    # writing (only .X + .obs are needed downstream).
    adata_out.raw = None

    if os.path.exists(OUT_H5AD):
        os.remove(OUT_H5AD)
    adata_out.write(OUT_H5AD)
    log_msg(f"Wrote {OUT_H5AD} ({adata_out.n_obs} cells x {adata_out.n_vars} genes).")


if __name__ == "__main__":
    main()
