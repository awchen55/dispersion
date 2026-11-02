#!/usr/bin/env python

import os
import scanpy as sc

REPO_DIR = "/project/gilad/brendan/dispersion/pilot/cHDC_data"
INPUT_H5AD = os.path.join(
    REPO_DIR, "analysis/epithelial_subcluster24_excl_sensitivity/outputs/",
    "epithelial_excl_cluster24_raw_counts.h5ad")
OUT_DIR = os.path.join(
    REPO_DIR, "figure_reproduction/outputs/coexpression_plasticity_lowfidelity")
OUT_TXT = os.path.join(OUT_DIR, "epithelial_excl_cluster24_barcodes.txt")

# Sanity check from epithelial_subcluster24_excl_sensitivity's own script 01
# log (2026-08-17 run): 26,914 cells retained after excluding sub-cluster 24
# (2,361 excluded from an original 29,275). Any deviation means the upstream
# h5ad changed and this number needs re-confirming, not silently trusting a
# new count.
EXPECTED_N_CELLS = 26914


def log_msg(msg):
    print(f"[00a_extract_barcodes] {msg}", flush=True)


def main():
    if not os.path.exists(INPUT_H5AD):
        raise FileNotFoundError(f"Missing input h5ad: {INPUT_H5AD}")
    os.makedirs(OUT_DIR, exist_ok=True)

    log_msg(f"Reading {INPUT_H5AD}")
    adata = sc.read(INPUT_H5AD)
    n_cells = adata.n_obs
    log_msg(f"Loaded {n_cells} cells x {adata.n_vars} genes.")
    if n_cells != EXPECTED_N_CELLS:
        raise ValueError(
            f"Expected {EXPECTED_N_CELLS} cells (per "
            "epithelial_subcluster24_excl_sensitivity's 2026-08-17 run log), "
            f"got {n_cells} -- the upstream h5ad may have changed; "
            "investigate before reusing this barcode list.")

    barcodes = list(adata.obs_names)
    if len(set(barcodes)) != len(barcodes):
        raise ValueError("Duplicate cell barcodes found -- aborting.")

    with open(OUT_TXT, "w") as f:
        f.write("\n".join(barcodes) + "\n")
    log_msg(f"Wrote {len(barcodes)} barcodes to {OUT_TXT}")


if __name__ == "__main__":
    main()
