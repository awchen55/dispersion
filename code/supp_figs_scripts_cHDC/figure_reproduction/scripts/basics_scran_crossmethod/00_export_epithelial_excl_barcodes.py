#!/usr/bin/env python

import os
import pandas as pd
import scanpy as sc

REPO_DIR = "/project/gilad/brendan/dispersion/pilot/cHDC_data"
INPUT_H5AD = os.path.join(
    REPO_DIR, "analysis/epithelial_subcluster24_excl_sensitivity/outputs",
    "epithelial_excl_cluster24_raw_counts.h5ad")
OUT_CSV = os.path.join(
    REPO_DIR, "figure_reproduction/outputs/basics_scran_crossmethod",
    "epithelial_excl_cluster24_barcodes.csv")

if not os.path.exists(INPUT_H5AD):
    raise FileNotFoundError(f"Missing input h5ad: {INPUT_H5AD}")

adata = sc.read(INPUT_H5AD)
print(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes from {INPUT_H5AD}")
if adata.n_obs == 0:
    raise ValueError("h5ad has zero cells -- aborting")

os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
pd.Series(adata.obs_names, name="barcode").to_csv(OUT_CSV, index=False)
print(f"Wrote {adata.n_obs} barcodes to {OUT_CSV}")
