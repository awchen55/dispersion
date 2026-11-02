#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import anndata as ad

REPO_DIR = "/project/gilad/brendan/dispersion/pilot/cHDC_data"
FULL_H5AD = os.path.join(
    REPO_DIR, "analysis/subcluster_dispersion_robustness/outputs/subcluster_raw_counts.h5ad")
EXCL_H5AD = os.path.join(
    REPO_DIR, "analysis/epithelial_subcluster24_excl_sensitivity/outputs/"
              "epithelial_excl_cluster24_raw_counts.h5ad")
ORIG_CSV = os.path.join(
    REPO_DIR, "cellranger_cluster-mode_trial/analysis/datasets/"
              "aggregated_metadata_by_celltype_cellline_new_pre_id.csv")
OUT_CSV = os.path.join(REPO_DIR, "figure_reproduction/data/aggregated_metadata_excl_cluster24.csv")

CELL_TYPE_COL = "new.predicted.id"
INDIVIDUAL_COL = "vireo.individual"
SUBCLUSTER_COL = "subcluster_res_2"
TARGET_CELL_TYPE = "Epithelial_cells"
EXCLUDED_SUBCLUSTER = "Epithelial_cells::24"
INDIVIDUALS = ["NA19099", "NA19093", "NA18522"]
RECOMPUTED_COLS = ["nCount_RNA", "nFeature_RNA", "percent.mt", "cell_num"]
REL_TOL = 1e-6  # validation tolerance for recomputed means vs original CSV


def per_cell_qc(adata):
    """Per-cell nCount_RNA / nFeature_RNA from X, checked against obs; percent.mt from ^MT- genes."""
    X = adata.X.tocsr()
    n_count = np.asarray(X.sum(axis=1)).ravel()
    n_feature = np.asarray((X > 0).sum(axis=1)).ravel()
    if not np.allclose(n_count, adata.obs["nCount_RNA"].to_numpy(), rtol=1e-9):
        raise ValueError("X row sums != obs nCount_RNA -- X is not the full RNA count matrix")
    if not np.array_equal(n_feature, adata.obs["nFeature_RNA"].to_numpy().astype(int)):
        raise ValueError("X nonzero counts != obs nFeature_RNA -- X is not the full RNA count matrix")
    mt_mask = adata.var_names.str.startswith("MT-")
    print(f"[qc] {mt_mask.sum()} MT- genes")
    if mt_mask.sum() != 13:
        raise ValueError(f"expected 13 MT- genes, found {mt_mask.sum()}")
    mt_count = np.asarray(X[:, mt_mask].sum(axis=1)).ravel()
    return pd.DataFrame({
        INDIVIDUAL_COL: adata.obs[INDIVIDUAL_COL].astype(str).to_numpy(),
        "nCount_RNA": n_count,
        "nFeature_RNA": n_feature,
        "percent.mt": 100 * mt_count / n_count,
    }, index=adata.obs_names)


def aggregate(qc):
    """Per-individual means (aggregateToPseudoBulk aggr_means) + cell_num."""
    agg = qc.groupby(INDIVIDUAL_COL)[["nCount_RNA", "nFeature_RNA", "percent.mt"]].mean()
    agg["cell_num"] = qc.groupby(INDIVIDUAL_COL).size()
    return agg.loc[INDIVIDUALS]


# ---- Load full Epithelial_cells population ----------------------------------
full = ad.read_h5ad(FULL_H5AD, backed="r")
epi_mask = (full.obs[CELL_TYPE_COL] == TARGET_CELL_TYPE).to_numpy()
print(f"[load] {full.n_obs} cells in full h5ad; {epi_mask.sum()} {TARGET_CELL_TYPE}")
epi = full[epi_mask].to_memory()
qc = per_cell_qc(epi)
subcl = epi.obs[SUBCLUSTER_COL].astype(str)

# ---- Validate method against the original CSV --------------------------------
orig = pd.read_csv(ORIG_CSV)
print(f"[load] original CSV: {orig.shape[0]} samples (all cell types, not just the 6 analysed)")
if orig["sample_id"].duplicated().any():
    raise ValueError("original CSV has duplicated sample_id values")
orig_epi = orig[orig[CELL_TYPE_COL] == TARGET_CELL_TYPE].set_index(INDIVIDUAL_COL).loc[INDIVIDUALS]
full_agg = aggregate(qc)
print("[validate] recomputed full-Epithelial covariates:\n", full_agg)
print("[validate] original CSV Epithelial rows:\n", orig_epi[RECOMPUTED_COLS])
if not (full_agg["cell_num"].to_numpy() == orig_epi["cell_num"].to_numpy()).all():
    raise ValueError("recomputed cell_num does not match original CSV")
for col in ["nCount_RNA", "nFeature_RNA", "percent.mt"]:
    if not np.allclose(full_agg[col].to_numpy(), orig_epi[col].to_numpy(), rtol=REL_TOL):
        raise ValueError(f"recomputed {col} does not reproduce original CSV (rtol={REL_TOL})")
print("[validate] OK -- method reproduces the original Epithelial covariates")

# ---- Exclude sub-cluster 24 and check against Memento's barcodes ------------
keep = (subcl != EXCLUDED_SUBCLUSTER).to_numpy()
print(f"[excl] dropping {(~keep).sum()} {EXCLUDED_SUBCLUSTER} cells; keeping {keep.sum()}")
excl_obs_names = ad.read_h5ad(EXCL_H5AD, backed="r").obs_names
if set(qc.index[keep]) != set(excl_obs_names):
    raise ValueError("retained barcodes differ from epithelial_excl_cluster24_raw_counts.h5ad")
excl_agg = aggregate(qc[keep])
print("[excl] sub-cluster-24-excluded covariates:\n", excl_agg)

# ---- Write patched covariate CSV -------------------------------------------
out = orig.copy()
other_cols = [c for c in out.columns
              if c not in RECOMPUTED_COLS + [CELL_TYPE_COL, INDIVIDUAL_COL, "sample_id"]]
for ind in INDIVIDUALS:
    row = (out[CELL_TYPE_COL] == TARGET_CELL_TYPE) & (out[INDIVIDUAL_COL] == ind)
    if row.sum() != 1:
        raise ValueError(f"expected 1 row for {ind} {TARGET_CELL_TYPE}, found {row.sum()}")
    for col in RECOMPUTED_COLS:
        out.loc[row, col] = excl_agg.loc[ind, col]
    out.loc[row, other_cols] = np.nan
out.to_csv(OUT_CSV, index=False)
print(f"[done] wrote {OUT_CSV} ({out.shape[0]} rows)")
