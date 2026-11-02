#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import scanpy as sc
import memento

REPO_DIR = "/project/gilad/brendan/dispersion/pilot/cHDC_data"
INPUT_H5AD = os.path.join(
    REPO_DIR, "analysis/epithelial_subcluster24_excl_sensitivity/outputs",
    "epithelial_excl_cluster24_raw_counts.h5ad")
OUT_DIR = os.path.join(REPO_DIR, "figure_reproduction/outputs/basics_scran_crossmethod")
MEM_OUT_DIR = os.path.join(OUT_DIR, "memento_1000subset")
BARCODES_CSV = os.path.join(OUT_DIR, "Epithelial_cells_1000subset_barcodes.csv")

INDIVIDUAL_COL = "vireo.individual"
CELL_TYPE_COL = "new.predicted.id"
INDIVIDUALS = ["NA19099", "NA19093", "NA18522"]

THRESHOLD = 1000
RANDOM_SEED = 42          # same as reproduce_1000subset.py / sweep_brendan_cell_count.py
MIN_PERC_GROUP = 0.95     # same as reproduce_1000subset.py (matched-1000cell comparator convention)
CAPTURE_RATE = 0.25       # same as reproduce_1000subset.py

os.makedirs(MEM_OUT_DIR, exist_ok=True)


def downsample_to(adata_ct, n, rng):
    """Same logic as reproduce_1000subset.py's downsample_to."""
    keep = []
    for indiv in adata_ct.obs[INDIVIDUAL_COL].unique():
        pos = np.where(adata_ct.obs[INDIVIDUAL_COL] == indiv)[0]
        if len(pos) < n:
            continue
        keep.extend(rng.choice(pos, size=n, replace=False).tolist())
    assert len(keep) > 0, f"No individuals survived downsampling to {n} cells"
    return adata_ct[keep].copy()


def rename_memento_cols(df):
    if "gene" in df.columns:
        df = df.set_index("gene")
    new_cols = []
    for col in df.columns:
        parts = col.split("^")
        assert len(parts) == 2, f"Unexpected Memento column format: {col!r}"
        new_cols.append(parts[1])
    df.columns = new_cols
    return df


def main():
    assert os.path.exists("CLAUDE.md"), "must run from project root"
    adata = sc.read(INPUT_H5AD)
    print(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")
    assert INDIVIDUAL_COL in adata.obs.columns

    n_before = {ind: int((adata.obs[INDIVIDUAL_COL] == ind).sum()) for ind in INDIVIDUALS}
    for ind, n in n_before.items():
        print(f"  {ind}: {n} cells before downsampling")
        assert n >= THRESHOLD, f"{ind} has fewer than {THRESHOLD} cells after cluster-24 exclusion"

    # Fresh RNG (seed=42) -- this population is Epithelial_cells only, so there
    # is no shared-CELL_TYPES-order RNG-advance concern the original script had
    # (that concern only applies when multiple cell types share one RNG stream
    # from one AnnData load; here there is exactly one cell type).
    rng = np.random.default_rng(RANDOM_SEED)
    adata_sub = downsample_to(adata, THRESHOLD, rng)
    assert adata_sub.shape[0] == 3 * THRESHOLD, (
        f"expected {3 * THRESHOLD} cells, got {adata_sub.shape[0]}")
    print(f"Downsampled to {adata_sub.shape[0]} cells")

    pd.Series(adata_sub.obs_names, name="barcode").to_csv(BARCODES_CSV, index=False)
    print(f"Wrote barcodes: {BARCODES_CSV}")

    adata_sub.obs["capture_rate"] = CAPTURE_RATE
    memento.setup_memento(adata_sub, q_column="capture_rate")
    memento.create_groups(adata_sub, label_columns=[INDIVIDUAL_COL])
    memento.compute_1d_moments(adata_sub, min_perc_group=MIN_PERC_GROUP)
    _, dispersion_new, _ = memento.get_1d_moments(adata_sub)
    dispersion_new = rename_memento_cols(dispersion_new)
    dispersion_new.columns = [f"{ind}_Epithelial_cells" for ind in dispersion_new.columns]
    dispersion_new.index.name = "gene"

    out_csv = os.path.join(MEM_OUT_DIR, "dispersion_by_ind_Epithelial_cells.csv")
    dispersion_new.reset_index().to_csv(out_csv, index=False)
    print(f"Wrote {dispersion_new.shape[0]} genes x {dispersion_new.shape[1]} individuals to {out_csv}")


if __name__ == "__main__":
    main()
