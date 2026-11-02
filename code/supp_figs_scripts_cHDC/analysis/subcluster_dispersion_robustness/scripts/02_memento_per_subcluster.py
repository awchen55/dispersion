#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import scanpy as sc
import memento
from scipy.sparse import csr_matrix, issparse

# ---- Named constants --------------------------------------------------------
ANALYSIS_DIR = "/project/gilad/brendan/dispersion/pilot/cHDC_data/analysis/subcluster_dispersion_robustness"
BASE_OUT_DIR = os.path.join(ANALYSIS_DIR, "outputs")
H5AD_PATH = os.path.join(BASE_OUT_DIR, "subcluster_raw_counts.h5ad")

INDIVIDUAL_COL = "vireo.individual"

# Sub-cluster obs columns written by 01_subcluster_and_export.R, one per
# resolution. Keys are human-readable resolution labels used in output paths.
RESOLUTION_COLUMNS = {
    "0.5": "subcluster_res_0p5",
    "1.0": "subcluster_res_1",
    "2.0": "subcluster_res_2",
}

# Memento parameters (match memento_pip_forloop.py).
CAPTURE_RATE = 0.25
MIN_PERC_GROUP = 0.95

# Minimum cells per (individual x sub-cluster) group for a group to be kept.
# Originally set to 1000 (Brendan 2026-07-20): "start at 1000 and see what we
# get" -- this matches the manuscript's >=1000-cells-per-individual inclusion
# threshold for the 6 main cell types, so per-sub-cluster dispersion estimates
# are on the same footing as the published cell-type-level estimates.
# Configurable via env var since the "Minimum sub-cluster cell count"
# diagnostic in SUBCLUSTER_DISPERSION_ROBUSTNESS.md motivated running more
# than one threshold (Brendan 2026-08-11: 1000 and 700, both cross-checked
# against explore_cell_count_thresholds.Rmd before picking).
MIN_CELLS_PER_GROUP = int(os.environ.get("SUBCLUSTER_MIN_CELLS_PER_GROUP", "1000"))

# Minimum number of individuals (of 3) that must clear MIN_CELLS_PER_GROUP for
# a sub-cluster to be kept at all. A sub-cluster where only 1 individual has
# enough cells can't support even a 2-individual comparison, so it is dropped
# in full (see the two-stage filter in main()). Brendan 2026-08-11: both
# configs run so far use minind=2.
MIN_INDIVIDUALS_PER_SUBCLUSTER = int(os.environ.get("SUBCLUSTER_MIN_INDIVIDUALS", "2"))

CONFIG_LABEL = f"mincells{MIN_CELLS_PER_GROUP}_minind{MIN_INDIVIDUALS_PER_SUBCLUSTER}"
OUT_DIR = os.path.join(BASE_OUT_DIR, CONFIG_LABEL)


def log_msg(msg):
    print(f"[02_memento] {msg}", flush=True)


def assert_raw_counts(adata):
    """Fail loudly unless adata.X is an integer-valued, non-negative CSR matrix."""
    # ANALYSIS_OK[layer-choice]: .X is the raw counts written by
    # 01_subcluster_and_export.R (Seurat counts slot -> AnnData X); this
    # function exists to prove that (integer/non-negative/CSR checks below).
    x = adata.X
    if not issparse(x):
        raise ValueError(
            "adata.X is dense; expected a sparse raw-counts matrix. Aborting."
        )
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
            "The Seurat->AnnData conversion did not preserve raw counts. Aborting."
        )
    log_msg(
        f"Raw-count assertion passed: CSR, non-negative, integer "
        f"(max deviation {max_int_err:.3g}); {adata.n_obs} cells x {adata.n_vars} genes."
    )


def main():
    if not os.path.exists(H5AD_PATH):
        raise FileNotFoundError(
            f"Missing input h5ad: {H5AD_PATH}. Run 01_subcluster_and_export.R first."
        )

    os.makedirs(OUT_DIR, exist_ok=True)
    log_msg(
        f"Config: MIN_CELLS_PER_GROUP={MIN_CELLS_PER_GROUP}, "
        f"MIN_INDIVIDUALS_PER_SUBCLUSTER={MIN_INDIVIDUALS_PER_SUBCLUSTER} "
        f"(outputs -> {OUT_DIR})"
    )

    log_msg(f"Reading {H5AD_PATH}")
    adata = sc.read(H5AD_PATH)
    log_msg(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes.")

    # Raw-count safety net (the whole point of the gotcha).
    assert_raw_counts(adata)

    if INDIVIDUAL_COL not in adata.obs.columns:
        raise KeyError(f"Missing obs column: {INDIVIDUAL_COL}")

    dropped_records = []

    for res_label, sc_col in RESOLUTION_COLUMNS.items():
        log_msg(f"=== Resolution {res_label} (obs column '{sc_col}') ===")
        if sc_col not in adata.obs.columns:
            raise KeyError(
                f"Missing sub-cluster obs column '{sc_col}' for resolution "
                f"{res_label}. Available: {list(adata.obs.columns)}"
            )

        # Stage 1: count cells per (individual x sub-cluster) and flag groups
        # that clear the per-group cell threshold.
        grp = (
            adata.obs.groupby([INDIVIDUAL_COL, sc_col], observed=True)
            .size()
            .reset_index(name="n_cells")
        )
        grp["clears_threshold"] = grp["n_cells"] >= MIN_CELLS_PER_GROUP

        # Stage 2: a sub-cluster is kept only if >= MIN_INDIVIDUALS_PER_SUBCLUSTER
        # individuals cleared stage 1; groups in a sub-cluster that doesn't meet
        # this bar are dropped even if that individual's own count is high enough.
        n_clearing_per_subcluster = grp.groupby(sc_col, observed=True)["clears_threshold"].sum()
        grp["n_individuals_clearing"] = grp[sc_col].map(n_clearing_per_subcluster)
        keep_mask = grp["clears_threshold"] & (
            grp["n_individuals_clearing"] >= MIN_INDIVIDUALS_PER_SUBCLUSTER
        )
        kept = grp[keep_mask]
        dropped = grp[~keep_mask]

        log_msg(
            f"  {len(grp)} (individual x sub-cluster) groups; "
            f"{keep_mask.sum()} kept (>= {MIN_CELLS_PER_GROUP} cells AND "
            f">= {MIN_INDIVIDUALS_PER_SUBCLUSTER} individuals clearing per sub-cluster), "
            f"{(~keep_mask).sum()} dropped."
        )
        for _, r in dropped.iterrows():
            if not r["clears_threshold"]:
                reason = f"< {MIN_CELLS_PER_GROUP} cells"
            else:
                reason = (
                    f">= {MIN_CELLS_PER_GROUP} cells but sub-cluster has only "
                    f"{int(r['n_individuals_clearing'])} individual(s) clearing "
                    f"(< {MIN_INDIVIDUALS_PER_SUBCLUSTER} required)"
                )
            log_msg(
                f"    DROP resolution={res_label} individual={r[INDIVIDUAL_COL]} "
                f"subcluster={r[sc_col]} n_cells={r['n_cells']} ({reason})"
            )
            dropped_records.append(
                {
                    "resolution": res_label,
                    "individual": r[INDIVIDUAL_COL],
                    "subcluster": r[sc_col],
                    "n_cells": int(r["n_cells"]),
                    "reason": reason,
                }
            )

        if len(kept) == 0:
            raise ValueError(
                f"No (individual x sub-cluster) groups survive the "
                f"{MIN_CELLS_PER_GROUP}-cell / {MIN_INDIVIDUALS_PER_SUBCLUSTER}-individual "
                f"filter at resolution {res_label}."
            )

        # Keep only cells belonging to surviving (individual x sub-cluster) pairs.
        kept_pairs = set(
            zip(kept[INDIVIDUAL_COL].astype(str), kept[sc_col].astype(str))
        )
        pair_series = list(
            zip(adata.obs[INDIVIDUAL_COL].astype(str), adata.obs[sc_col].astype(str))
        )
        cell_mask = np.array([p in kept_pairs for p in pair_series])
        log_msg(
            f"  Retaining {int(cell_mask.sum())} of {adata.n_obs} cells after "
            f"group-size filtering."
        )
        adata_sub = adata[cell_mask].copy()

        # Re-assert raw counts on the subset (conversions/copies can surprise).
        assert_raw_counts(adata_sub)

        # ---- Memento pipeline (grouped by individual x sub-cluster) ----------
        adata_sub.obs["capture_rate"] = CAPTURE_RATE
        memento.setup_memento(adata_sub, q_column="capture_rate")
        memento.create_groups(adata_sub, label_columns=[INDIVIDUAL_COL, sc_col])
        memento.compute_1d_moments(adata_sub, min_perc_group=MIN_PERC_GROUP)

        mean, var, counts = memento.get_1d_moments(adata_sub)
        log_msg(
            f"  Memento moments: mean {mean.shape}, var {var.shape} "
            f"(genes x groups)."
        )

        mean_path = os.path.join(OUT_DIR, f"mean_expression_by_subcluster_res_{res_label}.csv")
        var_path = os.path.join(OUT_DIR, f"res_var_by_subcluster_res_{res_label}.csv")
        mean.to_csv(mean_path, index=True)
        var.to_csv(var_path, index=True)
        log_msg(f"  Wrote {mean_path}")
        log_msg(f"  Wrote {var_path}")

    dropped_path = os.path.join(OUT_DIR, "dropped_subclusters.csv")
    pd.DataFrame(
        dropped_records,
        columns=["resolution", "individual", "subcluster", "n_cells", "reason"],
    ).to_csv(dropped_path, index=False)
    log_msg(f"Wrote dropped-group log: {dropped_path} ({len(dropped_records)} rows).")
    log_msg("DONE. Next: scripts/03_robustness_correlations.R")


if __name__ == "__main__":
    main()
