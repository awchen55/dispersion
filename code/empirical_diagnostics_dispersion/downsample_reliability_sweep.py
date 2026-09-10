#!/usr/bin/env python3
"""
Downsample reliability sweep for memento dispersion estimates.

Takes the largest cell type (Cardiomyocytes, by total pooled cells across
individuals), downsamples each individual's cells to a series of fixed
depths, reruns memento at each depth with multiple independent replicate
subsamples, and saves the resulting mean/dispersion CSVs. The "full" depth
(no downsampling) is run once as the reference.

This produces the raw inputs for a reliability-vs-depth analysis in
dispersion_diagnostics.Rmd: comparing each replicate's per-gene dispersion
estimate against the full-data reference (bias) and across replicates at a
fixed depth (technical noise).

Reads raw counts read-only from Brendan's directory; writes only under this
project's own data/ directory.

Run with the memento conda environment:
    /project/gilad/awchen55/envs/memento/bin/python downsample_reliability_sweep.py
"""

import sys
import gc
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import memento
from pathlib import Path

# ── Configuration ──────────────────────────────────────────────────────────────

ADATA_PATH = Path(
    "/project/gilad/brendan/dispersion/pilot/cHDC_data"
    "/cellranger_cluster-mode_trial/analysis/datasets"
    "/100k_full_for_adata_conv_raw_counts_only.h5ad"
)
OUTPUT_BASE = Path(__file__).resolve().parent.parent / "data" / "downsample_reliability"

CELL_TYPE_COL  = "new.predicted.id"
INDIVIDUAL_COL = "vireo.individual"

# Largest cell type by total pooled cells across the 3 individuals.
CELL_TYPE = "Cardiomyocytes"

# Depths span well below the ~1000-cell inclusion threshold used elsewhere in
# this project, up to the largest available individual pool (~17k cells).
# None = full data for that individual, no downsampling (the reference).
DEPTHS = [100, 250, 500, 1000, 2500, 5000, None]

# Independent random subsamples per depth, to separate technical noise
# (spread across replicates) from bias (drift of the replicate average away
# from the full-data reference). "full" has no randomness, so only one run.
REPLICATE_SEEDS = [42 + i for i in range(30)]

MIN_PERC_GROUP = 0.95  # matches the value used for the project's main dispersion CSVs
CAPTURE_RATE   = 0.25  # cHDC cardiac pilot was sequenced on the 10x GEM-X platform


# ── Logging ────────────────────────────────────────────────────────────────────

def setup_logging(log_path: Path) -> logging.Logger:
    log = logging.getLogger("downsample_reliability_sweep")
    log.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s  %(levelname)-8s  %(message)s")
    for handler in [logging.StreamHandler(sys.stdout), logging.FileHandler(log_path)]:
        handler.setFormatter(fmt)
        log.addHandler(handler)
    return log


# ── Helpers ────────────────────────────────────────────────────────────────────

def downsample_to(adata_ct: sc.AnnData, n: int, rng: np.random.Generator) -> sc.AnnData:
    """Downsample each individual to exactly n cells; exclude individuals with < n."""
    keep = []
    for indiv in adata_ct.obs[INDIVIDUAL_COL].unique():
        pos = np.where(adata_ct.obs[INDIVIDUAL_COL] == indiv)[0]
        if len(pos) < n:
            continue
        keep.extend(rng.choice(pos, size=n, replace=False).tolist())
    assert len(keep) > 0, f"No individuals survived downsampling to {n} cells"
    return adata_ct[keep].copy()


def rename_memento_cols(df: pd.DataFrame) -> pd.DataFrame:
    """Rename 'sg^NA18522^Cardiomyocytes' -> 'NA18522_Cardiomyocytes'."""
    if "gene" in df.columns:
        df = df.set_index("gene")
    new_cols = []
    for col in df.columns:
        parts = col.split("^")
        assert len(parts) == 3, f"Unexpected Memento column format: {col!r}"
        new_cols.append(f"{parts[1]}_{parts[2]}")
    df.columns = new_cols
    return df


def run_memento(adata_ct: sc.AnnData, log: logging.Logger) -> tuple[pd.DataFrame, pd.DataFrame]:
    work = adata_ct.copy()
    work.obs["capture_rate"] = CAPTURE_RATE
    memento.setup_memento(work, q_column="capture_rate")
    memento.create_groups(work, label_columns=[INDIVIDUAL_COL, CELL_TYPE_COL])
    memento.compute_1d_moments(work, min_perc_group=MIN_PERC_GROUP)
    mean, dispersion, _ = memento.get_1d_moments(work)
    assert mean.shape == dispersion.shape, (
        f"Shape mismatch: mean {mean.shape} vs dispersion {dispersion.shape}"
    )
    mean, dispersion = rename_memento_cols(mean), rename_memento_cols(dispersion)
    log.info(f"    {mean.shape[0]:,} genes x {mean.shape[1]} individuals")
    del work
    gc.collect()
    return mean, dispersion


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    OUTPUT_BASE.mkdir(parents=True, exist_ok=True)
    log = setup_logging(OUTPUT_BASE / "sweep.log")

    log.info(f"Cell type: {CELL_TYPE}")
    log.info(f"Depths: {DEPTHS}")
    log.info(f"Replicate seeds: {REPLICATE_SEEDS}")

    log.info("Loading AnnData...")
    adata_full = sc.read(str(ADATA_PATH))
    log.info(f"  {adata_full.shape[0]:,} cells x {adata_full.shape[1]:,} genes")

    assert CELL_TYPE_COL in adata_full.obs.columns, f"Missing obs column: {CELL_TYPE_COL!r}"
    assert INDIVIDUAL_COL in adata_full.obs.columns, f"Missing obs column: {INDIVIDUAL_COL!r}"

    adata_ct = adata_full[adata_full.obs[CELL_TYPE_COL] == CELL_TYPE].copy()
    del adata_full
    log.info(f"  {adata_ct.shape[0]:,} cells after subsetting to {CELL_TYPE}")
    log.info("  Cells per individual:")
    for indiv, n in adata_ct.obs[INDIVIDUAL_COL].value_counts().items():
        log.info(f"    {indiv}: {n:,}")

    summary_rows = []

    for depth in DEPTHS:
        label = "full" if depth is None else str(depth)
        depth_dir = OUTPUT_BASE / f"depth_{label}"
        depth_dir.mkdir(exist_ok=True)

        seeds = [None] if depth is None else REPLICATE_SEEDS
        log.info(f"\n{'='*55}\nDepth: {label} ({len(seeds)} replicate(s))")

        for rep_idx, seed in enumerate(seeds):
            rep_label = "0" if seed is None else str(rep_idx)
            log.info(f"  replicate {rep_label} (seed={seed})")

            if depth is None:
                adata_rep = adata_ct
            else:
                rng = np.random.default_rng(seed)
                adata_rep = downsample_to(adata_ct, depth, rng)

            n_indiv = adata_rep.obs[INDIVIDUAL_COL].nunique()
            n_cells_rep = adata_rep.shape[0]
            log.info(f"    {n_cells_rep:,} cells, {n_indiv} individuals")

            mean = dispersion = None
            try:
                mean, dispersion = run_memento(adata_rep, log)
                mean.to_csv(depth_dir / f"mean_rep{rep_label}.csv")
                dispersion.to_csv(depth_dir / f"dispersion_rep{rep_label}.csv")
                status = "success"
                error = ""
                n_genes = mean.shape[0]
            # ANALYSIS_OK[broad-except]: sweep catches any per-(depth, replicate)
            # failure (e.g. too few individuals surviving at a given depth)
            # without knowing in advance which exception memento may raise; all
            # failures are logged and recorded in the summary CSV.
            except Exception as e:
                log.warning(f"    FAILED: {e}")
                status = "failed"
                error = str(e)
                n_genes = None

            summary_rows.append({
                "cell_type": CELL_TYPE, "depth": label, "replicate": rep_label,
                "seed": seed, "n_cells": n_cells_rep, "n_individuals": n_indiv,
                "status": status, "n_genes": n_genes, "error": error,
            })

            # adata_rep is a fresh downsampled copy each rep (except the "full"
            # case, which reuses adata_ct directly) -- release it explicitly so
            # memory doesn't creep up across the ~70 iterations in this sweep.
            if depth is not None:
                del adata_rep
            del mean, dispersion
            gc.collect()

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUTPUT_BASE / "sweep_summary.csv", index=False)
    log.info(f"\nSummary written to {OUTPUT_BASE / 'sweep_summary.csv'}")
    log.info("\n" + summary.to_string())
    log.info("\nDone.")


if __name__ == "__main__":
    main()
