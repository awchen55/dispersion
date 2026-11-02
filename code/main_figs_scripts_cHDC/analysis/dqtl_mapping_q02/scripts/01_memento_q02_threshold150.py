#!/usr/bin/env python3

import sys
import logging
import pandas as pd
import scanpy as sc
import memento
from pathlib import Path

# ── Configuration ─────────────────────────────────────────────────────────────

DATA_PATH = Path(
    "/project/gilad/brendan/dispersion/pilot/cHDC_data"
    "/cellranger_cluster-mode_trial/analysis/datasets/memento/erik_anndata_objects"
)
OUTPUT_BASE = Path(__file__).resolve().parent.parent / "outputs" / "cell_count_threshold_sweep_q02"

# ANALYSIS_OK[magic-threshold]: unchanged from sweep_cell_count_threshold.py --
# same min_perc_group determined there by sweep_min_perc_group.py; q does not
# affect that choice.
MIN_PERC_GROUP = 0.95

# 2026-09-24: Brendan's requested rerun value, replacing the manuscript's
# canonical 0.25 (see .living/decisions.md, 2026-08-24 entry, for how 0.25
# was confirmed as canonical -- this script deliberately deviates from it).
CAPTURE_RATE = 0.2

CELL_TYPES = [
    "Cardiomyocytes",
    "Mesenchymal_cells",
    "Cardiac_fibroblasts",
    "Epithelial_cells",
    "Activated_fibroblasts",
    "Foregut_endodermal_cells",
]

# Scope pinned to threshold=150 only -- see module docstring.
CELL_COUNT_THRESHOLDS = [150]


# ── Logging ───────────────────────────────────────────────────────────────────

def setup_logging(log_path: Path) -> logging.Logger:
    log = logging.getLogger("memento_q02_threshold150")
    log.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s  %(levelname)-8s  %(message)s")
    for handler in [logging.StreamHandler(sys.stdout), logging.FileHandler(log_path)]:
        handler.setFormatter(fmt)
        log.addHandler(handler)
    return log


# ── Per-cell-type helpers ─────────────────────────────────────────────────────

def load_and_filter(cell_type: str, min_cells: int, log: logging.Logger) -> sc.AnnData:
    path = DATA_PATH / f"erik_{cell_type}_counts_for_adata.h5ad"
    assert path.exists(), f"H5AD not found: {path}"

    adata = sc.read(str(path))
    n_cells_raw, n_genes = adata.shape
    log.info(f"    loaded {n_cells_raw:,} cells x {n_genes:,} genes")

    indiv_counts = adata.obs["vireo.individual"].value_counts()
    valid = indiv_counts[indiv_counts >= min_cells].index
    n_before, n_after = len(indiv_counts), len(valid)
    log.info(f"    {n_after}/{n_before} individuals pass >= {min_cells} cell threshold")

    assert n_after > 0, f"No individuals pass cell-count threshold for {cell_type}"

    adata_f = adata[adata.obs["vireo.individual"].isin(valid)].copy()
    log.info(f"    {adata_f.shape[0]:,} cells retained")

    assert adata_f.shape[0] > 0, f"No cells remain after filtering for {cell_type}"
    return adata_f


def run_memento(
    adata: sc.AnnData, min_perc_group: float
) -> tuple[pd.DataFrame, pd.DataFrame]:
    adata.obs["capture_rate"] = CAPTURE_RATE
    memento.setup_memento(adata, q_column="capture_rate")
    memento.create_groups(adata, label_columns=["vireo.individual", "cell_type"])
    memento.compute_1d_moments(adata, min_perc_group=min_perc_group)
    mean, dispersion, _ = memento.get_1d_moments(adata)
    assert mean.shape == dispersion.shape, (
        f"Shape mismatch: mean {mean.shape} vs dispersion {dispersion.shape}"
    )
    return mean, dispersion


# ── Main sweep ────────────────────────────────────────────────────────────────

def main() -> None:
    assert MIN_PERC_GROUP is not None, (
        "MIN_PERC_GROUP is not set. Review sweep_min_perc_group results and set "
        "the value at the top of this script before running."
    )

    OUTPUT_BASE.mkdir(parents=True, exist_ok=True)
    log = setup_logging(OUTPUT_BASE / "sweep.log")

    log.info("=" * 60)
    log.info(f"Memento q=0.2 rerun, cell_count_threshold = {CELL_COUNT_THRESHOLDS}")
    log.info(f"min_perc_group fixed at {MIN_PERC_GROUP}")
    log.info(f"Cell types: {CELL_TYPES}")
    log.info("=" * 60)

    summary_rows = []

    for threshold in CELL_COUNT_THRESHOLDS:
        threshold_dir = OUTPUT_BASE / f"cell_count_{threshold}"
        threshold_dir.mkdir(exist_ok=True)
        log.info(f"\ncell_count_threshold = {threshold}")
        all_ok = True

        for ct in CELL_TYPES:
            log.info(f"  {ct}")
            try:
                adata = load_and_filter(ct, threshold, log)
                mean, dispersion = run_memento(adata, MIN_PERC_GROUP)

                n_genes, n_indivs = mean.shape
                log.info(f"    SUCCESS: {n_genes:,} genes x {n_indivs} individuals")

                mean.to_csv(threshold_dir / f"mean_by_ind_{ct}.csv")
                dispersion.to_csv(threshold_dir / f"dispersion_by_ind_{ct}.csv")

                summary_rows.append({
                    "cell_count_threshold": threshold, "cell_type": ct,
                    "status": "success", "n_genes": n_genes,
                    "n_individuals": n_indivs, "error": "",
                })

            # ANALYSIS_OK[broad-except]: intentional -- sweep catches any failure
            # (mirrors sweep_cell_count_threshold.py's identical waiver) without
            # knowing in advance which exception memento raises. All failures are
            # logged and recorded in the summary CSV for post-run inspection.
            except Exception as e:
                log.warning(f"    FAILED: {e}")
                all_ok = False
                summary_rows.append({
                    "cell_count_threshold": threshold, "cell_type": ct,
                    "status": "failed", "n_genes": None,
                    "n_individuals": None, "error": str(e),
                })

        if not all_ok:
            log.warning(f"  One or more cell types failed at threshold={threshold}")

    # ── Summary ───────────────────────────────────────────────────────────────
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUTPUT_BASE / "sweep_summary.csv", index=False)
    log.info(f"\nSummary written to {OUTPUT_BASE / 'sweep_summary.csv'}")

    pivot = summary.pivot(
        index="cell_count_threshold", columns="cell_type", values="status"
    )
    pivot_path = OUTPUT_BASE / "success_by_cell_type.csv"
    pivot.to_csv(pivot_path)
    log.info(f"Per-cell-type success table written to {pivot_path}")
    log.info("\n" + pivot.to_string())

    n_indiv_pivot = summary.pivot(
        index="cell_count_threshold", columns="cell_type", values="n_individuals"
    )
    n_indiv_path = OUTPUT_BASE / "n_individuals_by_cell_type.csv"
    n_indiv_pivot.to_csv(n_indiv_path)
    log.info(f"Per-cell-type individual counts written to {n_indiv_path}")
    log.info("\n" + n_indiv_pivot.to_string())


if __name__ == "__main__":
    main()
