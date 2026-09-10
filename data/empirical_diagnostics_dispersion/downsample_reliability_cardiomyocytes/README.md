# Downsample Reliability — Cardiomyocytes

Data for assessing the sensitivity of memento's dispersion (residual
variance) estimates to per-individual cell count. Cardiomyocytes is the
largest cell type in the cHDC cardiac pilot dataset (30,423 cells pooled
across 3 individuals: NA19099, NA18522, NA19093), chosen so that even the
largest downsampling depths tested here stay well within the available
data.

For each of 6 fixed depths (100, 250, 500, 1000, 2500, 5000 cells per
individual), each individual is independently downsampled to that many
cells (individuals with fewer cells than the target depth are excluded for
that depth), memento is rerun on the downsampled data, and the resulting
per-gene mean/dispersion estimates are saved. This is repeated for 30
independent random replicates per depth to separate technical sampling
noise (spread across replicates at a fixed depth) from systematic bias
(drift of the replicate average away from the full-data reference). A
`depth_full` version (no downsampling, one run) is included as the
reference every downsampled estimate is compared against.

See `code/empirical_diagnostics_dispersion/README.md` for the script that
generates this data and how to regenerate it.

## Contents

- `depth_{100,250,500,1000,2500,5000}/` — 30 replicates each:
  `mean_rep{0-29}.csv`, `dispersion_rep{0-29}.csv`. Columns are
  `{individual}_Cardiomyocytes` (one column per surviving individual at
  that depth); rows are genes.
- `depth_full/` — the single full-data reference: `mean_rep0.csv`,
  `dispersion_rep0.csv`.
- `sweep_summary.csv` — one row per (depth, replicate): cell/individual
  counts, gene count, and success/failure status. All 61 runs in this
  dataset succeeded (0 failures).
- `sweep.log` — full run log (timing, per-replicate cell/individual/gene
  counts).

## Notes

- Capture rate assumption: 0.25 (10x Chromium v3, matching the rest of this
  project). `min_perc_group = 0.95` (memento's `compute_1d_moments`),
  matching the value used for the project's main dispersion CSVs.
- At smaller depths, individuals with fewer cells than the target depth are
  dropped entirely for that depth/replicate (not padded or oversampled) —
  `n_individuals` in `sweep_summary.csv` shows how many of the 3 individuals
  actually contributed at each depth.
- Column names are `{individual}_{celltype}` (e.g. `NA18522_Cardiomyocytes`),
  reformatted from memento's native `sg^{individual}^{celltype}` naming.
