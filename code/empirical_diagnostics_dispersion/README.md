# Empirical Diagnostics: Dispersion

Code for empirically characterizing how reliable memento's dispersion
(residual variance) estimates are as a function of per-individual cell
count.

## `downsample_reliability_sweep.py`

Takes Cardiomyocytes (the largest cell type in the cHDC cardiac pilot
dataset by total pooled cells across individuals), downsamples each
individual's cells to a series of fixed depths, reruns memento at each
depth with 30 independent replicate subsamples, and saves the resulting
per-gene mean/dispersion estimates. The full-data (no downsampling) case is
run once as the reference. This produces the raw inputs for a
reliability-vs-depth analysis: comparing each replicate's dispersion
estimate against the full-data reference (bias) and across replicates at a
fixed depth (technical noise).

Reads raw counts read-only from the cHDC pilot h5ad
(`/project/gilad/brendan/dispersion/pilot/cHDC_data/...`); writes only
under this project's own `data/` directory. See
`data/empirical_diagnostics_dispersion/downsample_reliability_cardiomyocytes/README.md`
for the resulting data layout.

```
/project/gilad/awchen55/envs/memento/bin/python downsample_reliability_sweep.py
```

Key parameters (hardcoded at the top of the script):
- `CELL_TYPE = "Cardiomyocytes"`
- `DEPTHS = [100, 250, 500, 1000, 2500, 5000, None]` (`None` = full data)
- `REPLICATE_SEEDS = [42..71]` (30 replicates per depth)
- `MIN_PERC_GROUP = 0.95`, `CAPTURE_RATE = 0.25` (the cHDC cardiac pilot was
  sequenced on the 10x GEM-X platform)

## Software versions

Python environment: `/project/gilad/awchen55/envs/memento`. Same
package/version conventions as
`code/de_dd_snakemake_mt_rb_mean_filtered/README.md` (memento-de, numpy,
pandas, scanpy, anndata) — see that file for the full version table.
