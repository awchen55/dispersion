# Radial Glia Downsampling Code (Recapitulation Rate)

Scripts that generate the downsampling data in
`../downsample_recapitulation_rate_radia_glia/`, used to measure how DD
(differential dispersion) recapitulation rate scales with allotetraploid
cell count. Orchestrated by `../submit_downsample_jobs.sh` (SLURM array
job, one task per cell-count x replicate combination -- see that script for
the exact grid and replicate count actually run).

## `downsample_memento_dd.py`

For one allotetraploid cell type (Radial glia), downsamples to a fixed
number of physical hybrid cells and reruns the production memento DD test
(`binary_test_1d`) on just that downsampled cell type. Mirrors
`de_dd_snakemake_mt_rb_mean_filtered/scripts/memento_tetraploid_dd_ct.py`
exactly (same `q`, `min_perc_group`, `mean_expr_cutoff`, mito/ribosomal
gene filtering) except for the downsampling step itself.

Each physical hybrid cell contributes two rows in the source h5ad -- a
`human_<barcode>` row (human-allele counts) and a `chimp_<barcode>` row
(chimp-allele counts) sharing the same barcode suffix -- so downsampling
selects whole physical cells (both rows together), not rows independently;
otherwise memento's paired-species test design breaks.

```
python downsample_memento_dd.py \
  --h5ad <tetraploid Radial_glia h5ad> --celltype Radial_glia \
  --n-cells <N> --seed <seed> \
  --q 0.2 --min-perc-group 0.95 --mean-expr-cutoff -6 \
  --num-cpus 8 --num-boot 5000 \
  --out-csv <path>
```

## `merge_mash_recap.R`

Takes one downsampled run's DD result CSV (from the script above) and:
1. Rebuilds the full 21-cell-type allotetraploid logFC/SE/p-value matrices
   (mirroring `tetraploid_dd_merge.R`), substituting in the downsampled
   run's CSV for the target cell type and reading the unchanged production
   CSVs for every other cell type.
2. Refits MASH on those matrices (mirroring `tetraploid_dd_mash.R`),
   restricted to the same gene list used in production
   (`dd/diploid_de_dd_results_filtered.csv`).
3. Computes the recapitulation rate for the target cell type: the fraction
   of the (fixed, production, unchanged) diploid DD gene set that is also
   DD-significant with matching effect-size sign in this newly-refit
   allotetraploid MASH result.

Diploid data is never touched -- only the target cell type's allotetraploid
arm is downsampled -- so the diploid DD gene set (the denominator of the
recapitulation rate) is held fixed across every point on the resulting
curve.

```
Rscript merge_mash_recap.R \
  --target-ct Radial_glia --target-csv <path from downsample_memento_dd.py> \
  --n-cells <N> --replicate <R> \
  --out-csv <path> \
  [--prod-dir <production q0.2 results dir, defaults to the standard path>] \
  [--pca-components 5] [--lfsr-threshold 0.05] [--save-mash]
```

## Software versions

Same Python environment (`/project/gilad/awchen55/envs/memento`) and R
environment (`module load R/4.4.1 gsl hdf5/1.12.0`,
`R_LIBS_USER=/project/gilad/awchen55/Rlibs/R4-4-1`) as
`code/de_dd_snakemake_mt_rb_mean_filtered/README.md` -- see that file for
the full Python/R package version table. One additional R package used
only here:

| | Version |
|---|---|
| optparse | 1.7.5 |
