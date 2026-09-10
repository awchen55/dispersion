# DD/DE Interspecies Pipeline (mt/ribo + mean-expression filtered)

Snakemake pipeline that tests for differential dispersion (DD) and
differential expression (DE) between human and chimpanzee, in both diploid
cell lines and allotetraploid (fused) hybrid cells, across the same 20 cell
types in both ploidy contexts, at a sweep of `memento` capture-efficiency
values (`q`). Both
mitochondrial/ribosomal genes and a log-mean-expression cutoff are applied
to every `memento` variance/mean estimate (see "Gene filtering" below) --
that combination is what distinguishes this pipeline from the sibling
`de_dd_snakemake` and `de_dd_snakemake_mt_rb_filtered` trees, which apply
neither or only one of the two filters.

Final output: one `DD_DE_interspecies_results.csv` per `q` value, under
`data/dd_de_results/mt_rb_mean_filtered_q_thresholds/q{q}/`.

---

## Pipeline stages

Run via `snakemake --profile profile` (see `submit.sh`); `rule all` in the
`Snakefile` requests one master CSV per `q` in `config["memento_q"]`
(`[0.05, 0.1, 0.15, 0.2]`).

### 1. Memento variance/mean/DD estimation (`scripts/memento_diploid.py`, `scripts/memento_tetraploid_dd_ct.py`)

- **Diploid** (`memento_diploid_ct` rule): for each cell type, groups cells
  by `species x donor_id` and runs `memento.compute_1d_moments()` to get
  each donor's mean expression and (mean-corrected) dispersion per gene.
  Donor x cell-type groups with fewer than `memento_min_cells_per_group`
  (100) cells are dropped first (memento's moment estimator becomes
  numerically unstable below that). No DD test is run here -- that happens
  downstream in `diploid_dd_limma_ct.R`.
- **Tetraploid** (`memento_tetraploid_dd_ct` rule): for each cell type,
  groups cells by `species x celltype` (no donor/replicate level exists for
  the allotetraploid data) and runs memento's `binary_test_1d` directly to
  test for a species effect on dispersion (DD) and mean (DE) in one step.
- **Gene filtering** (both): mitochondrial (`MT-*`) and ribosomal protein
  (`RPL*`/`RPS*`, excluding the `RPS6K*` kinase family, which are kinases
  not ribosomal proteins despite the name) genes are dropped, and any gene
  with log mean expression above `memento_mean_expr_cutoff` (-6) in *any*
  individual species x sample group is also dropped (memento's mean
  correction is unreliable at the high-expression tail -- see
  `dispersion_diagnostics.Rmd` Sections 6-7 in the `differentialDispersion`
  project). Both filters are applied *after* memento's own moment
  computation and mean-variance trend fit, not before -- i.e. these
  excluded genes still inform the fitted trend used to score every
  retained gene's dispersion, but are then dropped from the output.

### 2. Diploid DE (`diploid_pseudobulk_ct.R` -> `diploid_de_limma_ct.R` -> `diploid_de_merge.R`)

Pseudobulk (summed counts) per `donor_id x species`, filtered to the same
`memento_min_cells_per_group` threshold for consistency with the DD side,
then a standard `limma`/`dream` pipeline: `voomWithDreamWeights`
(cyclicloess normalization) -> `dream` (human vs. chimp contrast) ->
`eBayes`. Merged across cell types into logFC/SE/p-value/adjusted-p-value
matrices. Independent of `q` (memento isn't involved in DE).

### 3. Diploid DD (`diploid_dd_limma_ct.R` -> `diploid_dd_merge.R` -> `diploid_gene_filter.R` -> `diploid_dd_mash.R`)

Runs `limma::normalizeQuantiles()` on the per-cell-type dispersion matrix
(genes x donors, quantile-normalized *across donor columns* to remove
donor-to-donor scale/shape differences before testing), then the same
`dream`/`eBayes` species contrast used for DE. Cell types are merged into
matrices, then filtered (`diploid_gene_filter.R`, see "Gene filter for
MASH" below) before being passed to `mashr` (`mash_1by1` -> `cov_pca` +
`cov_ed` -> `mash`) for multivariate shrinkage/significance (LFSR) across
cell types.

### 4. Tetraploid DD (`tetraploid_dd_merge.R` -> `diploid_gene_filter.R` -> `tetraploid_dd_mash.R`)

Per-cell-type DD results from step 1 (memento's `binary_test_1d`) are
merged into matrices, gene-filtered, and passed through `mashr` the same
way as the diploid side.

### 5. Gene filter for MASH (`diploid_gene_filter.R`)

Applied once, shared by both diploid and tetraploid MASH runs. A gene must
be "tested" (`logFC != 0`) in at least `gene_filter_threshold` (70%) of
cell types -- tetraploid filter applied first, then the diploid filter
applied to the tetraploid-passing gene set. Final gene list = intersection
of both, with untested genes imputed (logFC=0, SE=1e16, p=1) so `mashr`
still receives complete matrices.

### 6. Tetraploid DE (`tetraploid_de_wilcoxon_ct.R`)

Paired Wilcoxon signed-rank test (`matrixTests::row_wilcoxon_paired`) on
raw counts, human vs. chimp, per cell type. Independent of `q`.

### 7. Integration (`integrate_dd_de.R`)

Combines all of the above into one master
`DD_DE_interspecies_results.csv` per `q`, restricted to genes tested for DD
in both diploid and allotetraploid. Also derives, per cell type common to
both ploidy contexts:
- `dd.recapitulated.same.sign.<ct>` -- DD-significant (LFSR < 0.05) with
  matching effect-size sign in both diploid and allotetraploid.
- `dd.cis.proportion.<ct>` -- the allotetraploid (cis-only, since the two
  parental genomes share one nucleus/trans environment) effect size as a
  fraction of the total diploid effect size; the remainder is attributed to
  trans effects.
- `dd.cis.category.<ct>` / `dd.cis.category.overall` /
  `dd.cis.category.overall.mean.definition` -- cis (>0.70) / trans (<0.30) /
  neither classification, per cell type and pooled across cell types (by
  per-celltype vote or by the mean cis proportion across diploid-DD cell
  types).

---

## Standalone analysis (not part of `rule all`)

### `memento_variance_decomposition.py` (+ `submit_variance_decomposition.sh`)

Not wired into the Snakefile DAG -- run manually via the submit script. For
every (diploid cell type x species x donor) and (tetraploid cell type x
species) group, re-runs memento's own moment computation with identical
parameters to the production pipeline above, but pulls the raw 2nd-moment
variance out of `adata.uns['memento']['1d_moments']` before memento's
normal API (`get_1d_moments()`) discards it -- that raw variance is never
written to any of the production `_var_*.csv` files, which store only
`log(residual variance)` (dispersion). Outputs
`all_celltypes_variance_genes.csv` / `all_celltypes_variance_summary.csv`
under `data/dd_de_results/mt_rb_mean_filtered_q_thresholds/q0.2/variance_decomposition/`,
used by the "Total variability decomposition" section of
`differentialDispersion/analysis/dispersion_diagnostics/dispersion_diagnostics.Rmd`
to characterise how much of a gene's total variability memento's dispersion
estimate captures, independent of the production DD significance test (which
runs on quantile-normalized dispersion -- see that Rmd section for the
distinction).

---

## Software versions

### Python (conda env `/project/gilad/awchen55/envs/memento`, used by all `memento_*.py` scripts)

| | Version |
|---|---|
| Python | 3.11.15 |
| memento-de | 0.1.2 |
| numpy | 2.4.6 |
| pandas | 2.3.3 |
| scipy | 1.17.1 |
| scanpy | 1.11.5 |
| anndata | 0.12.18 |

### R (`module load R/4.4.1 gsl hdf5/1.12.0`; `R_LIBS_USER=/project/gilad/awchen55/Rlibs/R4-4-1`, used by all `*.R` scripts)

| | Version |
|---|---|
| R | 4.4.1 (2024-06-14) |
| Matrix | 1.7.3 |
| Seurat | 5.3.0 |
| SeuratObject | 5.2.0 |
| dplyr | 1.1.4 |
| tidyr | 1.3.1 |
| dreamlet | 1.4.1 |
| variancePartition | 1.36.3 |
| edgeR | 4.4.2 |
| limma | 3.62.2 |
| mashr | 0.2.79 |
| matrixTests | 0.2.3.1 |

Versions above were queried directly from the installed environments
(`python3 --version`, `packageVersion()` per package) on 2026-09-10 --


---

## Key parameters (`config.yaml`)

| Parameter | Value | Notes |
|---|---|---|
| `memento_q` | 0.05, 0.1, 0.15, 0.2 | assumed UMI capture efficiency; one full pipeline run per value |
| `memento_min_perc_group` | 0.95 | fraction of groups a gene must be detected in to be kept by memento |
| `memento_mean_expr_cutoff` | -6 | log-mean-expression cutoff (per individual sample group) |
| `memento_min_cells_per_group` | 100 | diploid only; smaller species x donor groups are dropped before memento runs |
| `memento_num_cpus` / `memento_num_boot` | 8 / 5000 | tetraploid `binary_test_1d` parallelism / bootstrap replicates |
| `limma_voomspan` | 0.5 | `voomWithDreamWeights` LOWESS span |
| `limma_min_count` | 5 | `filterByExpr` minimum count (DE only) |
| `limma_ddf` | Satterthwaite | `dream` degrees-of-freedom method |
| `mash_pca_components` | 5 | `cov_pca` components for `mashr` |
| `mash_lfsr_threshold` | 0.05 | significance threshold used both for `mashr`'s "strong" subset and for all downstream DD/DE calls |
| `gene_filter_threshold` | 0.70 | fraction of cell types a gene must be tested in (both ploidy contexts) to enter MASH |

Cell type lists (`celltypes_diploid` and `celltypes_tetraploid`, the same 20
cell types in the same order -- `Early_mesoderm_cells`, which has no
diploid counterpart, was previously included in `celltypes_tetraploid` only
and has been removed for consistency) are in `config.yaml`.
