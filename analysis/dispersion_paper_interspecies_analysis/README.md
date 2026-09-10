# Dispersion Paper — Interspecies Analysis

Paper-ready figures characterizing differential dispersion (DD) and differential expression (DE) between human and chimpanzee, in diploid cell lines and allotetraploid (fused) hybrid cells, across 20 cell types.


------------------------------------------------------------------------

## `dd_de_interspecies_analysis.Rmd`

**Input:** `data/dd_de_results/mt_rb_mean_filtered_q_thresholds/q0.2/DD_DE_interspecies_results.csv` (the master DD/DE table; see `code/de_dd_snakemake_mt_rb_mean_filtered/README.md` for how it's produced) — mitochondrial/ribosomal genes and a log-mean expression cutoff of -6 are already applied upstream of this table.

**DD definition used throughout:** `dd.lfsr.diploid.<celltype> < 0.05`. **Recapitulated:** also `dd.lfsr.allotetraploid.<celltype> < 0.05` with the same effect-size sign.

### Summary Tables / Plot of DD Genes

Number of DD genes per cell type (diploid), split by recapitulated vs. diploid-only, ordered by descending total.

### Recapitulation rate

-   Bar chart of recapitulation rate per cell type, annotated with each cell type's allotetraploid cell count, to see whether cell types with fewer cells also show lower rates.
-   Scatter of recapitulation rate vs. allotetraploid cell number across cell types (confounds cell number with other cell-type-specific biology).
-   **Downsampling within one cell type** (Radial glia, the largest allotetraploid cell type): isolates the effect of cell number alone by downsampling to 10 grid points (100-40,000 cells, 10 replicates each), rerunning the production memento DD test + MASH refit at each size, and recomputing recapitulation rate against the fixed diploid DD gene set. See `analysis/dd_de_analysis/downsample_recapitulation/` for the underlying data/code and their own READMEs.

### Examples of DD Genes

-   **DD but not DE** (*SESN3*, Neuroblast): density plot showing overlapping location (no DE) but differing shape (DD) between species.
-   **Mean expression and dispersion by species** (*SETD2*, Stromal cells): a second, independent DD-not-DE example, shown directly via memento's per-donor mean/dispersion estimates (box plots) rather than a density plot.
-   **DD and DE opposite signs** (*MAML2*, ENS glia): a gene significant for both DD and DE, where the direction of the two effects differs between species.

### Effect Size Plot

Scatter of diploid vs. allotetraploid DD effect size (logFC) across all cell types, for genes DD at a stricter LFSR \< 0.001 cutoff, highlighting genes recapitulated with the same sign. Quadrant percentages show how often the two effects agree vs. disagree in sign. - **Supplementary:** percent of DD genes (significant in both diploid and allotetraploid) with the same effect-size sign, per cell type.

### Cis vs Trans Analysis

Classifies each diploid-DD gene's *cis* proportion (the allotetraploid, same-nucleus effect as a fraction of the total diploid effect) per cell type, using cutoffs *cis* \> 0.7 / *trans* \< 0.3. - Number of *cis*/*trans* DD genes per cell type. - ECDF comparing absolute effect size between *cis* and *trans* genes (K-S test). - Mean *cis* proportion vs. number of cell types a gene is DD in, plus a robustness check that absolute effect size doesn't itself trend with that same number-of-cell-types axis (ruling out an effect-size confound). - **LOEUF** (loss-of-function intolerance) and **percent coding identity** compared between *cis* and *trans* genes, each with a supplementary sweep across six *cis*/*trans* cutoff choices (0.6/0.4 through 0.85/0.15) to check the main 0.7/0.3 result isn't an artifact of that specific cutoff. - **Supplementary: Determining cis/trans cutoffs** -- the analysis behind choosing 0.7/0.3 in the first place: number of DD genes retained at each candidate cutoff (separately for *cis* and *trans*), and how many cell types each gene remains *cis*/*trans* in across the same cutoff sweep. - **Supplementary: Cis/Trans Replication at a Stricter Mean Expression Cutoff (-7)** -- re-derives the main *cis*/*trans* panel (ECDF, mean *cis* proportion, LOEUF, percent coding identity) after dropping genes that only pass the pipeline's -6 cutoff by a narrow margin, to check the main results hold at a stricter -7 cutoff.

### Supplementary: Cell Counts

-   Cell counts per cell type, species, and ploidy context (human, chimp, allotetraploid), excluding the same low-count donor x cell-type combinations dropped upstream by the DD/DE pipeline's 100-cell filter.
-   Cell counts by individual donor.
-   Cell type proportions by individual donor (what fraction of each cell type's cells come from each donor, per species) -- a QC check for whether any single donor disproportionately drives a cell type's result.

### Packages and versions

| R 4.4.1      | Version |
|--------------|---------|
| Seurat       | 5.3.0   |
| SeuratObject | 5.2.0   |
| dplyr        | 1.1.4   |
| tidyverse    | 2.0.0   |
| ggplot2      | 3.5.2   |
| ggrepel      | 0.9.6   |
| Matrix       | 1.7.3   |
| purrr        | 1.1.0   |
| ggpubr       | 0.6.3   |
| stringr      | 1.5.1   |
| scales       | 1.4.0   |
| patchwork    | 1.3.2   |

Run under `module load R/4.4.1 gsl hdf5/1.12.0`, `R_LIBS_USER=/project/gilad/awchen55/Rlibs/R4-4-1`.

------------------------------------------------------------------------

## `q_threshold_sensitivity.Rmd`

**Input:** same master table as above, but loaded once per `q` value from `data/dd_de_results/mt_rb_mean_filtered_q_thresholds/q{0.05,0.1,0.15,0.2}/DD_DE_interspecies_results.csv`. `q` is the assumed UMI capture efficiency passed to `memento`; this document asks how sensitive the diploid DD and recapitulated DD gene calls are to that choice.

### 1. Gene Universe Size per q

Number of genes in the master table at each `q`.

### 2. Diploid DD Gene Sets

Number of diploid DD genes per cell type, at each `q`, and how that count trends with `q`.

### 3. Recapitulated DD Gene Sets

Number of recapitulated DD genes per cell type, at each `q`, and how that count trends with `q`.

### 4. Retention of q = 0.2 DD Genes as q Decreases

For each baseline DD gene set (using `q = 0.2`, the most permissive value, as baseline), tracks what fraction of genes remain significant, become non-significant, or drop out of the gene universe entirely as `q` decreases to 0.05. Repeated with `q = 0.05` as the baseline instead, tracking the same three outcomes plus genes *gained* at more permissive q. Retention/fate plots are faceted by cell type, with wrapped (max two-line) cell type labels so long names stay legible in every panel.

### Packages and versions

| R 4.4.1 | Version |
|---------|---------|
| dplyr   | 1.1.4   |
| tidyr   | 1.3.1   |
| ggplot2 | 3.5.2   |
| purrr   | 1.1.0   |
| tibble  | 3.3.0   |

Same R environment as above (`module load R/4.4.1 gsl hdf5/1.12.0`, `R_LIBS_USER=/project/gilad/awchen55/Rlibs/R4-4-1`).

------------------------------------------------------------------------

## Notes for future edits

-   `results_path` (`/project/gilad/awchen55/mechanism/results/dispersion_paper/`) is where figures are saved via commented-out `ggsave()` calls throughout -- uncomment the relevant line(s) to regenerate a specific figure file rather than re-running the whole document.
-   Several plotting helper functions (`lof_cutoff()`, `pci_cutoff()`, `round_up_signif()`) are defined once and reused across sections with different arguments -- if editing one, check all call sites.
