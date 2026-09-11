# Empirical Diagnostics: Dispersion

The purpose of this study is to further characterize **total variability** and **dispersion** as memento estimates them, and to empirically diagnose how reliable those estimates are. Total variability is a gene's raw (mean-uncorrected) variance; dispersion is the residual left after regressing out the mean-variance trend across genes, i.e. the portion of a gene's variability not explained by its mean expression level alone. This document works through: how memento actually computes those two quantities, how much of a gene's total variability its dispersion estimate captures, whether the dispersion estimate still depends on mean expression after that correction, how reproducible it is across individuals/donors, and how sensitive it is to per-sample cell count.

> **Note:** the sentence introducing this README was cut off in the request that produced it ("...the purpose of this study was to further define total variability and dispersion, where"). The paragraph above is my best reconstruction from the document's own content (see the "Total variability decomposition" section) — please correct it if it doesn't match what you intended.

## `dispersion_diagnostics.Rmd`

Two datasets are used throughout:

- **cHDC cardiac pilot** (`/project/gilad/brendan/dispersion/pilot/cHDC_data/...`): 6 cell types x 3 individuals (NA19099, NA19093, NA18522). Used for most sections.
- **mechanism project** diploid (human/chimp) + allotetraploid memento output (`data/dd_de_results/mt_rb_mean_filtered_q_thresholds/q0.2/`): used only in the sections that say so explicitly (3b, 5b).

Three filters are used throughout wherever gene-level filtering is applied: mitochondrial/ribosomal gene removal, a log-mean-expression cutoff of -6, and (for the mechanism-project sections) a minimum of 100 cells per donor x cell-type group.

### 1. Total Variability Decomposition

What fraction of a gene's total (raw, mean-uncorrected) variance is "dispersion" vs. explained by its mean expression level, computed per gene x cell type x ploidy x species x donor by `memento_variance_decomposition.py`.

**How memento estimates mean and dispersion.** Memento estimates each gene's mean and variance with a hypergeometric method-of-moments estimator (Kim et al. 2024, PMC11556465, Equation 8), modeling observed UMI counts as a hypergeometric sample of a cell's true transcript pool:

$$\hat{\mu}_g = \frac{1}{n_{cells}} \sum_c \frac{Y_{cg}}{N_c \, q}$$

$$\hat{\sigma}_g^2 = \frac{1}{n_{cells}} \sum_c \frac{Y_{cg}^2 - Y_{cg}(1 - q)}{N_c^2 \, q^2} - \hat{\mu}_g^2$$

where, for gene $g$:

- $Y_{cg}$ -- the observed (sequenced) UMI count for gene $g$ in cell $c$.
- $N_c$ -- the total UMI count of cell $c$.
- $q$ -- the capture efficiency: the assumed probability a true transcript is captured and sequenced (`memento_q` in the pipeline config, the same value for every cell in a run).
- $n_{cells}$ -- the number of cells in the group being estimated (e.g. one species x donor x cell-type sample).
- $\hat{\mu}_g$ and $\hat{\sigma}_g^2$ -- the resulting mean and variance estimates.

The $Y_{cg}(1-q)$ term corrects the raw second moment for the extra sampling variance introduced purely by imperfect capture, so $\hat{\sigma}_g^2$ isn't a capture-noise-inflated estimate of the gene's true variance.

Memento's actual code (`_hyper_1d_relative`) computes an algebraically equivalent quantity in *relative* units -- without the $q$ / $q^2$ divisions -- which is an exact constant rescaling that leaves everything below unchanged; see the Rmd for the derivation.

**Total variability vs. dispersion.** Three quantities per gene:

-   **Total variability**, $V = \hat{\sigma}_g^2$ -- the gene's raw, mean-uncorrected variance (Memento's internal second moment).

-   **Mean component of variance**, $E$ -- the variance this gene would have if it behaved exactly like a typical gene at its expression level. Across all genes in a sample, memento fits a mean-variance trend:

    $$\log(V) \sim \text{poly}(\log(\hat{\mu}_g), \text{degree}=2)$$

    ($\log(E)$ is that trend's fitted value at gene $g$'s mean, computed by `_fit_mv_regressor`) -- the mean-driven baseline that the gene's actual variance gets compared against.

-   **Dispersion** -- simply the log-ratio of the two above:

    $$r_g = \log(\text{dispersion}) = \log(V) - \log(E) = \log(V / E)$$

    $r_g = 0$ means the gene sits exactly on the mean-variance trend (the linear ratio $\exp(r_g) = V/E$ would be 1).

These satisfy the exact identity $\log(V) = \log(E) + r_g$ (i.e. `log.total.variability = log.mean.component.of.variance + log.dispersion`). 


### 2. Mean Dependence

Whether the raw (unfiltered) dispersion estimate still correlates with mean expression -- memento corrects for this internally, so little residual association is expected.

### 3. LOESS Mean Dependence After Expression Cutoff

Repeats the mean-dependence check restricted to genes passing the -6 log-mean cutoff with mito/ribo genes removed, with a per-panel LOESS fit (cell type x individual).

### 3b. Same Diagnostic, Mechanism Project

Repeats the mean-dependence LOESS diagnostic on the mechanism project's human/chimp diploid and allotetraploid data instead of the cHDC pilot data, including the 100-cell-per-donor filter and why it's needed (a low-cell-count numerical artifact, also documented in `q_threshold_sensitivity.Rmd` Section 8b).

### 4. Reproducibility Across Individuals

Pairwise Spearman correlation of dispersion estimates between individuals within each cell type, plus a scatter for one individual pair per cell type.

### 5. Genes Excluded by Expression Threshold

How many genes per cell type are excluded by the -6 log-mean cutoff.

### 5b. DD Genes Removed at Stricter Expression Cutoffs (Mechanism Project)

How many of the mechanism project's current DD genes would be lost if the log-mean cutoff were tightened from -6 to -6.5 or -7.

### 6. Downsampling Reliability Sweep (Cardiomyocytes)

Compares dispersion estimates from downsampled Cardiomyocyte data (100-5,000 cells per individual, 30 replicates per depth) against the full-data reference, to isolate the effect of cell count on estimate reliability from between-sample biological variation. See `code/empirical_diagnostics_dispersion/README.md` and `data/empirical_diagnostics_dispersion/downsample_reliability_cardiomyocytes/README.md` for how this data was generated.

## Packages and versions

| R 4.4.1   | Version |
|-----------|---------|
| dplyr     | 1.1.4   |
| tidyr     | 1.3.1   |
| ggplot2   | 3.5.2   |
| ggrepel   | 0.9.6   |
| purrr     | 1.1.0   |
| tibble    | 3.3.0   |
| stringr   | 1.5.1   |
| gridExtra | 2.3     |

Run under `module load R/4.4.1 gsl hdf5/1.12.0`, `R_LIBS_USER=/project/gilad/awchen55/Rlibs/R4-4-1`.
