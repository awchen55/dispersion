# main_figs_scripts_cHDC

Scripts that produce the cHDC paper's main figures (resubmission, assembled 2026-10-09).
Copied from the cHDC_data repo, keeping their repo-relative paths so `source()` calls between scripts still resolve. Hard-coded paths still point at `/project/gilad/brendan/dispersion/pilot/cHDC_data/`. Repoint them to the deposited data before running elsewhere.

Figures 1D-4 and 6 come from `figure_reproduction/` (Epithelial_cells sub-cluster 24 excluded; per-cell-type > -6 mean exclusion). The exact run order is fixed by its Snakemake DAG (`figure_reproduction/Snakefile`, not copied); numbered prefixes give the order within a folder.

## Shared inputs and helpers

- `figure_reproduction/memento_provenance/`: builds `figure_reproduction/data/{res_var,mean_expression}_excl_cluster24.csv` and `aggregated_metadata_excl_cluster24.csv` (Epithelial sub-cluster 24 removed, Memento rerun, spliced into the published matrices).
- `figure_reproduction/scripts/01_exclude_topexpressed_genes.R`, `02_rerun_ranked_dispersion.R`, `03_rerun_dd_rank_analyses.R`, `04_replicate_fig2_panels.R`: retained genes, dispersion quintiles, DD/DE dream fits used by Figs 2-4.
- Helpers: `figure_reproduction/scripts/utils.R`, `analysis/epithelial_subcluster24_excl_sensitivity/scripts/utils.R`, `analysis/dqtl_genotype_dosage_plots/scripts/utils.R`.

## Upstream inputs NOT included

- Cell Ranger / demultiplexing / annotation that produced the Seurat object and the canonical Memento matrices (`cellranger_cluster-mode_trial/`).
- Genotypes and other dQTL-mapping inputs (`analysis/dqtl_mapping/inputs/`), the q = 0.2 Memento h5ad, and `/project/gilad/brendan/scripts/conda_init.sh` (sourced by the sbatch files).
- External annotation tables (TWN, PPI, TATA/TSS/enhancer/intron features, eGenes, MSigDB GMT).

## Figure -> scripts

### Figure 1A
_Schematic; no code._

### Figure 1B (UMAP) and 1C (cell counts)
_`figure_reproduction` rule `ddx_48` (script 48, ported from `main_figs_1-4_related_supp.Rmd` Figure 1 '## B'/'## C'). Reads the full Seurat object; Fig 1C counts include Epithelial sub-cluster 24, as published._

- `figure_reproduction/scripts/48_fig1bc_umap_and_cell_counts.R`

### Figure 1D (dispersion vs mean, Cardiomyocytes)
_rule `hc_24` (script 24 output is the panel itself)._

- `figure_reproduction/scripts/01_exclude_topexpressed_genes.R`
- `figure_reproduction/scripts/24_fig1b_cardiomyocytes_excl_topexpressed.R`

### Figure 2 (composite)
_rule `hc_38_fig2_composite`._

- `figure_reproduction/scripts/01_exclude_topexpressed_genes.R`
- `figure_reproduction/scripts/02_rerun_ranked_dispersion.R`
- `figure_reproduction/scripts/03_rerun_dd_rank_analyses.R`
- `figure_reproduction/scripts/04_replicate_fig2_panels.R`
- `figure_reproduction/scripts/05_replot_fig2_panels.R`
- `figure_reproduction/scripts/25_fig2b_varpart_excl_topexpressed.R`
- `figure_reproduction/scripts/32_fig2e_go_dag_semantic_clustering.R`
- `figure_reproduction/scripts/34_fig2e_godag_top10_representatives_plot.R`
- `figure_reproduction/scripts/35_fig2a_narrow_wrapped_plot.R`
- `figure_reproduction/scripts/36_fig2b_varpart_relabeled.R`
- `figure_reproduction/scripts/37_fig2c_heatmap_relabeled.R`
- `figure_reproduction/scripts/38_fig2_composite.R`
- `figure_reproduction/scripts/utils.R`

### Figure 3 (composite)
_rule `hc_28_fig3_composite`._

- `figure_reproduction/scripts/01_exclude_topexpressed_genes.R`
- `figure_reproduction/scripts/02_rerun_ranked_dispersion.R`
- `figure_reproduction/scripts/03_rerun_dd_rank_analyses.R`
- `figure_reproduction/scripts/27_fig3_panels_excl_topexpressed.R`
- `figure_reproduction/scripts/28_fig3_composite_excl_topexpressed.R`
- `figure_reproduction/scripts/utils.R`

### Figure 4 (composite)
_rule `hc_30_fig4_composite`._

- `figure_reproduction/scripts/01_exclude_topexpressed_genes.R`
- `figure_reproduction/scripts/02_rerun_ranked_dispersion.R`
- `figure_reproduction/scripts/03_rerun_dd_rank_analyses.R`
- `figure_reproduction/scripts/04_replicate_fig2_panels.R`
- `figure_reproduction/scripts/29_fig4_panels_excl_topexpressed.R`
- `figure_reproduction/scripts/30_fig4_composite_excl_topexpressed.R`
- `figure_reproduction/scripts/utils.R`

### Figure 5
_Alexander Chen's interspecies figure; not in this folder._

### Figure 6 (dQTL QQ plots + HADHB dosage)
_rule `hc_42_fig6_dqtl_composite` (script 42). Upstream: the q = 0.2 dQTL mapping run RUN_ID 20260925_143116 (`analysis/dqtl_mapping_q02/` + the unmodified `analysis/dqtl_mapping/` steps it calls; order in `analysis/dqtl_mapping_q02/DQTL_MAPPING_Q02.md`), listed in run order._

- `analysis/dqtl_mapping_q02/scripts/01_memento_q02_threshold150.py`
- `analysis/dqtl_mapping_q02/scripts/01_memento_q02_threshold150.sbatch`
- `analysis/dqtl_mapping_q02/scripts/00_reformat_memento_q02.R`
- `analysis/dqtl_mapping_q02/scripts/01_prepare_dispersion_inputs_q02.R`
- `analysis/dqtl_mapping_q02/scripts/01_prep_pcs_covopt_q02.sbatch`
- `analysis/dqtl_mapping/scripts/02_calculate_dispersion_pcs.R`
- `analysis/dqtl_mapping/scripts/03_optimize_covariates.R`
- `analysis/dqtl_mapping_q02/scripts/03c_libsize_q02.sbatch`
- `analysis/dqtl_mapping/scripts/03c_build_autopick_libsize_covariates.R`
- `analysis/dqtl_mapping/scripts/run_map_autopick_libsize.sbatch`
- `analysis/dqtl_mapping/scripts/04_map_dqtls.R`
- `analysis/dqtl_mapping/scripts/run_reformat_only.sbatch`
- `analysis/dqtl_mapping/scripts/05_reformat_for_mash.R`
- `analysis/dqtl_mapping/scripts/06a_run_fastqtl2mash.sbatch`
- `analysis/dqtl_mapping/scripts/run_mash.sbatch`
- `analysis/dqtl_mapping/scripts/06b_run_mash.R`
- `analysis/dqtl_mapping_q02/scripts/04_dosage_boxplots.R`
- `analysis/dqtl_genotype_dosage_plots/scripts/utils.R`
- `figure_reproduction/scripts/42_fig6_dqtl_composite.R`
