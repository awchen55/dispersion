# supp_figs_scripts_cHDC

Scripts that produce the cHDC paper's supplementary figures (resubmission, assembled 2026-10-09).
Copied from the cHDC_data repo, keeping their repo-relative paths so `source()` calls between scripts still resolve. Hard-coded paths still point at `/project/gilad/brendan/dispersion/pilot/cHDC_data/`. Repoint them to the deposited data before running elsewhere.

## Shared inputs and helpers

- `figure_reproduction/memento_provenance/`: builds `figure_reproduction/data/{res_var,mean_expression}_excl_cluster24.csv` and `aggregated_metadata_excl_cluster24.csv`, the inputs to almost every figure below (Epithelial sub-cluster 24 removed, Memento rerun, spliced into the published matrices).
- `figure_reproduction/scripts/01_exclude_topexpressed_genes.R` and `02_rerun_ranked_dispersion.R`: the -6 retained-gene matrix and per-cell-type dispersion ranking/quintiles used by most figures.
- Helpers sourced by the scripts: `figure_reproduction/scripts/utils.R`, `.../cross_cohort_dispersion_corr/utils.R`, `.../mean_matched_regulatory_revision/config.R`, `.../q_threshold_sensitivity/utils.R`, `analysis/subcluster_dispersion_robustness/scripts/utils.R`, `analysis/dd_sweep_reviewer_response/scripts/00_config.R`, `analysis/epithelial_subcluster24_excl_sensitivity/scripts/utils.R`.

## Upstream inputs NOT included (produced elsewhere)

- The original canonical Memento run (published `new_res_var/mean_expression_by_cell_type_ind_1000_cts_new_predicted_id.csv`) and cell-type annotation (`cellranger_cluster-mode_trial/`).
- Raw-count h5ad / Seurat objects that the Memento-rerunning scripts read.
- BASiCS MCMC chains and 1000-cell Memento subsets for the 5 non-Epithelial cell types (`analysis/basics_dispersion_reanalysis/`, `analysis/qtl_analyses/outputs/brendan_cell_count_sweep/`), symlinked by `basics_scran_crossmethod/setup_mirror_dirs.sh`.
- DD-sweep fits for the dd_gene_count figures (`analysis/qtl_analyses/outputs/dd_sweep/`, Erik McIntire's cohort fits).
- External data: Erik McIntire's MatrixEQTL/mash eQTL results, MSigDB C5 GO:BP GMT, TWN/PPI/regulatory annotation tables, Cell Ranger per-lane metrics (`sequencing_metrics_per_lane.csv`), resolution 0.5/0.8 Memento matrices.

## Figure -> scripts

Each list is the full set of scripts the figure depends on (from a Snakemake `--forceall` dry run of its output, or traced by hand where noted), sorted by path rather than run order. Numbered prefixes within a folder give the order.

### 07_retention_fate_by_contrast.pdf
_Snakemake `qthresh_synthesis`_

- `figure_reproduction/scripts/q_threshold_sensitivity/00_build_source_h5ad_excl_epithelial_cluster24.py`
- `figure_reproduction/scripts/q_threshold_sensitivity/01_memento_q_sweep.py`
- `figure_reproduction/scripts/q_threshold_sensitivity/02_dream_dd_refit_per_q.R`
- `figure_reproduction/scripts/q_threshold_sensitivity/03_dispersion_quintiles_per_q.R`
- `figure_reproduction/scripts/q_threshold_sensitivity/04_q_sensitivity_synthesis.R`
- `figure_reproduction/scripts/q_threshold_sensitivity/utils.R`

### 08_quintile_retention_fate_overall.pdf
_Snakemake `qthresh_quintile_retention_fate`_

- `figure_reproduction/scripts/q_threshold_sensitivity/00_build_source_h5ad_excl_epithelial_cluster24.py`
- `figure_reproduction/scripts/q_threshold_sensitivity/01_memento_q_sweep.py`
- `figure_reproduction/scripts/q_threshold_sensitivity/03_dispersion_quintiles_per_q.R`
- `figure_reproduction/scripts/q_threshold_sensitivity/05_quintile_retention_fate.R`
- `figure_reproduction/scripts/q_threshold_sensitivity/utils.R`

### across_subcluster_cor_heatmap_combined.pdf
_analysis/ only (`mincells700_minind2` config); original matrices by design_

- `analysis/subcluster_dispersion_robustness/scripts/01_subcluster_and_export.R`
- `analysis/subcluster_dispersion_robustness/scripts/01_subcluster_and_export.sbatch`
- `analysis/subcluster_dispersion_robustness/scripts/02_memento_per_subcluster.py`
- `analysis/subcluster_dispersion_robustness/scripts/02_memento_per_subcluster.sbatch`
- `analysis/subcluster_dispersion_robustness/scripts/03_robustness_correlations.R`
- `analysis/subcluster_dispersion_robustness/scripts/05_combined_correlation_heatmaps.R`

### alt_clust_go_af2.pdf (panels: res_05/res_08_consistently_low_disp_GO_enrich.pdf)
_Snakemake `ddx_47`; -6 exclusion applied; resolutions not rerun without sub-cluster 24_

- `figure_reproduction/scripts/47_figS_alt_cluster_go_dotplots.R`

### basics_crossmethod_combined.pdf
_Snakemake `cross_method.smk`_

- `figure_reproduction/scripts/01_exclude_topexpressed_genes.R`
- `figure_reproduction/scripts/02_rerun_ranked_dispersion.R`
- `figure_reproduction/scripts/basics_scran_crossmethod/00_export_epithelial_excl_barcodes.py`
- `figure_reproduction/scripts/basics_scran_crossmethod/01_build_epithelial_sce.R`
- `figure_reproduction/scripts/basics_scran_crossmethod/02_new_epithelial_1000subset_memento.py`
- `figure_reproduction/scripts/basics_scran_crossmethod/03_build_epithelial_1000subset_sce.R`
- `figure_reproduction/scripts/basics_scran_crossmethod/04_run_basics_mcmc_epithelial.R`
- `figure_reproduction/scripts/basics_scran_crossmethod/06_run_scran_all_celltypes_fullcell.R`
- `figure_reproduction/scripts/basics_scran_crossmethod/08_celltype_specificity_basics.R`
- `figure_reproduction/scripts/basics_scran_crossmethod/10_consistent_geneset_retention.R`
- `figure_reproduction/scripts/basics_scran_crossmethod/12_test_de_pairwise_celltypes.R`
- `figure_reproduction/scripts/basics_scran_crossmethod/13_crosscheck_memento_pairwise_dd.R`
- `figure_reproduction/scripts/basics_scran_crossmethod/16_fit_dream_dd_1000subset.R`
- `figure_reproduction/scripts/basics_scran_crossmethod/setup_mirror_dirs.sh`
- `figure_reproduction/scripts/basics_scran_crossmethod/11_basics_crossmethod_combined_figure.R`

### cohort_scatter.pdf, matched_vs_mismatched.pdf
_Snakemake `cross_cohort_correlate` (CCDC_CELL_COUNT_THRESHOLD=1000)_

- `figure_reproduction/scripts/cross_cohort_dispersion_corr/01_load_and_align_dispersion.R`
- `figure_reproduction/scripts/cross_cohort_dispersion_corr/02_correlate_cohorts.R`
- `figure_reproduction/scripts/cross_cohort_dispersion_corr/utils.R`

### combined_loess_no_mt_highlight_grid.pdf
_Snakemake `ddx_21`_

- `figure_reproduction/scripts/01_exclude_topexpressed_genes.R`
- `figure_reproduction/scripts/06_loess_mean_dispersion_diagnostic.R`
- `figure_reproduction/scripts/21_combined_loess_no_mt_highlight_grid.R`

### ct_spec_TWNs_supp.pdf (= FigureS10_composite.pdf)
_run standalone (no rule)_

- `figure_reproduction/scripts/01_exclude_topexpressed_genes.R`
- `figure_reproduction/scripts/43_figS10_twn_cutoff_robustness.R`

### dd_effect_vs_artifact_histograms-2.pdf
_Snakemake `ddx_15`_

- `figure_reproduction/scripts/01_exclude_topexpressed_genes.R`
- `figure_reproduction/scripts/03_rerun_dd_rank_analyses.R`
- `figure_reproduction/scripts/14_dd_trend_corrected_effects.R`
- `figure_reproduction/scripts/15_dd_effect_vs_artifact_histograms.R`

### dd_gene_count_comparison_brendan_vs_erik.pdf
_analysis/ only; run from analysis/dd_sweep_reviewer_response/ (scripts use relative `scripts/00_config.R`); no Epithelial cell types involved_

- `analysis/dd_sweep_reviewer_response/scripts/01_refit_brendan_full_reference_no_random_effect.R`
- `analysis/dd_sweep_reviewer_response/scripts/02_build_dd_gene_count_comparison.R`

### dd_gene_count_erik_by_threshold.pdf
_analysis/ only; as above_

- `analysis/dd_sweep_reviewer_response/scripts/05_erik_dd_gene_count_by_threshold.R`

### excl_top_mods.pdf (= regulatory_network_features_percelltype_module_exclusion_sweep_csore_universe.pdf)
_Snakemake `func_enrich_regulatory_features_plot`_

- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/00_subset_epithelial_excl_cluster24_seurat.R`
- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/00a_extract_epithelial_excl_cluster24_barcodes.py`
- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/01_memento_per_celltype_epithelial.py`
- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/02_prepare_cscore_inputs_epithelial.R`
- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/03_run_cscore_wgcna_epithelial.R`
- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/04_floor_sweep_first_crossing_epithelial.R`
- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/06_consensus_wgcna_epithelial.R`
- `figure_reproduction/scripts/functional_enrichment_percelltype/03_tata_enrichment_percelltype_csore_universe.R`
- `figure_reproduction/scripts/functional_enrichment_percelltype/04_tata_enrichment_module_exclusion_sweep_csore_universe.R`
- `figure_reproduction/scripts/functional_enrichment_percelltype/05_regulatory_network_features_percelltype_module_exclusion_sweep_csore_universe.R`
- `figure_reproduction/scripts/functional_enrichment_percelltype/05b_plot_regulatory_network_features_percelltype_module_exclusion_sweep_csore_universe.R`
- `figure_reproduction/scripts/module_dispersion_profile/02_module_dispersion_profile_consensus.R`

### FigureS8_matrixeqtl_sliding_window.pdf, FigureS9_mash_sliding_window.pdf
_Snakemake `ddx_45`_

- `figure_reproduction/scripts/01_exclude_topexpressed_genes.R`
- `figure_reproduction/scripts/45_figS8_S9_egene_sliding_window.R`

### geneset_robustness_boxplot.pdf
_analysis/ only (`mincells700_minind2`)_

- `analysis/subcluster_dispersion_robustness/scripts/01_subcluster_and_export.R`
- `analysis/subcluster_dispersion_robustness/scripts/01_subcluster_and_export.sbatch`
- `analysis/subcluster_dispersion_robustness/scripts/02_memento_per_subcluster.py`
- `analysis/subcluster_dispersion_robustness/scripts/02_memento_per_subcluster.sbatch`
- `analysis/subcluster_dispersion_robustness/scripts/03_robustness_correlations.R`

### hl_shift_by_resolution.pdf
_run standalone (no rule)_

- `figure_reproduction/scripts/01_exclude_topexpressed_genes.R`
- `figure_reproduction/scripts/02_rerun_ranked_dispersion.R`
- `analysis/subcluster_dispersion_robustness/scripts/01_subcluster_and_export.R`
- `analysis/subcluster_dispersion_robustness/scripts/01_subcluster_and_export.sbatch`
- `analysis/subcluster_dispersion_robustness/scripts/02_memento_per_subcluster.py`
- `analysis/subcluster_dispersion_robustness/scripts/02_memento_per_subcluster.sbatch`
- `figure_reproduction/scripts/44_subcluster_retained_genesets_hl_by_resolution.R`

### mean_matched_dists.pdf (= balance_check_barplot.pdf)
_Snakemake `mm_mean_match` (script 02 writes the barplot)_

- `figure_reproduction/scripts/01_exclude_topexpressed_genes.R`
- `figure_reproduction/scripts/mean_matched_regulatory_revision/01_build_dispersion_gene_sets.R`
- `figure_reproduction/scripts/mean_matched_regulatory_revision/02_mean_match_gene_sets.R`
- `figure_reproduction/scripts/mean_matched_regulatory_revision/config.R`

### mean_res_var_heatmaps.pdf (= FigureS6_composite.pdf)
_Snakemake `ddx_41`_

- `figure_reproduction/scripts/01_exclude_topexpressed_genes.R`
- `figure_reproduction/scripts/41_figS6_mean_and_mean_vs_dispersion_heatmaps.R`

### module_plots.pdf (= consensus_qc_combined.pdf)
_Snakemake `consensus_qc_combined_figure`_

- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/00_subset_epithelial_excl_cluster24_seurat.R`
- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/00a_extract_epithelial_excl_cluster24_barcodes.py`
- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/02_prepare_cscore_inputs_epithelial.R`
- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/03_run_cscore_wgcna_epithelial.R`
- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/04_floor_sweep_first_crossing_epithelial.R`
- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/05a_effect_size_sweep.R`
- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/06_consensus_wgcna_epithelial.R`
- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/08_consensus_tom_and_cohesion.R`
- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/10_consensus_module_count_plot.R`
- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/11_consensus_qc_combined_figure.R`
- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/12_consensus_module_dendrograms.R`

### permut_cellid_CorrelationHeatmap.pdf
_Snakemake `cp_01`, `cp_02`_

- `figure_reproduction/scripts/celltype_permutation/01_permute_memento.py`
- `figure_reproduction/scripts/celltype_permutation/02_figS7_permuted_heatmap.R`
- `figure_reproduction/scripts/q_threshold_sensitivity/00_build_source_h5ad_excl_epithelial_cluster24.py`

### QC_plots_af2.pdf (panels: <metric>_barplot.pdf)
_Snakemake `ddx_46`_

- `figure_reproduction/scripts/46_figS_qc_sequencing_barplots.R`

### ranked_mods.pdf (= module_dispersion_profile_consensus_combined_positive_scale_multipanel.pdf)
_Snakemake `module_dispersion_profile_multipanel`_

- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/00_subset_epithelial_excl_cluster24_seurat.R`
- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/00a_extract_epithelial_excl_cluster24_barcodes.py`
- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/01_memento_per_celltype_epithelial.py`
- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/02_prepare_cscore_inputs_epithelial.R`
- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/03_run_cscore_wgcna_epithelial.R`
- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/04_floor_sweep_first_crossing_epithelial.R`
- `figure_reproduction/scripts/coexpression_plasticity_lowfidelity/06_consensus_wgcna_epithelial.R`
- `figure_reproduction/scripts/module_dispersion_profile/02_module_dispersion_profile_consensus.R`
- `figure_reproduction/scripts/module_dispersion_profile/07_module_dispersion_profile_markers_consensus.R`
- `figure_reproduction/scripts/module_dispersion_profile/09_module_dispersion_profile_consensus_positive_scale.R`
- `figure_reproduction/scripts/module_dispersion_profile/12_cell_cycle_module_phase_composition_consensus.R`
- `figure_reproduction/scripts/module_dispersion_profile/14_module_dispersion_profile_consensus_combined.R`
- `figure_reproduction/scripts/module_dispersion_profile/16_module_dispersion_profile_consensus_combined_positive_scale_multipanel.R`

### regulatory_network_composite_excl_cellstate.pdf
_Snakemake `rnx_02_composite`_

- `figure_reproduction/scripts/01_exclude_topexpressed_genes.R`
- `figure_reproduction/scripts/02_rerun_ranked_dispersion.R`
- `figure_reproduction/scripts/03_rerun_dd_rank_analyses.R`
- `figure_reproduction/scripts/04_replicate_fig2_panels.R`
- `figure_reproduction/scripts/27_fig3_panels_excl_topexpressed.R`
- `figure_reproduction/scripts/29_fig4_panels_excl_topexpressed.R`
- `figure_reproduction/scripts/39_dd_refit_excl_cellstate_genes.R`
- `figure_reproduction/scripts/regulatory_network_excl_cellstate/01_panels_excl_cellstate.R`
- `figure_reproduction/scripts/regulatory_network_excl_cellstate/02_composite_excl_cellstate.R`
- `figure_reproduction/scripts/utils.R`

### regulatory_network_composite_matched.pdf
_Snakemake `mm_09_composite`_

- `figure_reproduction/scripts/01_exclude_topexpressed_genes.R`
- `figure_reproduction/scripts/03_rerun_dd_rank_analyses.R`
- `figure_reproduction/scripts/04_replicate_fig2_panels.R`
- `figure_reproduction/scripts/mean_matched_regulatory_revision/01_build_dispersion_gene_sets.R`
- `figure_reproduction/scripts/mean_matched_regulatory_revision/02_mean_match_gene_sets.R`
- `figure_reproduction/scripts/mean_matched_regulatory_revision/05_fig3b_3d_4a_4b_matched.R`
- `figure_reproduction/scripts/mean_matched_regulatory_revision/06_fig3c_percelltype_matched.R`
- `figure_reproduction/scripts/mean_matched_regulatory_revision/08_remaining_panels_matched.R`
- `figure_reproduction/scripts/mean_matched_regulatory_revision/09_composite_matched.R`
- `figure_reproduction/scripts/mean_matched_regulatory_revision/config.R`
- `figure_reproduction/scripts/utils.R`

### subcluster_vs_original_rho_combined.pdf
_analysis/ only (`mincells700_minind2`)_

- `analysis/subcluster_dispersion_robustness/scripts/01_subcluster_and_export.R`
- `analysis/subcluster_dispersion_robustness/scripts/01_subcluster_and_export.sbatch`
- `analysis/subcluster_dispersion_robustness/scripts/02_memento_per_subcluster.py`
- `analysis/subcluster_dispersion_robustness/scripts/02_memento_per_subcluster.sbatch`
- `analysis/subcluster_dispersion_robustness/scripts/03_robustness_correlations.R`
- `analysis/subcluster_dispersion_robustness/scripts/04_combined_subcluster_vs_original_rho.R`
