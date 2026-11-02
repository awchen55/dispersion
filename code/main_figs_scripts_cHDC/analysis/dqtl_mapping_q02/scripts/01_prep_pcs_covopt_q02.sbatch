#!/bin/bash
#SBATCH --job-name=dqtl_prep_pcs_covopt_q02
#SBATCH --output=analysis/dqtl_mapping/outputs/covariate_optimization/logs/%x_%j.out
#SBATCH --error=analysis/dqtl_mapping/outputs/covariate_optimization/logs/%x_%j.err
#SBATCH --partition=gilad-hm
#SBATCH --account=pi-gilad
#SBATCH --qos=gilad-hm
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=22
#SBATCH --mem=400G

set -euo pipefail
# $0 is unreliable here -- Slurm copies the submitted script to a spool
# directory before executing it (see .living/learnings.md, 2026-07-27). Use
# $SLURM_SUBMIT_DIR instead.
cd "$SLURM_SUBMIT_DIR"
echo "Working directory: $(pwd)"
mkdir -p analysis/dqtl_mapping/outputs/covariate_optimization/logs

export DQTL_CELL_COUNT_THRESHOLD=150
export DQTL_SD_CUTOFF_FRAC=0.30
export DQTL_CELL_TYPES="Activated_fibroblasts,Cardiac_fibroblasts,Cardiomyocytes,Foregut_endodermal_cells"
export DQTL_GENE_SELECTION_STRATEGY=union

export DQTL_RUN_ID=$(date +%Y%m%d_%H%M%S)
echo "RUN_ID: $DQTL_RUN_ID"
echo "$DQTL_RUN_ID" > "analysis/dqtl_mapping/outputs/covariate_optimization/logs/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.run_id"

echo "Resolved config: cell_count_threshold=$DQTL_CELL_COUNT_THRESHOLD sd_cutoff_frac=$DQTL_SD_CUTOFF_FRAC gene_selection_strategy=$DQTL_GENE_SELECTION_STRATEGY cell_types=$DQTL_CELL_TYPES"

module unload node/24.2
module load gcc/10.2.0
module load R/4.4.1

echo "Using Rscript: $(which Rscript)"

echo "[0/3] Reformatting q=0.2 Memento output..."
Rscript analysis/dqtl_mapping_q02/scripts/00_reformat_memento_q02.R

echo "[1/3] Preparing dispersion inputs (q=0.2)..."
Rscript analysis/dqtl_mapping_q02/scripts/01_prepare_dispersion_inputs_q02.R

echo "[2/3] Calculating dispersion PCs..."
Rscript analysis/dqtl_mapping/scripts/02_calculate_dispersion_pcs.R

echo "[3/3] Optimizing covariates..."
Rscript analysis/dqtl_mapping/scripts/03_optimize_covariates.R

echo "Done. RUN_ID $DQTL_RUN_ID ready for 03c_libsize_q02.sbatch."
