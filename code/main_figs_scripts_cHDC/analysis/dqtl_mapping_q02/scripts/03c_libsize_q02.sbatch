#!/bin/bash
#SBATCH --job-name=dqtl_03c_libsize_q02
#SBATCH --output=analysis/dqtl_mapping/outputs/covariate_optimization/logs/%x_%j.out
#SBATCH --error=analysis/dqtl_mapping/outputs/covariate_optimization/logs/%x_%j.err
#SBATCH --partition=caslake
#SBATCH --account=pi-gilad
#SBATCH --qos=caslake
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
echo "Working directory: $(pwd)"
mkdir -p analysis/dqtl_mapping/outputs/covariate_optimization/logs

: "${PREP_JOB_ID:?Set PREP_JOB_ID explicitly to the 01_prep_pcs_covopt_q02.sbatch job ID}"

PREP_RECEIPT="analysis/dqtl_mapping/outputs/covariate_optimization/logs/dqtl_prep_pcs_covopt_q02_${PREP_JOB_ID}.run_id"
if [ ! -f "$PREP_RECEIPT" ]; then
  echo "ERROR: no RUN_ID receipt at $PREP_RECEIPT -- was PREP_JOB_ID=$PREP_JOB_ID really a 01_prep_pcs_covopt_q02.sbatch job?" >&2
  exit 1
fi
export DQTL_SOURCE_RUN_ID=$(cat "$PREP_RECEIPT")
echo "PREP_JOB_ID: $PREP_JOB_ID -> DQTL_SOURCE_RUN_ID: $DQTL_SOURCE_RUN_ID"

export DQTL_CELL_TYPES="Activated_fibroblasts,Cardiac_fibroblasts,Cardiomyocytes,Foregut_endodermal_cells"

export DQTL_RUN_ID=$(date +%Y%m%d_%H%M%S)
echo "Final RUN_ID (for steps 04-06b): $DQTL_RUN_ID"
echo "$DQTL_RUN_ID" > "analysis/dqtl_mapping/outputs/covariate_optimization/logs/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.run_id"

module unload node/24.2
module load gcc/10.2.0
module load R/4.4.1

echo "Using Rscript: $(which Rscript)"

Rscript analysis/dqtl_mapping/scripts/03c_build_autopick_libsize_covariates.R

echo "Done. RUN_ID $DQTL_RUN_ID ready for step 04 (04_map_dqtls.R via run_map_autopick_libsize.sbatch):"
echo "  DQTL_RUN_ID=$DQTL_RUN_ID DQTL_CELL_TYPES=$DQTL_CELL_TYPES sbatch analysis/dqtl_mapping/scripts/run_map_autopick_libsize.sbatch"
