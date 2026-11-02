#!/bin/bash
#SBATCH --job-name=dqtl_reformat_only
#SBATCH --output=analysis/dqtl_mapping/outputs/covariate_optimization/logs/%x_%j.out
#SBATCH --error=analysis/dqtl_mapping/outputs/covariate_optimization/logs/%x_%j.err
#SBATCH --partition=caslake
#SBATCH --account=pi-gilad
#SBATCH --qos=caslake
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
echo "Working directory: $(pwd)"

: "${DQTL_RUN_ID:?Set DQTL_RUN_ID explicitly}"
echo "RUN_ID: $DQTL_RUN_ID"

module unload node/24.2
module load gcc/10.2.0
module load R/4.4.1
module load gsl

echo "Using Rscript: $(which Rscript)"

Rscript analysis/dqtl_mapping/scripts/05_reformat_for_mash.R

echo "Done."
