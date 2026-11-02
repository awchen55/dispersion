#!/bin/bash
#SBATCH --job-name=dqtl_run_mash
#SBATCH --output=analysis/dqtl_mapping/outputs/mash/logs/%x_%j.out
#SBATCH --error=analysis/dqtl_mapping/outputs/mash/logs/%x_%j.err
#SBATCH --partition=caslake
#SBATCH --account=pi-gilad
#SBATCH --qos=caslake
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
echo "Working directory: $(pwd)"
mkdir -p analysis/dqtl_mapping/outputs/mash/logs

: "${DQTL_RUN_ID:?Set DQTL_RUN_ID explicitly to the RUN_ID 06a_run_fastqtl2mash.sbatch just processed}"
export DQTL_RUN_ID
echo "DQTL_RUN_ID: $DQTL_RUN_ID"

module unload node/24.2
module load gcc/10.2.0
module load R/4.4.1
module load gsl

echo "Using Rscript: $(which Rscript)"

Rscript analysis/dqtl_mapping/scripts/06b_run_mash.R

echo "Done."
