#!/bin/bash
#SBATCH --job-name=dqtl_map_autopick_libsize
#SBATCH --output=analysis/dqtl_mapping/outputs/covariate_optimization/logs/%x_%j.out
#SBATCH --error=analysis/dqtl_mapping/outputs/covariate_optimization/logs/%x_%j.err
#SBATCH --partition=caslake
#SBATCH --account=pi-gilad
#SBATCH --qos=caslake
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=22
#SBATCH --mem=150G

set -euo pipefail
# $0 is unreliable here -- Slurm copies the submitted script to a spool
# directory before executing it, so a $0-relative cd resolves against
# /var/spool/..., not the real submission directory (see .living/learnings.md,
# 2026-07-27). Use $SLURM_SUBMIT_DIR (set by Slurm to the directory sbatch was
# invoked from) instead.
cd "$SLURM_SUBMIT_DIR"
echo "Working directory: $(pwd)"
mkdir -p analysis/dqtl_mapping/outputs/covariate_optimization/logs

: "${DQTL_RUN_ID:?Set DQTL_RUN_ID explicitly}"
: "${DQTL_CELL_TYPES:?Set DQTL_CELL_TYPES explicitly (no default -- 04_map_dqtls.R defaults to the original 3-cell-type set, wrong for these runs)}"

module unload node/24.2
module load gcc/10.2.0
module load R/4.4.1

echo "Using Rscript: $(which Rscript)"
echo "RUN_ID: $DQTL_RUN_ID"
echo "CELL_TYPES: $DQTL_CELL_TYPES"

Rscript analysis/dqtl_mapping/scripts/04_map_dqtls.R

echo "Done."
