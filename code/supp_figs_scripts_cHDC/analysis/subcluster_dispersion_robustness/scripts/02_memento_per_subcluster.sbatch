#!/bin/bash
#SBATCH --job-name=subcluster_memento
#SBATCH --output=analysis/subcluster_dispersion_robustness/outputs/logs/%x_%A_%a.out
#SBATCH --error=analysis/subcluster_dispersion_robustness/outputs/logs/%x_%A_%a.err
#SBATCH --partition=caslake
#SBATCH --account=pi-gilad
#SBATCH --qos=caslake
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --array=0-1

set -euo pipefail
# $0 is unreliable under sbatch (Slurm copies the script to a spool dir before
# running it) -- use $SLURM_SUBMIT_DIR instead (see .living/learnings.md, 2026-07-27).
cd "$SLURM_SUBMIT_DIR"
echo "Working directory: $(pwd)"
mkdir -p analysis/subcluster_dispersion_robustness/outputs/logs

# Array index -> (MIN_CELLS_PER_GROUP, MIN_INDIVIDUALS_PER_SUBCLUSTER).
MIN_CELLS_BY_TASK=(1000 700)
MIN_IND_BY_TASK=(2 2)
export SUBCLUSTER_MIN_CELLS_PER_GROUP="${MIN_CELLS_BY_TASK[$SLURM_ARRAY_TASK_ID]}"
export SUBCLUSTER_MIN_INDIVIDUALS="${MIN_IND_BY_TASK[$SLURM_ARRAY_TASK_ID]}"
echo "Array task $SLURM_ARRAY_TASK_ID: SUBCLUSTER_MIN_CELLS_PER_GROUP=$SUBCLUSTER_MIN_CELLS_PER_GROUP SUBCLUSTER_MIN_INDIVIDUALS=$SUBCLUSTER_MIN_INDIVIDUALS"

/project/gilad/brendan/miniforge3/envs/my_jupyter_env/bin/python \
  analysis/subcluster_dispersion_robustness/scripts/02_memento_per_subcluster.py

echo "Done."
