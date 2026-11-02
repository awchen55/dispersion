#!/bin/bash
#SBATCH --job-name=epi_excl24_memento
#SBATCH --output=analysis/epithelial_subcluster24_excl_sensitivity/outputs/logs/%x_%j.out
#SBATCH --error=analysis/epithelial_subcluster24_excl_sensitivity/outputs/logs/%x_%j.err
#SBATCH --partition=caslake
#SBATCH --account=pi-gilad
#SBATCH --qos=caslake
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=80G

set -euo pipefail
# $0 is unreliable under sbatch (Slurm copies the script to a spool dir before
# running it) -- use $SLURM_SUBMIT_DIR instead (see .living/learnings.md, 2026-07-27).
cd "$SLURM_SUBMIT_DIR"
echo "Working directory: $(pwd)"
mkdir -p analysis/epithelial_subcluster24_excl_sensitivity/outputs/logs

PY=/project/gilad/brendan/miniforge3/envs/my_jupyter_env/bin/python

"$PY" analysis/epithelial_subcluster24_excl_sensitivity/scripts/01_build_excl_cluster24_epithelial_h5ad.py
"$PY" analysis/epithelial_subcluster24_excl_sensitivity/scripts/02_memento_epithelial_excl_cluster24.py

echo "Done."
