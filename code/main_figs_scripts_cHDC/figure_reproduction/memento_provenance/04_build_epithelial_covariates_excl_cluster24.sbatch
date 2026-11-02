#!/bin/bash
#SBATCH --job-name=epi_excl24_covariates
#SBATCH --output=figure_reproduction/outputs/logs/%x_%j.out
#SBATCH --error=figure_reproduction/outputs/logs/%x_%j.err
#SBATCH --partition=gilad-hm
#SBATCH --account=pi-gilad
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64G

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
echo "Working directory: $(pwd)"

PY=/project/gilad/brendan/miniforge3/envs/my_jupyter_env/bin/python
"$PY" figure_reproduction/memento_provenance/04_build_epithelial_covariates_excl_cluster24.py

echo "Done."
