#!/bin/bash
#SBATCH --job-name=subcluster_export
#SBATCH --output=analysis/subcluster_dispersion_robustness/outputs/logs/%x_%j.out
#SBATCH --error=analysis/subcluster_dispersion_robustness/outputs/logs/%x_%j.err
#SBATCH --partition=gilad-hm
#SBATCH --account=pi-gilad
#SBATCH --qos=gilad-hm
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=250G

set -euo pipefail
# $0 is unreliable under sbatch (Slurm copies the script to a spool dir before
# running it) -- use $SLURM_SUBMIT_DIR instead (see .living/learnings.md, 2026-07-27).
cd "$SLURM_SUBMIT_DIR"
echo "Working directory: $(pwd)"
mkdir -p analysis/subcluster_dispersion_robustness/outputs/logs

module unload node/24.2
module load gcc/10.2.0
module load R/4.4.1
# SeuratDisk -> hdf5r needs libhdf5_hl.so.310/libhdf5.so.310 at runtime, not
# provided by any module in the standard R/4.4.1 sequence above. Must be
# loaded AFTER gcc/10.2.0 + R/4.4.1 (loading it first pulls in its own default
# gcc/13.2.0 + node/24.2, which then breaks R/4.4.1) -- see .living/learnings.md, 2026-07-30.
module load hdf5/1.14.3

Rscript analysis/subcluster_dispersion_robustness/scripts/01_subcluster_and_export.R

echo "Done."
