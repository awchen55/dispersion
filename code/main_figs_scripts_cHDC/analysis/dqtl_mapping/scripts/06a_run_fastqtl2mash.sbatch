#!/bin/bash
#SBATCH --job-name=dqtl_fastqtl2mash
#SBATCH --output=analysis/dqtl_mapping/outputs/reformatted_dqtls/logs/fastqtl2mash_%j.out
#SBATCH --error=analysis/dqtl_mapping/outputs/reformatted_dqtls/logs/fastqtl2mash_%j.err
#SBATCH --partition=caslake
#SBATCH --account=pi-gilad
#SBATCH --qos=caslake
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G

set -euo pipefail

REFORMAT_DIR=/project/gilad/brendan/dispersion/pilot/cHDC_data/analysis/dqtl_mapping/outputs/reformatted_dqtls
SIF=/project/gilad/emcintire/chemo_project_2023/pipeline/reformat_eqtls/hdf5tools.sif

: "${DQTL_RUN_ID:?Set DQTL_RUN_ID explicitly to the RUN_ID script 05 just reformatted}"
RUN_ID="$DQTL_RUN_ID"
RUN_DIR="$REFORMAT_DIR/$RUN_ID"

echo "fastqtl2mash: preparing isolated workspace for RUN_ID=$RUN_ID"

# Fail loudly if the run or its script-05 outputs are missing.
[ -d "$RUN_DIR" ]            || { echo "ERROR: run dir not found: $RUN_DIR" >&2; exit 1; }
[ -d "$RUN_DIR/list_files" ] || { echo "ERROR: no list_files/ in $RUN_DIR" >&2; exit 1; }
ls "$RUN_DIR"/dQTL_*_AllChroms.out.gz >/dev/null 2>&1 \
  || { echo "ERROR: no dQTL_*_AllChroms.out.gz in $RUN_DIR" >&2; exit 1; }
[ -f "$REFORMAT_DIR/fastqtl_to_mash.ipynb" ] \
  || { echo "ERROR: shared fastqtl_to_mash.ipynb not found: $REFORMAT_DIR/fastqtl_to_mash.ipynb" >&2; exit 1; }
[ -f "$REFORMAT_DIR/run_fastqtl2mash_LOCAL.sh" ] \
  || { echo "ERROR: shared run_fastqtl2mash_LOCAL.sh not found: $REFORMAT_DIR/run_fastqtl2mash_LOCAL.sh" >&2; exit 1; }

# ---- Make $RUN_DIR self-contained for apptainer (idempotent: safe to re-run) ----
mkdir -p "$RUN_DIR/ready_for_mash"
cp -f "$REFORMAT_DIR/fastqtl_to_mash.ipynb" "$RUN_DIR/fastqtl_to_mash.ipynb"
cp -f "$REFORMAT_DIR/run_fastqtl2mash_LOCAL.sh" "$RUN_DIR/run_fastqtl2mash_LOCAL.sh"
echo "fastqtl2mash: workspace ready at $RUN_DIR"
# -----------------------------------------------------------------------------

module load apptainer

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"

apptainer exec --cleanenv \
  -B "$RUN_DIR:/work" \
  -B "/scratch/midway3/$USER:/scratch" \
  --pwd /work "$SIF" \
  bash run_fastqtl2mash_LOCAL.sh
