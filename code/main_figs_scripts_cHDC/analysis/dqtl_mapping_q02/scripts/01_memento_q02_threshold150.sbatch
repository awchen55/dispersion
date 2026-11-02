#!/bin/bash
#SBATCH --job-name=memento_q02_t150
#SBATCH --output=/project/gilad/brendan/dispersion/pilot/cHDC_data/slurm_output/%x_%A.out
#SBATCH --error=/project/gilad/brendan/dispersion/pilot/cHDC_data/slurm_output/%x_%A.err
#SBATCH --partition=gilad-hm
#SBATCH --account=pi-gilad
#SBATCH --mem=200G
#SBATCH --time=12:00:00

source /project/gilad/brendan/scripts/conda_init.sh
mamba activate /project/gilad/brendan/miniforge3/envs/my_jupyter_env

cd /project/gilad/brendan/dispersion/pilot/cHDC_data
python analysis/dqtl_mapping_q02/scripts/01_memento_q02_threshold150.py
