#!/bin/bash
#SBATCH --job-name=dd_de_q_multi_mt_rb_mean
#SBATCH --account=pi-gilad
#SBATCH --partition=caslake
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=36:00:00
#SBATCH --output=logs/snakemake_%j.out
#SBATCH --error=logs/snakemake_%j.err

source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
export PATH="/software/python-anaconda-2022.05-el8-x86_64/bin:$PATH"
export R_LIBS_USER=/project/gilad/awchen55/Rlibs/R4-4-1
module load R/4.4.1 gsl hdf5/1.12.0

SNAKEMAKE=/project/gilad/awchen55/envs/snakemake/bin/snakemake

cd /project/gilad/awchen55/mechanism/code/de_dd_snakemake_mt_rb_mean_filtered

$SNAKEMAKE --profile profile \
    --default-resources slurm_partition=caslake slurm_account=pi-gilad mem_mb=8000 cpus_per_task=1 \
    "$@"
