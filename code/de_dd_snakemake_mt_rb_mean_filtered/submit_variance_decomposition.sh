#!/bin/bash
# Submit memento_variance_decomposition.py as a SLURM array job across every
# diploid and tetraploid cell type, mirroring the production
# memento_diploid_ct / memento_tetraploid_dd_ct rules' inputs and params
# (q=0.2, min_perc_group=0.95, mean_expr_cutoff=-6, min_cells_per_group=100 --
# see config.yaml). One task per (ploidy, celltype) combination.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="/project/gilad/awchen55/mechanism/data/dd_de_results/mt_rb_mean_filtered_q_thresholds/q0.2/variance_decomposition"
LOG_DIR="$SCRIPT_DIR/logs_variance_decomp"
mkdir -p "$OUT_DIR" "$LOG_DIR"

DIPLOID_DIR="/project/gilad/awchen55/mechanism/data/diploid_crossfilt_data/celltype"
TETRAPLOID_DIR="/project/gilad/awchen55/mechanism/data/tetraploid/celltype"

CELLTYPES_DIPLOID=(
  Bronchiolar_and_alveolar_epithelial_cells Cardiac_mesoderm_cells Cardiomyocytes
  Ciliated_epithelial_cells ENS_glia Epicardial_fat_cells Glioblast Goblet_cells
  Hepatoblasts Intestinal_epithelial_cells Metanephric_cells Neuroblast Neuron
  Neuronal_IPC Radial_glia Retinal_neurons Squamous_epithelial_cells Stromal_cells
  Vascular_endothelial_cells iPSCs
)
# Kept identical to CELLTYPES_DIPLOID (Early_mesoderm_cells removed -- no
# diploid counterpart, matching config.yaml's celltypes_tetraploid).
CELLTYPES_TETRAPLOID=(
  Bronchiolar_and_alveolar_epithelial_cells Cardiac_mesoderm_cells Cardiomyocytes
  Ciliated_epithelial_cells ENS_glia Epicardial_fat_cells
  Glioblast Goblet_cells Hepatoblasts Intestinal_epithelial_cells Metanephric_cells
  Neuroblast Neuron Neuronal_IPC Radial_glia Retinal_neurons Squamous_epithelial_cells
  Stromal_cells Vascular_endothelial_cells iPSCs
)

JOBS_TSV="$LOG_DIR/job_list.tsv"
: > "$JOBS_TSV"
for ct in "${CELLTYPES_DIPLOID[@]}"; do
  printf "diploid\t%s\t%s\n" "$ct" "$DIPLOID_DIR/combined_1kFeatures_${ct}.h5ad" >> "$JOBS_TSV"
done
for ct in "${CELLTYPES_TETRAPLOID[@]}"; do
  printf "tetraploid\t%s\t%s\n" "$ct" "$TETRAPLOID_DIR/tetraploid_combined_variablfeatures_filtered_${ct}.h5ad" >> "$JOBS_TSV"
done
N_JOBS=$(wc -l < "$JOBS_TSV")

sbatch --account=pi-gilad \
       --partition=gilad-hm \
       --job-name=memento_variance_decomp \
       --array=1-"${N_JOBS}" \
       --cpus-per-task=4 \
       --mem=96G \
       --time=02:00:00 \
       --output="${LOG_DIR}/%x_%A_%a.out" \
       --error="${LOG_DIR}/%x_%A_%a.err" \
       --wrap="
set -euo pipefail
LINE=\$(sed -n \"\${SLURM_ARRAY_TASK_ID}p\" '${JOBS_TSV}')
PLOIDY=\$(echo \"\$LINE\" | cut -f1)
CT=\$(echo \"\$LINE\" | cut -f2)
H5AD=\$(echo \"\$LINE\" | cut -f3)

/project/gilad/awchen55/envs/memento/bin/python3 '${SCRIPT_DIR}'/scripts/memento_variance_decomposition.py \
  --h5ad \"\${H5AD}\" \
  --celltype \"\${CT}\" \
  --ploidy \"\${PLOIDY}\" \
  --q 0.2 --min-perc-group 0.95 --mean-expr-cutoff -6 --min-cells-per-group 100 \
  --out-csv-genes '${OUT_DIR}'/\${PLOIDY}_variance_genes_\${CT}.csv \
  --out-csv-summary '${OUT_DIR}'/\${PLOIDY}_variance_summary_\${CT}.csv
"

echo "Submitted ${N_JOBS} array tasks (job list: ${JOBS_TSV})"
