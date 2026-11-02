#!/bin/bash

set -euo pipefail
cd "$(dirname "$0")/../../.."   # repo root
REPO_DIR="$(pwd)"
BASE="figure_reproduction/outputs/basics_scran_crossmethod"

OTHER_CTS="Cardiomyocytes Mesenchymal_cells Cardiac_fibroblasts Activated_fibroblasts Foregut_endodermal_cells"

mkdir -p "$BASE/sce" "$BASE/epsilon" "$BASE/memento_1000subset"

for ct in $OTHER_CTS; do
  ln -sf "$REPO_DIR/analysis/basics_dispersion_reanalysis/outputs/sce/${ct}_sce_raw.rds" \
    "$BASE/sce/${ct}_sce_raw.rds"
  ln -sf "$REPO_DIR/analysis/basics_dispersion_reanalysis/outputs/sce/${ct}_1000subset_sce_raw.rds" \
    "$BASE/sce/${ct}_1000subset_sce_raw.rds"
  ln -sf "$REPO_DIR/analysis/basics_dispersion_reanalysis/outputs/epsilon/${ct}_1000subset_epsilon.csv" \
    "$BASE/epsilon/${ct}_1000subset_epsilon.csv"
  ln -sf "$REPO_DIR/analysis/qtl_analyses/outputs/brendan_cell_count_sweep/cell_count_1000/dispersion_by_ind_${ct}.csv" \
    "$BASE/memento_1000subset/dispersion_by_ind_${ct}.csv"
done

echo "Mirror symlinks created for: $OTHER_CTS"
echo "sce/:"; ls -la "$BASE/sce/"
echo "epsilon/:"; ls -la "$BASE/epsilon/"
echo "memento_1000subset/:"; ls -la "$BASE/memento_1000subset/"
