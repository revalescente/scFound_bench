#!/bin/bash
#SBATCH --job-name=test_qc_single
#SBATCH --time=00:30:00
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --output=/home/%u/scFound_bench/outputs/test_single-%j.out
#SBATCH --error=/home/%u/scFound_bench/logs/test_single-%j.err

set -euo pipefail

module purge
module load micromamba/bioc-3.23

# File singolo da testare
H5AD_FILE="/projects/shared/cellxgene_split_by_dataset/blood/01209dce-3575-4bed-b1df-129f57fbc031.h5ad"
R_SCRIPT="/home/vreffo/scFound_bench/py_scripts/cellxgene_processing/qc_pipeline.R"

# Cancella un eventuale file .rds corrotto generato in precedenza
rm -f "/projects/shared/cellxgene_split_by_dataset/blood/qc_results/01209dce-3575-4bed-b1df-129f57fbc031_qc_metrics.rds"

echo "=== Test Singolo File: ${H5AD_FILE} ==="

Rscript "${R_SCRIPT}" "${H5AD_FILE}"

echo "=== Test completato con successo ==="
