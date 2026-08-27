#!/bin/bash
#SBATCH --job-name=cxg_qc_arr
#SBATCH --time=01:00:00
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=/home/%u/scFound_bench/outputs/qc_array_%A_%a.out
#SBATCH --error=/home/%u/scFound_bench/logs/qc_array_%A_%a.err

set -euo pipefail

module purge
module load micromamba/bioc-3.23


DATA_DIR="/projects/shared/cellxgene_split_by_dataset"
R_SCRIPT="/home/vreffo/scFound_bench/py_scripts/cellxgene_processing/qc_pipeline.R"

# Scansione dinamica dei file .h5ad
mapfile -t FILES < <(find "${DATA_DIR}" -name "*.h5ad" | sort)

IDX=$((SLURM_ARRAY_TASK_ID - 1))
H5AD_FILE="${FILES[${IDX}]}"

echo "=== Task ID: ${SLURM_ARRAY_TASK_ID} / ${#FILES[@]} ==="
echo "=== File in elaborazione: ${H5AD_FILE} ==="

# 🚀 Usa direttamente il wrapper gestito dal modulo del cluster
Rscript "${R_SCRIPT}" "${H5AD_FILE}"

echo "=== Task ${SLURM_ARRAY_TASK_ID} completato ==="


# To run it on all files in the folders
# sbatch --array=1-$(find /projects/shared/cellxgene_split_by_dataset -name "*.h5ad" | wc -l)%8 py_scripts/cellxgene_processing/submit_qc_array.sh
