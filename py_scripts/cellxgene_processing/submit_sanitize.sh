#!/bin/bash
#SBATCH --job-name=sanitize_cxg
#SBATCH --time=02:00:00
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --output=/home/%u/scFound_bench/outputs/sanitize-%j.out
#SBATCH --error=/home/%u/scFound_bench/logs/sanitize-%j.err

set -euo pipefail

module purge
module load micromamba
micromamba activate nw-scGPTv2

# Percorso principale dove risiedono i dataset divisi
TARGET_DIR="/projects/shared/cellxgene_split_by_dataset"

echo "Avvio della bonifica globale degli indici su: ${TARGET_DIR}"

python /home/vreffo/scFound_bench/py_scripts/cellxgene_processing/sanitize_h5ad.py \
    --base-dir "${TARGET_DIR}" \
    --n-workers ${SLURM_CPUS_PER_TASK}

echo "=== Bonifica Completata ==="
