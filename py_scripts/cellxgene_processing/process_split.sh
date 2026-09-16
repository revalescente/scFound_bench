#!/bin/bash
#SBATCH --job-name=cxg_split
#SBATCH --time=04:00:00
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --array=1-9
#SBATCH --output=/home/%u/scFound_bench/outputs/split-%A_%a.out
#SBATCH --error=/home/%u/scFound_bench/logs/split-%A_%a.err

set -euo pipefail

module purge
module load micromamba
micromamba activate nw-scGPTv2

BASE_DIR="/projects/shared/cellxgene"
OUTPUT_DIR="/projects/shared/cellxgene_split_by_dataset"

# Lista delle cartelle da elaborare (escludendo faiss_index e cartelle di output)
TISSUES=($(ls -d ${BASE_DIR}/*/ | grep -v 'faiss_index' | xargs -n1 basename))

# Seleziona il tessuto per questo array job
TISSUE_NAME="${TISSUES[$((SLURM_ARRAY_TASK_ID-1))]}"
TISSUE_PATH="${BASE_DIR}/${TISSUE_NAME}"

echo "Avvio split per il tessuto: ${TISSUE_NAME}"

python /home/vreffo/scFound_bench/py_scripts/cellxgene_processing/split_datasets.py \
    --tissue-dir "${TISSUE_PATH}" \
    --output-dir "${OUTPUT_DIR}" \
    --n-workers ${SLURM_CPUS_PER_TASK}
