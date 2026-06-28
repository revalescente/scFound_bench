#!/bin/bash
#SBATCH --job-name=extract_BE1
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=05:00:00
#SBATCH --output=/home/vreffo/scFound_bench/outputs/extract-%A_%a.out
#SBATCH --error=/home/vreffo/scFound_bench/logs/extract-%A_%a.err
#SBATCH --array=1-16

set -euo pipefail

BASE_DIR="/projects/shared/intronic_bam/BE1"
LIST_FILE="${BASE_DIR}/sra_list.txt"

SRR_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST_FILE" | tr -d '\r ')
if [ -z "$SRR_ID" ]; then exit 1; fi

FASTQ_DIR="${BASE_DIR}/${SRR_ID}"
R1_CR="${FASTQ_DIR}/${SRR_ID}_S1_L001_R1_001.fastq.gz"
R2_CR="${FASTQ_DIR}/${SRR_ID}_S1_L001_R2_001.fastq.gz"
SRA_FILE="${FASTQ_DIR}/${SRR_ID}.sra"

module purge
module load sra-tools/3.0.3
module load pigz/2.8 2>/dev/null || echo "pigz non trovato"

if [[ ! -s "$R1_CR" || ! -s "$R2_CR" ]]; then
    if [ -f "$SRA_FILE" ]; then
        cd "$FASTQ_DIR"
        fasterq-dump --threads "${SLURM_CPUS_PER_TASK}" --mem 28G --split-files "${SRR_ID}.sra"

        if command -v pigz &> /dev/null; then
            pigz -p "${SLURM_CPUS_PER_TASK}" "${SRR_ID}_1.fastq" &
            pigz -p "${SLURM_CPUS_PER_TASK}" "${SRR_ID}_2.fastq" &
            wait
        else
            gzip "${SRR_ID}_1.fastq" &
            gzip "${SRR_ID}_2.fastq" &
            wait
        fi

        mv "${SRR_ID}_1.fastq.gz" "$R1_CR"
        mv "${SRR_ID}_2.fastq.gz" "$R2_CR"
        rm -f "${SRA_FILE}"
    fi
fi
echo "✅ Estrazione completata per $SRR_ID"
