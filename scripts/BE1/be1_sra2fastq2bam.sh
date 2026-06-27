#!/bin/bash
#SBATCH --job-name=cellranger_BE1
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=192G
#SBATCH --time=02-00:00:00
#SBATCH --output=home/%u/scFound_bench/outputs/be1-%A.out
#SBATCH --error=home/%u/scFound_bench/logs/be1-%A.err
#SBATCH --array=1-8

set -euo pipefail

# --- PATHS ---
BASE_DIR="/projects/shared/intronic_bam/BE1"
OUT_BASE="${BASE_DIR}/cellranger_out"
LIST_FILE="${BASE_DIR}/sra_list.txt"
CELLRANGER_REF="/projects/shared/intronic_bam/ref/ref/refdata-cellranger-hg19-3.0.0"

# Estrazione dell'ID del campione corrente dalla riga corrispondente all'ID dell'array
SRR_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST_FILE" | tr -d '\r ')

if [ -z "$SRR_ID" ]; then
    echo "❌ ERRORE: Nessun SRR_ID trovato alla riga ${SLURM_ARRAY_TASK_ID} di $LIST_FILE"
    exit 1
fi

FASTQ_DIR="${BASE_DIR}/${SRR_ID}"
R1_CR="${FASTQ_DIR}/${SRR_ID}_S1_L001_R1_001.fastq.gz"
R2_CR="${FASTQ_DIR}/${SRR_ID}_S1_L001_R2_001.fastq.gz"
SRA_FILE="${FASTQ_DIR}/${SRR_ID}.sra"

mkdir -p "${OUT_BASE}"
mkdir -p "${BASE_DIR}/logs"

# --- CARICAMENTO MODULI ---
module purge
module load cellranger/10.0  # <--- Usa la versione disponibile sul tuo cluster
module load sra-tools/3.0.3
module load samtools/1.21

if ! command -v pigz &> /dev/null; then
    module load pigz/2.8 2>/dev/null || echo "pigz non trovato, userò gzip"
fi

echo "====================================================="
echo "JobID      : ${SLURM_JOB_ID} (Task: ${SLURM_ARRAY_TASK_ID})"
echo "Sample     : ${SRR_ID}"
echo "CPUs       : ${SLURM_CPUS_PER_TASK}"
echo "Start      : $(date)"
echo "====================================================="

# -----------------------------------------------------------
# FASE 1: ESTRAZIONE DA SINGOLO SRA A FASTQ ACCOPPIATI
# -----------------------------------------------------------
if [[ ! -s "$R1_CR" || ! -s "$R2_CR" ]]; then
    if [ -f "$SRA_FILE" ]; then
        echo ">>> Step 1: Estrazione e splitting dal file singolo ${SRR_ID}.sra..."
        cd "$FASTQ_DIR"

        # --split-files si occupa di dividere l'SRA unico nei due FASTQ nativi
        fasterq-dump --threads "${SLURM_CPUS_PER_TASK}" --mem 60G --split-files --progress "${SRR_ID}.sra"

        echo ">>> Step 2: Compressione parallela dei file FASTQ estratti..."
        if command -v pigz &> /dev/null; then
            pigz -p "${SLURM_CPUS_PER_TASK}" "${SRR_ID}_1.fastq" &
            pigz -p "${SLURM_CPUS_PER_TASK}" "${SRR_ID}_2.fastq" &
            wait
        else
            gzip "${SRR_ID}_1.fastq" &
            gzip "${SRR_ID}_2.fastq" &
            wait
        fi

        echo ">>> Step 3: Formattazione nomi per Cell Ranger..."
        mv "${SRR_ID}_1.fastq.gz" "$R1_CR"
        mv "${SRR_ID}_2.fastq.gz" "$R2_CR"

        echo ">>> Step 4: Validazione..."
        if [ -s "$R1_CR" ] && [ -s "$R2_CR" ]; then
            echo "✅ Validazione superata. I file FASTQ sono pronti."
            echo "🗑️ Rimozione del file originale .sra..."
            rm "${SRR_ID}.sra"
        else
            echo "❌ ERRORE: Estrazione o compressione fallita!"
            exit 1
        fi
        cd "$BASE_DIR"
    else
        echo "❌ ERRORE: FASTQ mancanti e nessun file .sra trovato in $FASTQ_DIR!"
        exit 1
    fi
else
    echo "ℹ️ I file FASTQ estratti esistono già. Salto la Fase 1."
fi

# -----------------------------------------------------------
# FASE 2: ALLINEAMENTO CELL RANGER (METODO PAPER)
# -----------------------------------------------------------
echo ">>> Avvio cellranger count per $SRR_ID (Con letture introniche)..."
cd "$OUT_BASE"

cellranger count --id="${SRR_ID}" \
                 --transcriptome="${CELLRANGER_REF}" \
                 --fastqs="${FASTQ_DIR}" \
                 --sample="${SRR_ID}" \
                 --localcores="${SLURM_CPUS_PER_TASK}" \
                 --localmem=60 \
                 --include-introns=true \
                 --chemistry=auto

echo "====================================================="
echo "✅ Job completato per $SRR_ID"
echo "ℹ️ FASTQ mantenuti in: $FASTQ_DIR"
echo "End        : $(date)"
echo "====================================================="
