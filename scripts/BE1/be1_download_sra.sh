#!/bin/bash
#SBATCH --job-name=download_PRJNA1019356
#SBATCH --partition=cpu
#SBATCH --array=0-15
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=home/%u/scFound_bench/outputs/download-%A.out
#SBATCH --error=home/%u/scFound_bench/logs/download-%A.err

set -euo pipefail

# --- PATHS ---
BASE_DIR="/projects/shared/intronic_bam/BE1"
LISTA_SRA="${BASE_DIR}/sra_list.txt"


# --- CARICAMENTO MODULI ---
module purge
module load sra-tools/3.0.3

# --- SELEZIONE DEL CAMPIONE TRAMITE TASK ID ---
if [[ ! -f "$LISTA_SRA" ]]; then
    echo "❌ ERRORE: Il file $LISTA_SRA non esiste!"
    exit 1
fi

mapfile -t SRRS < "$LISTA_SRA"
SRR_ID="${SRRS[$SLURM_ARRAY_TASK_ID]}"

# Definiamo la cartella di destinazione per questo specifico campione
TARGET_DIR="${BASE_DIR}/${SRR_ID}"
mkdir -p "$TARGET_DIR"

echo "====================================================="
echo "JobID    : ${SLURM_JOB_ID}"
echo "Task ID  : ${SLURM_ARRAY_TASK_ID}"
echo "Download : ${SRR_ID}"
echo "Dest     : ${TARGET_DIR}"
echo "Start    : $(date)"
echo "====================================================="

# --- CONFIGURAZIONE SRA TOOLKIT (Opzionale ma sicura) ---
# Forza prefetch a scaricare i file nella cartella specifica anziché nella cache globale di default (~/ncbi)
export NCBI_SETTINGS=" "

# --- DOWNLOAD CON PREFETCH ---
echo ">>> Avvio del download di ${SRR_ID}.sra..."

# Usiamo l'opzione --output-directory per forzare il download dentro $TARGET_DIR
prefetch --output-directory "$TARGET_DIR" "$SRR_ID"

# Alcune versioni di prefetch creano una struttura sottocartelle del tipo: TARGET_DIR/SRR_ID/SRR_ID.sra
# Spostiamo il file direttamente in TARGET_DIR se necessario per uniformarlo ai tuoi script precedenti
if [ -f "${TARGET_DIR}/${SRR_ID}/${SRR_ID}.sra" ]; then
    mv "${TARGET_DIR}/${SRR_ID}/${SRR_ID}.sra" "${TARGET_DIR}/"
    rmdir "${TARGET_DIR}/${SRR_ID}"
fi

# --- VERIFICA ---
if [ -s "${TARGET_DIR}/${SRR_ID}.sra" ]; then
    echo "✅ Download completato con successo: ${TARGET_DIR}/${SRR_ID}.sra"
else
    echo "❌ ERRORE: Il download di $SRR_ID è fallito o il file è vuoto!"
    exit 1
fi

echo "End: $(date)"
echo "-----------------------------------------------------"
