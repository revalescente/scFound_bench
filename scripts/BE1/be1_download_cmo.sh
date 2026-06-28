#!/bin/bash
#SBATCH --job-name=download_CMO_BE1
#SBATCH --partition=cpu
#SBATCH --array=0-7
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=/home/vreffo/scFound_bench/outputs/download-%A_%a.out
#SBATCH --error=/home/vreffo/scFound_bench/logs/download-%A_%a.err

set -euo pipefail

# --- PATHS ---
BASE_DIR="/projects/shared/intronic_bam/BE1"
LISTA_SRA="${BASE_DIR}/sra_list_cmo.txt"

# --- CARICAMENTO MODULI ---
module purge
module load sra-tools/3.0.3

# --- SELEZIONE DEL CAMPIONE TRAMITE TASK ID ---
if [[ ! -f "$LISTA_SRA" ]]; then
    echo "❌ ERRORE: Il file $LISTA_SRA non esiste!"
    exit 1
fi

# Leggiamo il file e puliamo eventuali caratteri nascosti (\r) di Windows
mapfile -t SRRS < <(tr -d '\r' < "$LISTA_SRA")
SRR_ID="${SRRS[$SLURM_ARRAY_TASK_ID]}"

# Definiamo la cartella di destinazione per questo specifico campione
TARGET_DIR="${BASE_DIR}/${SRR_ID}"
mkdir -p "$TARGET_DIR"

echo "====================================================="
echo "JobID    : ${SLURM_JOB_ID}"
echo "Task ID  : ${SLURM_ARRAY_TASK_ID}"
echo "Download : ${SRR_ID} (CMO Library)"
echo "Dest     : ${TARGET_DIR}"
echo "Start    : $(date)"
echo "====================================================="

# --- CONFIGURAZIONE SRA TOOLKIT ---
# Evita che prefetch usi la cartella cache di default in home (~/.ncbi)
export NCBI_SETTINGS=" "

# --- DOWNLOAD CON PREFETCH ---
echo ">>> Avvio del download di ${SRR_ID}.sra..."

# Usiamo l'opzione --output-directory per forzare il download dentro $TARGET_DIR
prefetch --max-size 50G --output-directory "$TARGET_DIR" "$SRR_ID"

# Controllo e riposizionamento se prefetch crea la fastidiosa sottocartella nidificata
if [ -f "${TARGET_DIR}/${SRR_ID}/${SRR_ID}.sra" ]; then
    mv "${TARGET_DIR}/${SRR_ID}/${SRR_ID}.sra" "${TARGET_DIR}/"
    rmdir "${TARGET_DIR}/${SRR_ID}"
fi

# --- VERIFICA FINALE ---
if [ -s "${TARGET_DIR}/${SRR_ID}.sra" ]; then
    echo "✅ Download completato con successo: ${TARGET_DIR}/${SRR_ID}.sra"
else
    echo "❌ ERRORE: Il download di $SRR_ID è fallito o il file è vuoto!"
    exit 1
fi

echo "End: $(date)"
echo "-----------------------------------------------------"
