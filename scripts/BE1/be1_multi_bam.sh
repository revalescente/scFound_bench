#!/bin/bash
#SBATCH --job-name=multi_BE1
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=02-00:00:00
#SBATCH --output=/home/vreffo/scFound_bench/outputs/multi_be1-%J.out
#SBATCH --error=/home/vreffo/scFound_bench/logs/multi_be1-%J.err

set -euo pipefail

# --- PATHS ---
BASE_DIR="/projects/shared/intronic_bam/BE1"
OUT_BASE="${BASE_DIR}/cellranger_multi_out"
CONFIG_CSV="${BASE_DIR}/multi_config.csv"

# --- VERIFICA SE IL FILE ESISTE GIÀ ---
if [[ ! -f "$CONFIG_CSV" ]]; then
    echo "❌ ERRORE: Il file di configurazione $CONFIG_CSV non esiste nella directory!"
    exit 1
fi

mkdir -p "${OUT_BASE}"
cd "${OUT_BASE}"

# --- CARICAMENTO MODULI ---
module purge
module load cellranger/10.0

echo ">>> Avvio cellranger multi usando il file esistente: ${CONFIG_CSV}..."

# Lancio di cellranger multi puntando direttamente al tuo file statico
cellranger multi --id=BE1_demux \
                 --csv="${CONFIG_CSV}" \
                 --localcores="${SLURM_CPUS_PER_TASK}" \
                 --localmem=120

echo "====================================================="
echo "✅ Pipeline completata con successo!"
echo "====================================================="
