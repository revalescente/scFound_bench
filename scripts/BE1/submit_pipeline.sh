#!/bin/bash
set -euo pipefail

# Sostituisci con il percorso reale della tua cartella di lavoro se necessario
BASE_DIR="/projects/shared/intronic_bam/BE1"

cd "$BASE_DIR"

echo "=========================================================="
echo "      🚀 AVVIO AUTOMATIZZATO PIPELINE CELL RANGER         "
echo "=========================================================="

# 1. Sottomissione del Job Array per i singoli campioni
# Cattura l'output del comando sbatch (es: "Submitted batch job 123456")
ARRAY_OUTPUT=$(sbatch /home/vreffo/scFound_bench/scripts/BE1/be1_sra2fastq2bam.sh)
echo "$ARRAY_OUTPUT"

# Estrae solo il numero di Job ID dall'output usando grep e cut
ARRAY_JOB_ID=$(echo "$ARRAY_OUTPUT" | grep -oE '[0-9]+')

echo "📌 Job Array registrato con ID: $ARRAY_JOB_ID"
echo "⏳ I singoli campioni (1-8) stanno venendo elaborati in parallelo..."
echo "----------------------------------------------------------"

# 2. Sottomissione dello script di merge con dipendenza afterok
# L'aggregazione partirà SOLO se tutti i task dell'array terminano con successo (exit code 0)
MERGE_OUTPUT=$(sbatch --dependency=afterok:${ARRAY_JOB_ID} /home/vreffo/scFound_bench/scripts/BE1/be1_mergebam.sh)
echo "$MERGE_OUTPUT"

echo "----------------------------------------------------------"
echo "✅ Tutto impostato correttamente!"
echo "Il file BAM aggregato verrà creato automaticamente al termine dell'array."
echo "Puoi monitorare lo stato con il comando: squeue -u \$USER"
echo "=========================================================="
