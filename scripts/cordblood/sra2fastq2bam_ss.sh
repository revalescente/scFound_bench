#!/bin/bash
#SBATCH --job-name=starsolo_SRR5808751
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02-00:00:00
#SBATCH --output=/home/vreffo/scFound_bench/outputs/pipe-%A.out
#SBATCH --error=/home/vreffo/scFound_bench/logs/pipe-%A.err

set -euo pipefail

# --- PATHS ---
BASE_DIR="/projects/shared/intronic_bam/cite_seq_coord_blood"
OUT_BASE="${BASE_DIR}/starsolo_out"

# Definiamo staticamente il tuo campione
SRR_ID="SRR5808751"
FASTQ_DIR="${BASE_DIR}/${SRR_ID}"

# Riferimenti STAR e Whitelist
STAR_INDEX="/projects/shared/intronic_bam/ref/ref/STAR_hg19_10x_optimal"
WHITELIST="/projects/shared/intronic_bam/ref/ref/737K-august-2016.txt"

# Directory per i file FASTQ attesi
R1="${FASTQ_DIR}/${SRR_ID}_1.fastq.gz"
R2="${FASTQ_DIR}/${SRR_ID}_2.fastq.gz"
SRA_FILE="${FASTQ_DIR}/${SRR_ID}.sra"

# Creazione directory di log e output
mkdir -p "${OUT_BASE}"

# --- CARICAMENTO MODULI ---
module purge
module load star/2.7.11b
module load samtools/1.21
module load sra-tools/3.0.3

# Verifica se pigz è disponibile come modulo, altrimenti usa gzip standard
if ! command -v pigz &> /dev/null; then
    module load pigz/2.8 2>/dev/null || echo "pigz non trovato, userò gzip (più lento)"
fi

echo "====================================================="
echo "JobID  : ${SLURM_JOB_ID}"
echo "Sample : ${SRR_ID}"
echo "CPUs   : ${SLURM_CPUS_PER_TASK}"
echo "Mem    : ${SLURM_MEM_PER_NODE}MB"
echo "Start  : $(date)"
echo "====================================================="

# -----------------------------------------------------------
# FASE 1: ESTRAZIONE DA SRA A FASTQ (CON SUPER-PARALLELISMO)
# -----------------------------------------------------------
if [[ ! -s "$R1" || ! -s "$R2" ]]; then
    if [ -f "$SRA_FILE" ]; then
        echo ">>> Step 1: Estrazione FASTQ ultra-rapida da ${SRR_ID}.sra..."
        cd "$FASTQ_DIR"

        # Uso massimo dei thread e 60GB di cache RAM per fasterq-dump
        fasterq-dump --threads "${SLURM_CPUS_PER_TASK}" --mem 60G --split-files --progress "${SRR_ID}.sra"

        echo ">>> Step 2: Compressione parallela dei file FASTQ..."
        if command -v pigz &> /dev/null; then
            # pigz comprime usando tutte le CPU a disposizione
            pigz -p "${SLURM_CPUS_PER_TASK}" "${SRR_ID}_1.fastq" &
            pigz -p "${SLURM_CPUS_PER_TASK}" "${SRR_ID}_2.fastq" &
            wait # Aspetta che entrambe le compressioni parallele finiscano
        else
            gzip "${SRR_ID}_1.fastq" &
            gzip "${SRR_ID}_2.fastq" &
            wait
        fi

        echo ">>> Step 3: Validazione dei file generati..."
        if [ -s "$R1" ] && [ -s "$R2" ]; then
            echo "✅ Validazione superata."
            echo "🗑️ Eliminazione del file .sra originale..."
            rm "${SRR_ID}.sra"
        else
            echo "❌ ERRORE: Validazione fallita!"
            exit 1
        fi
        cd "$BASE_DIR"
    else
        echo "❌ ERRORE: File FASTQ mancanti E .sra non trovato!"
        exit 1
    fi
else
    echo "ℹ️ I file FASTQ compressi esistono già. Salto l'estrazione."
fi

# -----------------------------------------------------------
# FASE 2: PROFILING DELLA LIBRERIA
# -----------------------------------------------------------
echo ">>> Profiling rapido della libreria..."
TEST_DIR="${OUT_BASE}/${SRR_ID}_test_tmp"
mkdir -p "${TEST_DIR}"

set +e +o pipefail
zcat "$R1" 2>/dev/null | head -n 400000 | gzip > "${TEST_DIR}/_test_R1.fastq.gz"
zcat "$R2" 2>/dev/null | head -n 400000 | gzip > "${TEST_DIR}/_test_R2.fastq.gz"
set -e -o pipefail

STAR --runThreadN "${SLURM_CPUS_PER_TASK}" --genomeDir "$STAR_INDEX" \
     --readFilesIn "${TEST_DIR}/_test_R2.fastq.gz" \
     --readFilesCommand zcat --outSAMtype None \
     --outFileNamePrefix "${TEST_DIR}/_test_" >/dev/null 2>&1 || true

if [ -f "${TEST_DIR}/_test_Log.final.out" ]; then
    MAP_RATE=$(grep -m1 "Uniquely mapped reads %" "${TEST_DIR}/_test_Log.final.out" | awk -F'|' '{gsub(/ /,"",$2); print $2}' | tr -d '%')
else
    MAP_RATE=0
fi

MAP_RATE=${MAP_RATE:-0}
MAP_INT=${MAP_RATE%.*}
rm -rf "${TEST_DIR}"

if [[ "$MAP_INT" -lt 20 ]]; then
    echo "[SKIP] $SRR_ID: Mapping rate unico al ${MAP_RATE}%. Probabile ADT."
    exit 0
fi

# -----------------------------------------------------------
# FASE 3: ALLINEAMENTO COMPLETO (STARsolo) -> BAM
# -----------------------------------------------------------
echo "[RUN] $SRR_ID: Mapping rate al ${MAP_RATE}%. Avvio STARsolo a 40 core..."

OUT_DIR="${OUT_BASE}/${SRR_ID}"
mkdir -p "$OUT_DIR"

STAR --runThreadN "${SLURM_CPUS_PER_TASK}" \
      --genomeDir "$STAR_INDEX" \
      --readFilesIn "$R2" "$R1" \
      --readFilesCommand zcat \
      --soloType CB_UMI_Simple \
      --soloCBwhitelist "$WHITELIST" \
      --soloCBlen 16 \
      --soloUMIlen 9 \
      --soloBarcodeReadLength 0 \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
      --outFileNamePrefix "$OUT_DIR/"

echo ">>> Indicizzazione del file BAM con samtools (multi-thread)..."
# Anche samtools index ora usa più thread (-@)
samtools index -@ "${SLURM_CPUS_PER_TASK}" "$OUT_DIR/Aligned.sortedByCoord.out.bam"

echo "✅ Pipeline completata in modalità ULTRA per $SRR_ID."
echo "End: $(date)"
echo "-----------------------------------------------------"
