#!/bin/bash
#SBATCH --job-name=starsolo_citeblood
#SBATCH --partition=cpu
#SBATCH --array=0-10
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02-00:00:00
#SBATCH --output=/projects/shared/intronic_bam/cite_seq_coord_blood/logs/starsolo-%A_%a.out
#SBATCH --error=/projects/shared/intronic_bam/cite_seq_coord_blood/logs/starsolo-%A_%a.err

set -euo pipefail

# --- PATHS ---
BASE_DIR="/projects/shared/intronic_bam/cite_seq_coord_blood"
OUT_BASE="${BASE_DIR}/starsolo_out"
LOG_DIR="${BASE_DIR}/logs"

# Verify these match your actual shared reference paths on the cluster
STAR_INDEX="/projects/shared/intronic_bam/ref/ref/STAR_hg19_10x_optimal"
WHITELIST="/projects/shared/intronic_bam/ref/ref/737K-august-2016.txt"

# Create necessary directories
mkdir -p "${LOG_DIR}" "${OUT_BASE}"

# --- MODULES ---
module purge
module load star/2.7.11b
module load samtools/1.21

# --- SAMPLE SELECTION ---
# Read only directories starting with SRR into an array
mapfile -t SRRS < <(ls -1 "${BASE_DIR}" | grep '^SRR' | sort)
SRR_ID="${SRRS[$SLURM_ARRAY_TASK_ID]}"
FASTQ_DIR="${BASE_DIR}/${SRR_ID}"

R1="${FASTQ_DIR}/${SRR_ID}_1.fastq.gz"
R2="${FASTQ_DIR}/${SRR_ID}_2.fastq.gz"

echo "====================================================="
echo "JobID: ${SLURM_JOB_ID}"
echo "Task : ${SLURM_ARRAY_TASK_ID}"
echo "SRR  : ${SRR_ID}"
echo "Host : $(hostname)"
echo "Start: $(date)"
echo "====================================================="

if [[ ! -s "$R1" || ! -s "$R2" ]]; then
    echo "❌ [SKIP] $SRR_ID: Missing or empty fastq.gz files."
    exit 1
fi

echo ">>> Profiling $SRR_ID to determine library type (RNA vs ADT)..."

# Temporary isolated directory for the test alignment
TEST_DIR="${OUT_BASE}/${SRR_ID}_test_tmp"
mkdir -p "${TEST_DIR}"

# Temporarily disable strict error checking to prevent 'broken pipe' crashes during head
set +e +o pipefail
zcat "$R1" 2>/dev/null | head -n 400000 | gzip > "${TEST_DIR}/_test_R1.fastq.gz"
zcat "$R2" 2>/dev/null | head -n 400000 | gzip > "${TEST_DIR}/_test_R2.fastq.gz"
set -e -o pipefail

# Run a silent, fast test alignment (MAPPING R2 ONLY)
STAR --runThreadN "${SLURM_CPUS_PER_TASK}" --genomeDir "$STAR_INDEX" \
     --readFilesIn "${TEST_DIR}/_test_R2.fastq.gz" \
     --readFilesCommand zcat --outSAMtype None \
     --outFileNamePrefix "${TEST_DIR}/_test_" >/dev/null 2>&1 || true

# Safely check if the log file was created and extract map rate
if [ -f "${TEST_DIR}/_test_Log.final.out" ]; then
    MAP_RATE=$(grep -m1 "Uniquely mapped reads %" "${TEST_DIR}/_test_Log.final.out" | awk -F'|' '{gsub(/ /,"",$2); print $2}' | tr -d '%')
else
    echo "⚠️  WARNING: STAR failed to run properly on $SRR_ID. Check your STAR_INDEX path."
    MAP_RATE=0
fi

MAP_RATE=${MAP_RATE:-0}
MAP_INT=${MAP_RATE%.*}

# Cleanup test files immediately
rm -rf "${TEST_DIR}"

if [[ "$MAP_INT" -lt 20 ]]; then
    echo "[SKIP] $SRR_ID: Uniquely mapped rate is only ${MAP_RATE}%. This is likely an ADT/Protein library."
    echo "End: $(date)"
    exit 0
fi

echo "[RUN] $SRR_ID: Uniquely mapped rate is ${MAP_RATE}%. This is an RNA library. Starting full STARsolo alignment..."

OUT_DIR="${OUT_BASE}/${SRR_ID}"
mkdir -p "$OUT_DIR"

# Production STARsolo command
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

echo ">>> Indexing BAM file for $SRR_ID..."
samtools index "$OUT_DIR/Aligned.sortedByCoord.out.bam"

echo "✅ $SRR_ID complete. BAM ready for DropletQC."
echo "End: $(date)"
echo "-----------------------------------------------------"
