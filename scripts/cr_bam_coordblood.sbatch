#!/bin/bash
#SBATCH --job-name=cr_count_citeblood
#SBATCH --partition=cpu
#SBATCH --array=0-10
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02-00:00:00
#SBATCH --output=/projects/shared/intronic_bam/cite_seq_coord_blood/outputs/cellranger-%A_%a.out
#SBATCH --error=/projects/shared/intronic_bam/cite_seq_coord_blood/logs/cellranger-%A_%a.err

set -euo pipefail

FASTQ_BASE="/projects/shared/intronic_bam/cite_seq_coord_blood"
OUT_BASE="/projects/shared/intronic_bam/cite_seq_coord_blood/cellranger_out"

# Reference transcriptome
TRANSCRIPTOME="/projects/shared/intronic_bam/ref/ref/refdata-cellranger-hg19-3.0.0"

# Creazione delle cartelle di base
mkdir -p "${FASTQ_BASE}/logs" "${OUT_BASE}"

module purge
module load cellranger/10.0

# Lettura delle sole cartelle dei campioni (escludendo logs e cellranger_out)
mapfile -t SRRS < <(ls -1 "${FASTQ_BASE}" | grep '^SRR' | sort)

SRR="${SRRS[$SLURM_ARRAY_TASK_ID]}"
FASTQ_DIR="${FASTQ_BASE}/${SRR}"

# Creiamo una cartella temporanea per i link simbolici compatibili con Cell Ranger
TMP_FASTQ_DIR="${OUT_BASE}/${SRR}_tmp_fastqs"
mkdir -p "${TMP_FASTQ_DIR}"

# Creazione dei link con il formato standard accettato da Cell Ranger (_R1_001 e _R2_001)
ln -sf "${FASTQ_DIR}/${SRR}_1.fastq.gz" "${TMP_FASTQ_DIR}/${SRR}_S1_L001_R1_001.fastq.gz"
ln -sf "${FASTQ_DIR}/${SRR}_2.fastq.gz" "${TMP_FASTQ_DIR}/${SRR}_S1_L001_R2_001.fastq.gz"

echo "JobID: ${SLURM_JOB_ID}"
echo "Task : ${SLURM_ARRAY_TASK_ID}"
echo "SRR  : ${SRR}"
echo "FASTQ: ${TMP_FASTQ_DIR}"
echo "OUT  : ${OUT_BASE}"
echo "REF  : ${TRANSCRIPTOME}"
echo "Host : $(hostname)"
echo "Start: $(date)"

# Ci spostiamo nella directory di output generale
cd "${OUT_BASE}"

# Esecuzione di Cell Ranger puntando alla cartella dei link simbolici
cellranger count \
  --id="${SRR}" \
  --create-bam=true \
  --chemistry=SC3Pv2 \
  --transcriptome="${TRANSCRIPTOME}" \
  --fastqs="${TMP_FASTQ_DIR}" \
  --sample="${SRR}" \
  --localcores="${SLURM_CPUS_PER_TASK}" \
  --localmem=64 \
  --include-introns=true \
  --disable-ui

# Opzionale: Rimuove la cartella temporanea dei link simbolici a fine job per tenere pulito
rm -rf "${TMP_FASTQ_DIR}"

echo "End: $(date)"
