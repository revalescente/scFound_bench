#!/bin/bash
#SBATCH --job-name=merge_bams_BE1
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=home/%u/scFound_bench/outputs/be1merge-%A.out
#SBATCH --error=home/%u/scFound_bench/logs/be1merge-%A.err

set -euo pipefail

BASE_DIR="/projects/shared/intronic_bam/BE1"
OUT_BASE="${BASE_DIR}/cellranger_out"


module purge
module load samtools/1.21

echo "====================================================="
echo "JobID  : ${SLURM_JOB_ID}"
echo "Start  : $(date)"
echo "====================================================="

echo ">>> Recupero dei file BAM generati da Cell Ranger..."
BAM_LIST=()

for srr_dir in "${OUT_BASE}"/SRR*/outs; do
    if [ -f "${srr_dir}/possorted_genome_bam.bam" ]; then
        echo "Trovato: ${srr_dir}/possorted_genome_bam.bam"
        BAM_LIST+=("${srr_dir}/possorted_genome_bam.bam")
    fi
done

if [ ${#BAM_LIST[@]} -eq 0 ]; then
    echo "❌ ERRORE: Nessun file BAM trovato in ${OUT_BASE}!"
    exit 1
fi

OUTPUT_BAM="${OUT_BASE}/BE1_all_runs_merged.bam"

echo ">>> Unione di ${#BAM_LIST[@]} file BAM in corso..."
samtools merge -@ "${SLURM_CPUS_PER_TASK}" "${OUTPUT_BAM}" "${BAM_LIST[@]}"

echo ">>> Generazione dell'indice (.bai)..."
samtools index -@ "${SLURM_CPUS_PER_TASK}" "${OUTPUT_BAM}"

echo "====================================================="
echo "✅ Aggregazione completata!"
echo "ℹ️ File finale: ${OUTPUT_BAM}"
echo "End  : $(date)"
echo "====================================================="
