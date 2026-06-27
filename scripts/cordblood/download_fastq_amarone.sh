#!/bin/bash
#SBATCH --job-name=dl_fastqs
#SBATCH --partition=cpu
#SBATCH --time=5-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --output=scFound_bench/outputs/slurm-%j.out
#SBATCH --error=scFound_bench/logs/slurm-%j.err

set -euo pipefail

SRC_HOST="revalescente@amarone.stat.unipd.it"
SRC_BASE="~/cite_seq_coord_blood_raw_data"
DEST_BASE="/projects/shared/intronic_bam/cite_seq_coord_blood"

mkdir -p "$DEST_BASE/logs"

echo "Job ID: $SLURM_JOB_ID"
echo "Host: $(hostname)"
echo "Start time: $(date)"
echo "Destination: $DEST_BASE"
echo "Source: ${SRC_HOST}:${SRC_BASE}"
echo

# List SRR directories on the remote, then copy only *_1.fastq.gz and *_2.fastq.gz
ssh "$SRC_HOST" "ls -1d ${SRC_BASE}/SRR*/ 2>/dev/null" | while read -r remote_dir; do
  srr_id="$(basename "${remote_dir%/}")"
  echo "==> Copying FASTQs for ${srr_id}"

  mkdir -p "${DEST_BASE}/${srr_id}"

  # Copy only the fastq.gz files we want.
  # If either file is missing, warn and continue.
  if scp -p \
      "${SRC_HOST}:${SRC_BASE}/${srr_id}/${srr_id}_1.fastq.gz" \
      "${SRC_HOST}:${SRC_BASE}/${srr_id}/${srr_id}_2.fastq.gz" \
      "${DEST_BASE}/${srr_id}/"
  then
    echo "  [OK] ${srr_id}"
  else
    echo "  [WARN] Missing FASTQs or scp failed for ${srr_id}, skipping."
  fi
done

echo
echo "End time: $(date)"
echo "Done. FASTQs copied into: $DEST_BASE"
