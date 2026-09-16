#!/bin/bash
#SBATCH --job-name=scgpt_embeds
#SBATCH --time=02:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=/home/%u/scFound_bench/outputs/scgpt_embeds-%j.out
#SBATCH --error=/home/%u/scFound_bench/logs/scgpt_embeds-%j.err

set -euo pipefail

module purge
module load micromamba
micromamba activate nw-scGPTv2

echo "=== scGPT embedding extraction ==="
echo "Job ID:  ${SLURM_JOB_ID}"
echo "Node:    ${SLURMD_NODENAME}"
echo "GPU:     ${CUDA_VISIBLE_DEVICES:-none}"
echo "Start:   $(date)"

python /home/vreffo/scFound_bench/py_scripts/scGPT/embeds_extraction/all_embeds.py

echo "=== Completato: $(date) ==="
