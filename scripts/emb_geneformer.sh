#!/bin/bash
#SBATCH --job-name=geneformer_emb
#SBATCH --partition=gpu
#SBATCH --gres=gpu:nvidia_h200_nvl:1          # Richiede 1 GPU H200 intera (se preferisci un profilo MIG più grande usa: nvidia_h200_nvl_2g.70gb:1)
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --output=/home/%u/scFound_bench/outputs/geneformer-%j.out
#SBATCH --error=/home/%u/scFound_bench/logs/geneformer-%j.err

set -euo pipefail

echo "Job ID: $SLURM_JOB_ID"
echo "Host: $(hostname)"
echo "Start time: $(date)"

# Carica l'ambiente
module purge
module load micromamba/huggingface

# Esegui lo script Python
echo "Starting Geneformer embedding extraction..."
python /home/vreffo/scFound_bench/py_scripts/Geneformer/emb_extraction.py

echo "Job finished at: $(date)"
