#!/bin/bash
#SBATCH --job-name=geneformer_emb
#SBATCH --partition=gpu
#SBATCH --gres=gpu:nvidia_h200_nvl:1          # Richiede 1 GPU H200 intera (se preferisci un profilo MIG più grande usa: nvidia_h200_nvl_2g.70gb:1)
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --output=/home/%u/scFound_bench/outputs/geneformer-%j.out
#SBATCH --error=/home/%u/scFound_bench/logs/geneformer-%j.err

set -euo pipefail

echo "Job ID: $SLURM_JOB_ID"
echo "Host: $(hostname)"
echo "Start time: $(date)"

# Carica l'ambiente
module purge
module load micromamba/huggingface

# GPU monitoring
GPU_LOG="/home/vreffo/scFound_bench/logs/gpu_usage_geneformer_${SLURM_JOB_ID}.csv"

# every 30 seconds in background
nvidia-smi --query-gpu=timestamp,name,utilization.gpu,utilization.memory,memory.used,memory.total \
           --format=csv -l 30 > "$GPU_LOG" &
GPU_PID=$!

# Tokenization of the input data
echo "Tokenizing input data..."
python /home/vreffo/scFound_bench/py_scripts/Geneformer/tokenizer.py

# Embedding extraction
echo "Starting Geneformer embedding extraction..."
python /home/vreffo/scFound_bench/py_scripts/Geneformer/emb_extraction.py

# Stop GPU monitoring
kill $GPU_PID
echo "GPU log saved to: $GPU_LOG"

echo "Job finished at: $(date)"
