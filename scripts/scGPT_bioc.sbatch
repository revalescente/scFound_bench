#!/bin/bash
#SBATCH --job-name=scGPT_bioc
#SBATCH --partition=gpu
#SBATCH --gres=gpu:nvidia_h200_nvl_1g.35gb:1
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=70G
#SBATCH --output=/home/%u/scFound_bench/outputs/scGPT-%j.out
#SBATCH --error=/home/%u/scFound_bench/logs/scGPT-%j.err

# specify the paths of the out and err to have them created in a sensible position

set -euo pipefail

echo "Job ID: $SLURM_JOB_ID"
echo "Host: $(hostname)"
echo "Start time: $(date)"
echo "Starting Jupyter with scGPT environment"

# Load the required environment
module purge
module load micromamba/scGPT

# Start Jupyter on the compute node
jupyter lab --no-browser --ip 0.0.0.0 --port 8877
