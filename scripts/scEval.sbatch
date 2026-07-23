#!/bin/bash
#SBATCH --job-name=scEvals
#SBATCH --partition=gpu
#SBATCH --gres=gpu:nvidia_h200_nvl_1g.35gb:1
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=70G
#SBATCH --output=/home/%u/scFound_bench/outputs/scEvals-%j.out
#SBATCH --error=/home/%u/scFound_bench/logs/scEvals-%j.err

# specify the paths of the out and err to have them created in a sensible position

set -euo pipefail

echo "Job ID: $SLURM_JOB_ID"
echo "Host: $(hostname)"
echo "Start time: $(date)"
echo "Starting Jupyter with scEval defined environment"

MODEL=${MODEL:-scGPT}

# Load the required environment
module purge
apptainer exec --nv --overlay overlay.img --bind /projects --bind /software scEval.sif \
    conda run -n $MODEL jupyter lab --no-browser --ip 0.0.0.0 --port 8899 --ServerApp.token='stokenazzo'
