#!/bin/bash
#SBATCH --job-name=inspect_all_h5ad
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/home/%u/scFound_bench/outputs/inspect_all-%j.out
#SBATCH --error=/home/%u/scFound_bench/logs/inspect_all-%j.err

module purge
module load micromamba
micromamba activate nw-scGPTv2

python inspect_cxg_all.py
