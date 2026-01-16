#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.errset -euo pipefail

mkdir -p logs

python run.py
