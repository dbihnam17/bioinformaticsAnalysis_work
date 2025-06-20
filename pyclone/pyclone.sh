#!/bin/bash
#SBATCH --job-name=sample_pyclone
#SBATCH --partition=cpu-short
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=your.email@email.com
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=4G
#SBATCH --time=4-00:00:00
#SBATCH --chdir=./logs



# Source conda
source /path/to/miniconda3/etc/profile.d/conda.sh

# Activate pyclone-vi environment
conda activate pyclone-vi

# Run pyclone
pyclone-vi fit \
  -i /path/to/pyclone/sample_pycloneInput.tsv \
  -o /path/to/pyclone/results/sample_pycloneResults.h5 \
  -c 10 \
  --num-restarts 3 \
  --num-annealing-steps 5 \
  --max-iters 100000 \
  --print-freq 500 \
  --num-threads 16 \
  --seed 42
