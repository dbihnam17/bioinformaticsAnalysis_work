#!/bin/bash
#SBATCH --partition=cpu-short
#SBATCH --cpus-per-task=12
#SBATCH --mail-user=your.email@email.com
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=6G
#SBATCH --job-name=mutSigsFinDecomposeTest
#SBATCH --time=2-00:00:00
#SBATCH --chdir=./logs

### This SLURM script submits a job to run sigProfilerDecomposer on a SBS96 matrix

# Activate Conda environment
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate sigprofiler_env

# Run SigProfilerExtractor
python -c "
import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze
Analyze.cosmic_fit(
samples='/path/to/SBS96/Samples.txt', 
output='/path/to/SBS96/decomp_example',
input_type='matrix', 
genome_build='GRCh38', 
cosmic_version=3.4),
exclude_signature_subgroups=['Artifact_signatures']
"
