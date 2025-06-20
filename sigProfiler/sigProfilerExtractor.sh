#!/bin/bash
#SBATCH --job-name=yourJobName
#SBATCH --partition=cpu-short
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=your.email@email.com
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=2G
#SBATCH --time=4-00:00:00
#SBATCH --chdir=./logs

## This script is intended to run sigProfilerExtractor for WGS mutational signature
## extraction as a slurm job. Change the slurm instructions and file paths as needed.

## sigProfiler is not my own function, and can be found at https://github.com/AlexandrovLab/SigProfilerExtractor

## nmf_replicates=100 will run quickly and should be used for a quick look at the data
## nmf_replicates=1000 will take much more time and is better suited for finalized data

# Source conda
source /path/to/miniconda3/etc/profile.d/conda.sh

# Activate sigprofiler_env
conda activate sigprofiler_env

# SigProfilerExtractor usage
# sig.sigProfilerExtractor(fileType, outputDir, inputDir, reference_genome,
#                          minimum_signatures, maximum_signatures, nmf_replicates,
#                          cpu_count)

# Run SigProfilerExtractor
python -c "
from SigProfilerExtractor import sigpro as sig
sig.sigProfilerExtractor(
    'vcf',
    '/path/to/output', 
    '/path/to/input', 
    reference_genome='GRCh38', 
    minimum_signatures=1, 
    maximum_signatures=12, 
    nmf_replicates=1000,
    cpu=8
)
"
