#!/bin/bash
#SBATCH --job-name=fastqcAnalysis
#SBATCH --output=logs/fastQC/fastqc_%j.out
#SBATCH --error=logs/fastQC/fastqc_%j.err
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=cpu-short
#SBATCH --mail-user=your.email@email.edu
#SBATCH --mail-type=ALL

## This script is meant to submit individual fastq files for fastQC analysis
## in a slurm job

## Intended usage:
## sbatch slurmFastQC.sh /path/to/your/file.fq.gz

## The file extension does not matter as long as the data is in fastq format
## Script output will be in the work directory

# Load in fastQC
module load fastqc

# Check if arguments are properly provided
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <full/file/path.fq.gz>"
  exit 1
fi

# Assign user-specified file path
FILE=$1

# Check if user-specified file exists
if [ ! -f "$FILE" ]; then
  echo "Error: File $FILE does not exist"
  exit 1
fi

# Run fastQC with 4 threads
echo "FastQC on $FILE began at $(date)"
fastqc -t 4 "$FILE"
echo "FastQC on $FILE ended at $(date)"