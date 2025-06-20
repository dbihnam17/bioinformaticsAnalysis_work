#!/bin/bash
#SBATCH --partition=cpu-short
#SBATCH --cpus-per-task=16
#SBATCH --mem=100GB
#SBATCH --time=12:00:00
#SBATCH --output=cellranger_ref_build.log
#SBATCH --mail-user=your.email@email.com
#SBATCH --mail-type=ALL

### This SLURM script builds the transcriptome required to run cellRanger

# Change to the appropriate directory
cd /path/to/hg38

# Run the Cell Ranger reference build command
cellranger mkref --genome=hg38 \
    --fasta=./GRCh38.p13.genome.fa \
    --genes=./gencode.v38.annotation.gtf \
    --output-dir=./cellranger_ref
