#!/bin/bash
#SBATCH --job-name=Sample_ID_cellranger
#SBATCH --partition=cpu-short
#SBATCH --output=logs/cellRangerCount/cellrangerSample_ID_%j.out
#SBATCH --error=logs/cellRangerCount/cellrangerSample_ID_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=80GB
#SBATCH --time=24:00:00
#SBATCH --mail-user=your.email@email.com
#SBATCH --mail-type=ALL

### This SLURM script allows you to submit a job to run cellRanger on scRNA
### fastq files

# Increase process limit
ulimit -u 4125762

# Run Cell Ranger count
cellranger count --id=sample_ID \
    --fastqs=/path/to/sample/dir \
    --sample=sampleName \
    --transcriptome=/path/to/hg38/cellranger_ref \
    --create-bam=true \
    --expect-cells=5000 # Optional, allows you to set expected number of cells
