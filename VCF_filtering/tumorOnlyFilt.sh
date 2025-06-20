#!/bin/bash
#SBATCH --job-name=filterTumorOnlyVcfs
#SBATCH --output=logs/filterTumorOnlyVcfs_%A_%a.out
#SBATCH --error=logs/filterTumorOnlyVcfs_%A_%a.err
#SBATCH --partition=cpu-short
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-user=your.email@email.com
#SBATCH --mail-type=ALL
#SBATCH --array=0-74

## This SLURM script was written to filter a batch of tumor-only VCFs based off
## various metrics in an attempt to run mutational signatures in parallel with
## paired WGS samples. The filtering worked, the signatures were not so comparable

## Possibly play around with the filtering conditions to find a better medium

# Load bcftools
module load bcftools

# Set input/output directories
INPUT_DIR="/path/to/tumor_only.vcf"
OUTPUT_DIR="/path/to/output"

# Create the output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Filtering condition
FILTER_CONDITION='FILTER="PASS" && FORMAT/DP>=5 && CADD_v13_GRCh38.PHRED>=5'

# Get list of uncompressed VCF files in the input directory
vcf_files=(${INPUT_DIR}/*.vcf)

# Get the VCF file corresponding to the SLURM_ARRAY_TASK_ID
vcf="${vcf_files[${SLURM_ARRAY_TASK_ID}]}"

# Extract basename (e.g., sample1.vcf -> sample1)
base=$(basename "${vcf}" .vcf)

# Set the output file name and path
output_file="${OUTPUT_DIR}/${base}.filtered.vcf"

# Apply bcftools filter using the defined condition
bcftools filter -i "${FILTER_CONDITION}" "${vcf}" -o "${output_file}"
