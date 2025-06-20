#!/bin/sh
#SBATCH --partition=cpu-short
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=your.email@email.com
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=8G
#SBATCH -J bamFiltering
#SBATCH --export=USER
#SBATCH --time 1-00:00:00
#SBATCH --chdir=./logs

### This script is intended for batch trimming of WGS bam files
### Bam files are trimmed around a gene bed to only contain data in those regions
### Decreases WGS bam size, and allows for easy data transfer and interactive QC in IGV

## Originally based on an internal script shared by a colleague
## Modified for batch use and generalization

set -e
set -x

module load samtools/1.18

BAM_DIR=/path/to/bams/
WORK_DIR=/path/to/work/dir/
GENE_BED="${WORK_DIR}/genes.bed" # This is a gene bed of genes to include in the trimmed bam
OUT_DIR="${WORK_DIR}/outDir/"

# Create output directory if it doesn't exist already
mkdir -p "$OUT_DIR"

# Loop through each BAM file in the directory
for bam_file in "${BAM_DIR}"/*.bam; do
    # OPTIONAL: Check if the file contains '~' and skip it if true
    # Some bams containing '~' were erroneously included in a data delivery once
    if [[ "$bam_file" == *"~"* ]]; then
        echo "Skipping file with extra string: ${bam_file}"
        continue
    fi

    OUTNAME=$(basename "${bam_file}" .final.bam)

    #Extract regions, sort, and index the BAM file
    samtools view -h -b "${bam_file}" -o "${OUT_DIR}/${OUTNAME}.final.sub.temp.bam" -L "$GENE_BED"
    samtools sort "${OUT_DIR}/${OUTNAME}.final.sub.temp.bam" -O BAM -o "${OUT_DIR}/${OUTNAME}.final.sub.sorted.bam"
    samtools index "${OUT_DIR}/${OUTNAME}.final.sub.sorted.bam"

    #Remove the temporary file
    rm "${OUT_DIR}/${OUTNAME}.final.sub.temp.bam"
    
    echo "Processing file: ${bam_file}"
done


