#!/bin/bash
#SBATCH --partition=cpu-short
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=your.email@email.com
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=2G
#SBATCH -J proj_bam_idxstats
#SBATCH --export=USER
#SBATCH --time 4-00:00:00
#SBATCH --chdir=./logs

### This script is designed to automate the creation of idxstats reports for
### BAM files in a user-specified directory

# Exit on error
set -e

# Print commands
set -x

# Load in samtools module
module load samtools/1.18

BAM_DIR="$1" # Input bam directory
WORK_DIR="$2" # Output/working directory
OUT_DIR="${WORK_DIR}/bamStats"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Loop through each BAM file in the user-specified input directory
for bam_file in "${BAM_DIR}"/*.bam; do
    # Check if the file contains '~' and skip it if true
    # This was an issue in a particular Novogene WGS batch
    if [[ "$bam_file" == *"~"* ]]; then
        echo "Skipping file with extra string: ${bam_file}"
        continue
    fi

    # Generate output filename in the correct directory
    OUTNAME="${OUT_DIR}/$(basename "${bam_file}" .bam).idxstats"

    # Make idxstats report
    samtools idxstats "${bam_file}" > "${OUTNAME}"
done
