#!/bin/sh
#SBATCH --partition=cpu-short
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=your.email@email.com
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=2G
#SBATCH -J sample_ATAC_blacklistSortIndex
#SBATCH --export=USER
#SBATCH --time 4-00:00:00
#SBATCH --chdir=./logs

blacklistBed='/path/to/hg38-blacklist.bed'
ID="$1" # User supplied sample ID

# Load modules
module load samtools
module load bedtools

# Remove reads within the blacklist regions
bedtools intersect -nonamecheck -v -abam "$ID.filtered.bam" -b $blacklistBed > "$ID.tmp.bam"

# Sort and index the bam file
samtools sort -O bam -o "$ID.blacklist-filtered.bam" "$ID.tmp.bam"
samtools index "$ID.blacklist-filtered.bam"

# Remove temporary bam
rm "$ID.tmp.bam"
