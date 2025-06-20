#!/bin/sh
#SBATCH --partition=cpu-long
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=your.email@email.com
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=2G
#SBATCH -J ATACseq
#SBATCH --export=USER
#SBATCH --time 7-00:00:00
#SBATCH --chdir=./logs

### Due to development of ATACseq nextflow pipelines on GCP, production of this pipeline has been abandoned.

### This script is based on https://github.com/CebolaLab/ATAC-seq

### I was working on this piece-by-piece, then was going to connect it all in this script. Those scripts should also be in this repo folder.

# The purpose of this script is to complete ATAC-seq QC, alignment, peak calling, and visualization

# The pipeline is based on the pipeline described by the CebolaLab/ATAC-seq pipeline, but has been
# adapted into a function that can handle a variety of user inputs and be submitted for jobs with the SLURM scheduler

# Sample execution: bash ./ATACseq.sh sampleID /path/to/fastq_1 /path/to/fastq_2 /path/to/output/directory hg38/19

# TEMPORARY DEVELOPMENT NOTES: 
#   GENERAL:
#     - Make sure you are calling user-arguments correctly
#     - Test code chunks interactively before running full script
#     - Ask around to learn about the HPC resources necessary for a job like this
#   TO-DO:
#     - Make sure all required modules are loaded properly
#     - Check if we need to keep or remove multi-mapped reads: KEEP
#     - Check if we need to shift alignments: SHIFT
#     - Check MACS2 peak calling output file specification
#     - When we migrate to GCP, fix file paths for references and bed files

### Define user-specified input arguments as variables -------------------------

ID="$1" # Sample ID# to use for downstream file names
sample1="$2" # fastq_1
sample2="$3" # fastq_2
outDir="$4" # Work/output directory
hgRef="$5" # Human reference genome (hg38/hg19)


### PRE-ALIGNMENT QC -----------------------------------------------------------

# Load fastqc module
module load fastqc

# Generate initial QC reports
fastqc "$sample1" -d "$outDir" -o "$outDir"
fastqc "$sample2" -d "$outDir" -o "$outDir"

# Load fastp module
module load fastp

# Use fastp to remove adapter sequences
fastp -i "$sample1" -I "$sample2" -o "$ID"_1.trimmed.fastq.gz -O "$ID"_2.trimmed.fastq.gz --detect_adapter_for_pe -j "$ID".fastp.json -h "$ID".fastp.html


### ALIGNMENT ------------------------------------------------------------------

# Load bowtie2 module
module load bowtie2

# Load samtools module
module load samtools

# Align sample to user-specified reference genome (hg19/hg38)
if [ "$hgRef" == "hg38" ]; then
  # Set bt2idx variable to hg38 location
  bt2idx=/path/to/hg38
  blacklistBed=./hg38-blacklist.bed
  estGenSize=2862010428 # Assuming 150bp read length
  bowtie2 --local --very-sensitive --no-mixed --no-discordant -phred33 -I 10 -X 700 -p 32 -x $bt2idx -1 "$ID"_1.trimmed.fastq.gz -2 "$ID"_2.trimmed.fastq.gz | samtools view -bS - > "$ID".bam
  # bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 32 -x $bt2idx -1 sample_1.trimmed.fastq.gz -2 sample_2.trimmed.fastq.gz | samtools view -bS - > sample.bam
  # This is a code adapted from what a collaborator sent. I added in --phred33, changed -I to 10 (from 25) and added -p 32 to specify to use 32 cores
  # Switch to --end-to-end not --local
elif [ "$hgRef" == "hg19" ]; then
  # Set bt2idx variable to hg19 location
  bt2idx=/path/to/hg19
  blacklistBed=./hg19-blacklist.bed
  estGenSize=2827436883 # Assuming 150bp read length
  bowtie2 --local --very-sensitive --no-mixed --no-discordant -phred33 -X 700 -p 32 -x $bt2idx -1 "$ID"_1.trimmed.fastq.gz -2 "$ID"_2.trimmed.fastq.gz | samtools view -bS - > "$ID".bam
else
  echo "Error: Unknown reference genome. Please specify either 'hg38' or 'hg19'."
  exit 1
fi

# Sort output .bam file
samtools sort "$ID".bam -o "$ID"_sorted.bam

# Generate an index file
samtools index "$ID"_sorted.bam


### POST-ALIGNMENT QC ----------------------------------------------------------

# Generate idxstats report
samtools idxstats "$ID"_sorted.bam > "$ID"_sorted.idxstats

# Generate flagstat report
samtools flagstat "$ID"_sorted.bam > "$ID"_sorted.flagstat

# Remove mitochondrial reads
samtools view -h "$ID"_sorted.bam | grep -v chrM | samtools sort -O bam -o "$ID".rmChrM.bam -T "$outDir"

# Mark duplicates
module load picard

picard MarkDuplicates QUIET=true INPUT="$ID".rmChrM.bam OUTPUT="$ID".marked.bam METRICS_FILE="$ID".dup.metrics REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR="$outDir"

# Remove duplicates and low quality alignments
# samtools view -q 30 -c "$ID".marked.bam

# Remove multi-mapped reads
# samtools view -h -b -f 2 -F 1548 -q 30 "$ID".marked.bam | samtools sort -o "$ID".filtered.bam

# samtools index "$ID".filtered.bam

# To retain multi-mapped reads:
samtools view -h -b -f 2 -F 1548 "$ID".rmChrM.bam | samtools sort -n -o "$ID".filtered.bam 

# Remove ENCODE blacklist regions
# Remove reads within the blacklist regions
bedtools intersect -nonamecheck -v -abam "$ID".filtered.bam -b $blacklistBed > "$ID".tmp.bam

# Sort and index the bam file
samtools sort -O bam -o "$ID".blacklist-filtered.bam "$ID".tmp.bam
samtools index "$ID".blacklist-filtered.bam

rm "$ID".tmp.bam

module load deeptools

alignmentSieve --numberOfProcessors max --ATACshift --blackListFileName $blacklistBed -b "$ID".blacklist-filtered.bam -o "$ID".shifted.bam


# samtools index <output>.shifted.bam



### BAM VISUALIZATION ----------------------------------------------------------

# 2862010578 is the effective genome size for GRCh38 when using 150bp reads and including only regions which are uniquely mappable
# bam to bigwig
# Set your preferred number of processors
# bamCoverage --numberOfProcessors 8 --binSize 10 --normalizeUsing BPM \
 # --effectiveGenomeSize $estGenSize --bam "$ID".shifted.bam -o "$ID"_coverage_BPM.bw
  

### PEAK CALLING

# Convert the bam file to BEDPE
macs2 randsample -i "$ID".shifted.bam -f BAMPE -p 100 -o "$ID".bed

# Peak calling - MACS2
macs2 callpeak -f BEDPE --nomodel --shift -37 --extsize 73 -g $estGenSize -B --broad --keep-dup all --cutoff-analysis -n "$ID" -t "$ID".bed --outdir macs2/"$ID" 2> macs2.log

# 1. Fold-enrichment bigWig

# 2. -log10 p-value bigWig

# QC - fraction of reads in peak (FRiP score)

# Call peaks for pooled replicates

# Extract replicated peaks

### VISUALIZATION --------------------------------------------------------------

# Convert QC-ed bam file to bedGraph
bedtools bamCoverage --normalizeUsing BPM -b "$ID".filtered.bam > "$ID".bedGraph