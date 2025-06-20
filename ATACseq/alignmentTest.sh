#!/bin/sh
#SBATCH --partition=cpu-short
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=your.email@email.com
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=2G
#SBATCH -J sample_ATAC_AlignmentTest
#SBATCH --export=USER
#SBATCH --time 4-00:00:00
#SBATCH --chdir=./logs

bt2idx="/path/to/bowtie2Index/hg38"
file1='/path/to/fastq_1.trimmed.fastq.gz'
file2='/path/to/fastq_2.trimmed.fastq.gz'
output='/path/to/outDir'

module load samtools
module load bowtie2

bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 32 -x $bt2idx -1 $file1 -2 $file2 | samtools view -bS - > $output/ZWCT1465.bam
