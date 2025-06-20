#!/bin/bash
#SBATCH --partition=cpu-short
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=bihnam.daniel@mayo.edu
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=2G
#SBATCH -J G_2287v2_bam_flagstat
#SBATCH --export=USER
#SBATCH --time 0-10:00:00
#SBATCH --chdir=./logs

module load samtools/1.18

# Make flagstat report
samtools flagstat /research/bsi/archive/PI/Braggio_Esteban_m037525/secondary/Michelson/DNA/H202SC24051873/H202SC24051873.2/02.Bam/G_2287v2.final.bam > /research/labs/hematology/braggio-slager/m293215/G_2287v2.flagstat
