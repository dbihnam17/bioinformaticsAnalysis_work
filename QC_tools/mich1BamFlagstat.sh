#!/bin/bash
#SBATCH --partition=cpu-short
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=bihnam.daniel@mayo.edu
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=2G
#SBATCH -J mich1_bam_flagstat
#SBATCH --export=USER
#SBATCH --time 0-10:00:00
#SBATCH --chdir=./logs

module load samtools/1.18

# Make flagstat report
samtools flagstat /research/bsi/archive/PI/Braggio_Esteban_m037525/secondary/Michelson/DNA/H202SC23053360/02.Bam/D_3177_0917.final.bam > /research/labs/hematology/braggio-slager/m293215/D_3177_0917.flagstat

samtools flagstat /research/bsi/archive/PI/Braggio_Esteban_m037525/secondary/Michelson/DNA/H202SC23060010/02.Bam/D_3089_0715.final.bam > /research/labs/hematology/braggio-slager/m293215/D_3089_0715.flagstat
