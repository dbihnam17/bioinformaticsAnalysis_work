#!/bin/bash
#SBATCH --partition=cpu-short
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=bihnam.daniel@mayo.edu
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=2G
#SBATCH -J mich2_bam_flagstat
#SBATCH --export=USER
#SBATCH --time 0-10:00:00
#SBATCH --chdir=./logs

module load samtools/1.18

# Make flagstat report
samtools flagstat /research/bsi/archive/PI/Braggio_Esteban_m037525/secondary/Michelson/DNA/H202SC24051873/02.Bam/D_1488_0518.final.bam > /research/labs/hematology/braggio-slager/m293215/Michelson/batch2/D_1488_0518/D_1488_0518.flagstat

samtools flagstat /research/bsi/archive/PI/Braggio_Esteban_m037525/secondary/Michelson/DNA/H202SC24051873/02.Bam/G_1488.final.bam > /research/labs/hematology/braggio-slager/m293215/Michelson/batch2/G_1488/G_1488.flagstat

samtools flagstat /research/bsi/archive/PI/Braggio_Esteban_m037525/secondary/Michelson/DNA/H202SC24051873/02.Bam/D_2287_1215.final.bam > /research/labs/hematology/braggio-slager/m293215/Michelson/batch2/D_2287_1215/D_2287_1215.flagstat

samtools flagstat /research/bsi/archive/PI/Braggio_Esteban_m037525/secondary/Michelson/DNA/H202SC24051873/02.Bam/G_2287.final.bam > /research/labs/hematology/braggio-slager/m293215/Michelson/batch2/G_2287/G_2287.flagstat