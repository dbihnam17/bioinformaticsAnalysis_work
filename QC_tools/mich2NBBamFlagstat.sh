#!/bin/bash
#SBATCH --partition=cpu-short
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=bihnam.daniel@mayo.edu
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=2G
#SBATCH -J mich2NB_bam_flagstat
#SBATCH --export=USER
#SBATCH --time 0-10:00:00
#SBATCH --chdir=./logs

module load samtools/1.18

# Make flagstat report
samtools flagstat /research/bsi/archive/PI/Braggio_Esteban_m037525/secondary/Michelson/DNA/H202SC24051873/H202SC24051873.4/02.Bam/D1488_0518N6.final.bam > /research/labs/hematology/braggio-slager/m293215/Michelson/testSamples/03.Result_X202SC24051873-Z01-F010_cancer/flagstat/D1488_0518N6.flagstat

samtools flagstat /research/bsi/archive/PI/Braggio_Esteban_m037525/secondary/Michelson/DNA/H202SC24051873/H202SC24051873.4/02.Bam/G1488NEB6.final.bam > /research/labs/hematology/braggio-slager/m293215/Michelson/testSamples/03.Result_X202SC24051873-Z01-F010_cancer/flagstat/G1488NEB6.flagstat

samtools flagstat /research/bsi/archive/PI/Braggio_Esteban_m037525/secondary/Michelson/DNA/H202SC24051873/H202SC24051873.4/02.Bam/D2287_1215N6.final.bam > /research/labs/hematology/braggio-slager/m293215/Michelson/testSamples/03.Result_X202SC24051873-Z01-F010_cancer/flagstat/D2287_1215N6.flagstat

samtools flagstat /research/bsi/archive/PI/Braggio_Esteban_m037525/secondary/Michelson/DNA/H202SC24051873/H202SC24051873.4/02.Bam/G2287NEB6.final.bam > /research/labs/hematology/braggio-slager/m293215/Michelson/testSamples/03.Result_X202SC24051873-Z01-F010_cancer/flagstat/G2287NEB6.flagstat