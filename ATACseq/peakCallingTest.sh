#!/bin/sh
#SBATCH --partition=cpu-short
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=your.email@email.com
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=2G
#SBATCH -J ATAC_peakCallingTest
#SBATCH --export=USER
#SBATCH --time 4-00:00:00

ID='<sample_ID>'
estGenSize=2862010428

#Convert the bam file to BEDPE
/path/to/macs2/2.2.9.1/bin/macs2 randsample -i "$ID".blacklist-filtered.bam -f BAMPE -p 100 -o "$ID".bed

# Peak calling - MACS2
/path/to/macs2/2.2.9.1/bin/macs2 callpeak -f BEDPE --nomodel --shift -37 --extsize 73 -g $estGenSize -B --broad --keep-dup all --cutoff-analysis -n "$ID" -t "$ID".bed --outdir macs2/"$ID" 2> macs2.log
