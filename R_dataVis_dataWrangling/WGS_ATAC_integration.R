### This script was meant to intersect a subset of our WGS somatic
### mutation data with the ATAC peaks we had identified

## Load libraries --------------------------------------------------------------
library(GenomicRanges)
library(readxl)
library(openxlsx)

## Data wrangling --------------------------------------------------------------
# Load in annotated vcf
PROJECTTable <- read_excel("/path/to/PROJECT_noncodingAnalysis.xlsx",
                        sheet = "Gene Subset")

# Create a GRanges object for the VCF mutations
masterGR <- GRanges(
  seqnames = PROJECTTable$CHROM,
  ranges = IRanges(start = PROJECTTable$POS, end = PROJECTTable$POS)
)

# Load in ATAC peaks master BED
masterBED <- read.csv("/path/to/ATAC_masterPeak.bed", header = FALSE)
colnames(masterBED) <- c("CHROM", "start", "end")
masterGrBED <- GRanges(
  seqnames = masterBED$CHROM,
  ranges = IRanges(start = masterBED$start, end = masterBED$end)
)

# Load the up bed file
upBED <- read.table("/path/to/ATAC_up_bed.bed", header = FALSE)
colnames(upBED) <- c("CHROM", "start", "end")
upGrBED <- GRanges(
  seqnames = upBED$CHROM,
  ranges = IRanges(start = upBED$start, end = upBED$end)
)

# Load the down bed file
downBED <- read.table("/path/to/ATAC_down_bed.bed", header = FALSE)
colnames(downBED) <- c("CHROM", "start", "end")
downGrBED <- GRanges(
  seqnames = downBED$CHROM,
  ranges = IRanges(start = downBED$start, end = downBED$end)
)


# Initialize columns in main table
PROJECTTable$MasterPeak <- "."
PROJECTTable$UpPeak <- "."
PROJECTTable$DownPeak <- "."

# Find overlaps and update columns
PROJECTTable$MasterPeak[queryHits(findOverlaps(masterGR, masterGrBED))] <- "Chromatin open"
PROJECTTable$UpPeak[queryHits(findOverlaps(masterGR, upGrBED))] <- "Upregulated"
PROJECTTable$DownPeak[queryHits(findOverlaps(masterGR, downGrBED))] <- "Downregulated"

## Export data -----------------------------------------------------------------
write.xlsx(PROJECTTable, "/path/to/PROJECT_geneSubsetNoncodingATAC.xlsx", row.names = FALSE)
