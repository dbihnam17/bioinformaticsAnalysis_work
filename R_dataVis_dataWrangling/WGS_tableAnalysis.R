### This Rscript is meant to be used to analyze the <PROJECT> SNV/INDEL table we generated
### from all samples in the project

## Import libraries ------------------------------------------------------------
library(openxlsx)
library(dplyr)
library(ggplot2)

## Data wrangling --------------------------------------------------------------
# Load in data
projectTable <- read.xlsx('/path/to/WGS_mutations_annotation.xlsx')

geneList <- read.xlsx('/path/to/list of genes.xlsx', colNames = FALSE)

projectTableFiltered <- projectTable %>%
  filter(SYMBOL %in% geneList[[1]])

# noncoding <- c('list', 'of', 'noncoding', 'identifiers')

# projectTableFiltered <- ecogTable %>% 
#  filter(Consequence %in% noncoding)

write.xlsx(projectTableFiltered, '/path/to/WGS_drivers.xlsx')

## Prepare noncoding mutations -------------------------------------------------

# Load in noncoding data
noncodingMutations <- read.xlsx('/path/to/WGS_drivers.xlsx')

# Create a new column for unique sample information
noncodingMutations <- noncodingMutations %>%
  mutate(
    Sample_INFO = paste(Sample_ID, "(", Sample_GT, ":", Ref_AD, ",", Alt_AD, ":", Sample_VAF, ")", sep = "")
  )

# Group data by mutation-specific fields and summarize with sample details
topRegions <- noncodingMutations %>%
  group_by(CHROM, POS, REF, ALT, SYMBOL, Consequence, IMPACT, POP_AF,
           CAVA_GENE, CAVA_C, CAVA_P, CAVA_LOC, CAVA_TYPE, CAVA_SO, CAVA_TRANSCRIPT,
           gnomAD_exomes_AF, gnomAD_genomes_AF, gnomAD_genomes_NFE_AF,
           UK10K_AF, ExAC_Adj_AF, ClinVar_CLNSIG, Codons, DISTANCE, DP,
           EXON, Feature, Feature_type, Gene, HGVSc, HGVSp, INTRON, STRAND) %>%
  summarize(
    mutationCount = n(),
    samples = paste(unique(Sample_ID), collapse = ", "),
    Sample_INFO = paste(unique(Sample_INFO), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(mutationCount))

## Export data to Excel --------------------------------------------------------
write.xlsx(topRegions, "/path/to/WGS_driversFin.xlsx")
