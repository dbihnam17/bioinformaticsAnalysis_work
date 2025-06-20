### This script runs mutational signature analysis using Damiano Fantini's
### algorithm: https://github.com/dami82/mutSignatures
### And uses sample code from his documentation

### Extracts de novo signatures from a group of somatic WGS samples

### I ended up moving to using sigProfiler due to its strong integration
### with the COSMIC cancer signature database v3.4

## Import necessary packages ---------------------------------------------------
library("mutSignatures")
library("dplyr")
library("ggplot2")
library("BSgenome.Hsapiens.UCSC.hg38")


## Data pre-processing ---------------------------------------------------------

# Import VCF file data (as many samples as you would like)
x0 <- importVCFfiles(vcfFiles = c('/path/to/sample1.muTect2.somatic.snv.vcf',
                                  '/path/to/sample2.muTect2.somatic.snv.vcf',
                                  '/path/to/sample3.muTect2.somatic.snv.vcf',
                                  '/path/to/sample_n.muTect2.somatic.snv.vcf'))
                     
# Save copy of original imported file
x <- x0

# Remove INDELs and non-SNV mutations
# Specify 'REF' and 'ALT' columns
x <- filterSNV(dataSet = x, seq_colNames = c("REF", "ALT"))

# Remove mitochondrial data
x <- x %>% 
  filter(CHROM != "chrM")

# OPTIONAL:
# Remove variants with low quality, edit sequence names, adjust identifiers, etc
# Use DPLYR
# x <- x %>% dplyr::filter(QUAL >= 10)

# For muTect2
x <- x %>% dplyr::select(-matches("^G"))

# Retrieve codon at each SNV position
hg38 <- BSgenome.Hsapiens.UCSC.hg38

x <- attachContext(mutData = x, BSGenomeDb = hg38,
                   chr_colName = 'CHROM',
                   start_colName = 'POS', end_colName = 'POS',
                   nucl_contextN = 3, context_colName = 'context')

# Test
# distinctChrom <- x %>% 
  # distinct(CHROM)

# print(distinctChrom)


# Filter to remove SNVs when mismatches are detected between
# Expected and observed DNA sequences (sample vs hg38)
x <- removeMismatchMut(mutData = x, refMut_colName = "REF",
                       context_colName = "context", refMut_format = "N")

# Compute mutation types
x <- attachMutType(mutData = x, ref_colName = "REF",
                   var_colName = "ALT", context_colName = "context")

# Save copy of RevComp1 transformed data in final format
xRCT <- x

# Remove and mutTypes with N
x <- x %>% dplyr::filter(!grepl("N", mutType))

## Extracting mutational signatures --------------------------------------------

# Count three-nucleotide mutation types across samples
blca_counts <- countMutTypes(mutTable = x,
                             mutType_colName = 'mutType',
                             sample_colName = 'SAMPLEID')

# Define parameters for NMF algorithm
blca_params <- setMutClusterParams(num_processesToExtract = 4,
                                   num_totIterations = 50,
                                   num_parallelCores = 5)

# Extract mutational signatures
blca_analysis <- decipherMutationalProcesses(input = blca_counts,
                                             params = blca_params)

## Estimating activity of known mutational signatures --------------------------

# Count three-nucleotide mutations across samples
blca_counts <- countMutTypes(mutTable = x,
                             mutType_colName = 'mutType',
                             sample_colName = 'SAMPLEID')

# Obtain a mutationalSignatures object (de novo/COSMIC)
all_cosmic_v2 <- getCosmicSignatures()

# Select COSMIC mutational signatures
blca_cosmic <- all_cosmic_v2[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                               11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                               21, 22, 23, 24, 25, 26, 27, 28, 29, 30)]

# Resolve mutational signature activities
blca_activ_analysis <- resolveMutSignatures(mutCountData = blca_counts,
                                            signFreqData = blca_cosmic)

## Building plots --------------------------------------------------------------

# Visualize mutational signatures
de_novo_sigs <- blca_analysis$Results$signatures
msigPlot(de_novo_sigs, signature = 1, ylim = c(0, 0.16))
msigPlot(de_novo_sigs, signature = 2, ylim = c(0, 0.16))
msigPlot(de_novo_sigs, signature = 3, ylim = c(0, 0.16))
msigPlot(de_novo_sigs, signature = 4, ylim = c(0, 0.16))
msigPlot(de_novo_sigs, signature = 5, ylim = c(0, 0.16))

# Export a mutationSignatures object as a data.frame
sigs_df <- mutSignatures::as.data.frame(de_novo_sigs)

# Compare mutational signatures
blca_new_vs_cosmic <- matchSignatures(de_novo_sigs, blca_cosmic)

# Show heatmap
blca_new_vs_cosmic$plot

# Visualize signature activities
de_novo_activities <- blca_analysis$Results$exposure

# Set a cutom color palette
my_palette <- c('# e76f51','# f4a261','# e9c46a','# 2a9d8f','# 264653')

# Build plot
msigPlot(de_novo_activities, top = 50)+
  scale_fill_manual(values = my_palette)

# Export a mutSignExposures object as a data.frame
activ_df <- mutSignatures::as.data.frame(de_novo_activities,
                                         transpose = TRUE)

## EXTRA -----------------------------------------------------------------------

sampleCounts <- table(xRCT$SAMPLEID)

sampleCountsSorted <- sort(sampleCounts, decreasing = TRUE)

print(sampleCountsSorted)
