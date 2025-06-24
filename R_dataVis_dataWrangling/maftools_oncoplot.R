### This script was originally designed to plot oncoplots from the maftools library
### Other plots have been added, as well as a custom oncoplot

### Only run the first section of the code if your MAF file is missing the 
### End_Position column

## Add End_Position ------------------------------------------------------------
# Load dplyr
library(dplyr)

# Read in the MAF file (tab-delimited)
maf <- read.delim("/path/to/MAF.txt", stringsAsFactors = FALSE)

# Create End_Position column based on Variant_Type
maf <- maf %>%
  mutate(End_Position = case_when(
    Variant_Type == "SNP" ~ Start_Position,  # For SNPs, End_Position = Start_Position
    Variant_Type == "DEL" ~ Start_Position + nchar(Reference_Allele) - 1, # For Deletions, End_Position based on length of deletion
    Variant_Type == "INS" ~ Start_Position + nchar(Reference_Allele) - 1, # For Insertions, End_Position based on insertion length
    TRUE ~ NA_real_  # For other types, set to NA (or handle as needed)
  ))

# Check the updated MAF with End_Position
head(maf)

# OPTIONAL: Save the updated MAF to a new file
write.table(maf, "/path/to/MAF.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)


## Create oncoplot and other base maftools plots -------------------------------
# Load maftools library
library(maftools)

# Import data
maf <- read.maf("/path/to/MAF.txt")

# To import a maf file with non-standard variant classifications
# You must list ALL classifications you want to see in vc_nonSyn
# not just the new one you are adding
#maf <- read.maf("/path/to/MAF.txt", 
#                vc_nonSyn = c("Nonstop_Mutation", "Frame_Shift_Del", "Missense_Mutation", 
#                              "Nonsense_Mutation", "Splice_Site", "Frame_Shift_Ins", 
#                              "In_Frame_Del", "In_Frame_Ins", "Start_Codon_Del", "In_Frame_Indel",
#                              "etc..."))

# Set font and style
par(family = "Calibri", font = 1)

# Plot oncoplot
oncoplot(maf=maf, draw_titv = TRUE, showTumorSampleBarcodes = TRUE,
         top = 10000, fontSize = 0.7, barcode_mar = 5)

# Plot co-expression matrix
somaticInteractions(maf = maf, top = 10000,
                    pvalue = c(0.01, 0.05), fontSize = 0.5)

# Plot maf summary
plotmafSummary(maf = maf, rmOutlier = FALSE, addStat = 'mean',
               dashboard = TRUE, titvRaw = FALSE)

# Plot pathways
pathways(maf = maf, plotType = 'treemap')

# OPTIONAL: To run custom oncoplot with more decimal places
# To create oncoplot_custom.R, I extracted the source code for the oncoplot function
# and modified the original code to allow for more decimal places to show
# on the right side of the plot
source("/path/to/oncoplot_custom.R")
# Plot custom oncoplot
oncoplot_custom(maf=maf, draw_titv = FALSE, showTumorSampleBarcodes = FALSE,
                top = 10000, fontSize = 0.6, barcode_mar = 5)


## Create oncoplot with patient metadata ---------------------------------------
# Load in maftools if not already loaded
library(maftools)

# Clinical data table must have matching Tumor_Sample_Barcode to MAF file
# For this run, the disease stage was used as a binary for
# End stage and newly diagnosed patients
# Table format should be Tumor_Sample_Barcode | Diagnosis_Status
clinical <- read.csv("path/to/metadata.csv", stringsAsFactors = FALSE)


# Import patient metadata
clinical <- read.csv("/path/to/metadata.csv",
                     stringsAsFactors = FALSE)

# Load MAF as a dataframe for filtering
# This maf included CNV data derived from FISH entered in as "mutations"
# in order to structure plot nicely
maf <- read.delim("path/to/MAF.txt", stringsAsFactors = FALSE)

# Filter MAF rows to keep only samples present in clinical metadata
maf_filtered <- maf[maf$Tumor_Sample_Barcode %in% clinical$Tumor_Sample_Barcode, ]

# Add filler data to CNV "mutations" for proper plotting
maf_filtered$Chromosome[maf_filtered$Variant_Classification == "Translocation"] <- "chrUn"
maf_filtered$Start_Position[maf_filtered$Variant_Classification == "Translocation"] <- 0
maf_filtered$End_Position[maf_filtered$Variant_Classification == "Translocation"] <- 0
maf_filtered$Reference_Allele[maf_filtered$Variant_Classification == "Translocation"] <- "-"
maf_filtered$Tumor_Seq_Allele2[maf_filtered$Variant_Classification == "Translocation"] <- "-"


# Load MAF object and associate clinical metadata
maf <- read.maf(maf = maf_filtered,
                clinicalData = clinical,
                # Use custom variant classification list
                vc_nonSyn = c("Nonstop_Mutation", "Frame_Shift_Del", "Missense_Mutation", 
                              "Nonsense_Mutation", "Splice_Site", "Frame_Shift_Ins", 
                              "In_Frame_Del", "In_Frame_Ins", "Start_Codon_Del", "In_Frame_Indel",
                              "Nonstart_Mutation", "Gain", "Deletion", "Translocation",
                              "Hyperdiploid_Positive"))

# Multiple Myeloma driver gene list from previous studies
gene_list <- c(
  "ACTG1", "AKT1", "AKT2", "AKT3", "ATF4", "ATF6", "ATM", "ATR", "B2M", "BIRC2",
  "BIRC3", "BRAF", "BTG1", "CARD11", "CCND1", "CCNT1", "CD38", "CDK4", "CDK7",
  "CDKN1B", "CDKN2A", "CDKN2C", "DDIT3", "CRBN", "CUL4A", "CUL4B", "CXCR4",
  "CYLD", "DIS3", "DUSP2", "EGFR", "EGR1", "ERN1", "FAM46C", "FGFR3", "GRB2",
  "HSPA5", "IDH1", "IDH2", "IDH3A", "IFNGR2", "IGF1R", "IKZF1", "IKZF3", "IL6",
  "IL6R", "IL6ST", "IRF4", "JAK2", "KDM6A", "KRAS", "MAF", "MAFB", "MAX", "MYC",
  "MYD88", "NFKB2", "NFKBIA", "NFKBIB", "NR3C1", "NRAS", "EIF2AK3", "PIK3CA",
  "PIK3CG", "PIK3R1", "PIK3R2", "PIM1", "PIM2", "PIM3", "PRDM1", "PSMA1", "PSMA2",
  "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7", "PSMA8", "PSMB1", "PSMB10", "PSMB11",
  "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6", "PSMB7", "PSMB8", "PSMB9", "PSMC1",
  "PSMC2", "PSMC3", "PSMC3IP", "PSMC4", "PSMC5", "PSMC6", "PSMD1", "PSMD10",
  "PSMD11", "PSMD12", "PSMD13", "PSMD14", "PSMD2", "PSMD3", "PSMD4", "PSMD5",
  "PSMD6", "PSMD7", "PSMD8", "PSMD9", "PSME1", "PSME2", "PSME3", "PSME4", "PSMF1",
  "PSMG1", "PSMG2", "PSMG3", "PSMG4", "PTPN11", "RASA2", "RB1", "RIPK1", "RIPK4",
  "SHC1", "SP140", "STAT3", "TET2", "TGFBR2", "TLR4", "TNFRSF13B", "TNFRSF21",
  "TP53", "TRAF2", "TRAF3", "TRAF3IP1", "WHSC1", "XBP1", "ZFHX4"
)

# Get gene summary table from maf object
gene_summary <- getGeneSummary(maf)

# Filter gene summary to only genes of interest
gene_summary_filtered <- gene_summary[Hugo_Symbol %in% gene_list]

# Order by mutation count descending
gene_summary_filtered <- gene_summary_filtered[order(-MutatedSamples)]

# Take top 'n' genes from that filtered list
top30_genes <- head(gene_summary_filtered$Hugo_Symbol, 30)

# Insert CNV mutations manually above 'n' genes
top30_genes_FISH <- c("t11_14", "t4_14", "t14_20", "t14_16", "t6_14",
                      "Hyperdiploid", "1q", "17p", top30_genes)

# Plot with Arial Narrow font (Save space for large plot)
par(family = "Arial Narrow", font = 1)

# Alter plot margins (bottom, left, top, right)
par(mar = c(5, 5, 5, 5))

# Plot oncoplot with patient data and CNVs from FISH
# oncoplot_custom.R does not work here
# oncoplot(clinicalFeatures) is dependent on other functions from maftools
# This would require much more modification of maftools functions to implement
oncoplot(maf = maf,
         genes = top30_genes_FISH,
         sortByAnnotation = TRUE,
         #draw_titv = TRUE,
         fontSize = 0.75,
         legend_height = 5,
         annotationFontSize = 1.2,
         keepGeneOrder = TRUE,
         legendFontSize = 1.2,
         showTumorSampleBarcodes = FALSE,
         clinicalFeatures = "Disease.Stage")

