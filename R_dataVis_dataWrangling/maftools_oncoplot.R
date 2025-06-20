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
# Load maftools library
library(maftools)

# Clinical data table must have matching Tumor_Sample_Barcode to MAF file
# For this run, the disease stage was used as a binary for
# End stage and newly diagnosed patients
# Table format should be Tumor_Sample_Barcode | Diagnosis_Status
clinical <- read.csv("path/to/metadata.csv", header = TRUE)

clinical$Diagnosis_Status <- factor(clinical$Diagnosis_Status,
                                    levels = c(0, 1),
                                    labels = c("Newly Diagnosed", "End Stage"))

# Import data
maf <- read.maf("/path/to/MAF.txt", clinicalData = clinical)

# OPTIONAL: To run custom oncoplot with more decimal places
source("/path/to/oncoplot_custom.R")

# Plot custom oncoplot
oncoplot_custom(maf=maf, draw_titv = FALSE, showTumorSampleBarcodes = FALSE,
                top = 50, fontSize = 0.6, barcode_mar = 5)