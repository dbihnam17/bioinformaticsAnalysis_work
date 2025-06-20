### This script is intended to convert very large vcf files into partitioned
### Excel files

### Multiple output files are needed due to the maximum row number constraints
### of a typical Excel file

## Load libraries --------------------------------------------------------------
library(vcfR)
library(openxlsx)

## Data wrangling --------------------------------------------------------------
# Read the VCF file
vcf <- read.vcfR("/path/to/annotated.somatic.variants.vcf.gz")

# Extract the data from the VCF
vcfData <- as.data.frame(vcf@fix)  # Fixed fields
vcfGenotypes <- as.data.frame(vcf@gt)  # Genotype information

# Combine fixed fields and genotype data
vcfCombined <- cbind(vcfData, vcfGenotypes)

# Define the maximum number of rows per Excel file (Excel limit is 1,048,576 rows per sheet)
maxRows <- 1011796  # Use a number under the limit to ensure space for headers

# Calculate the number of Excel files needed
numFiles <- ceiling(nrow(vcfCombined) / maxRows)

## Output Excel files ----------------------------------------------------------
# Define output directory
outputDirectory <- "/path/to/output/"

# Split the data into chunks and write to separate Excel files, including headers
for (i in 1:numFiles) {
  startRow <- ((i - 1) * maxRows) + 1
  endRow <- min(i * maxRows, nrow(vcfCombined))
  
  # Subset the data for the current chunk
  chunkData <- vcfCombined[startRow:endRow, ]
  
  # Create a new Excel file for each chunk
  outputFile <- file.path(outputDirectory, paste0("annotatedVCF_", i, ".xlsx"))
  
  # Write the chunk to an Excel file with headers
  write.xlsx(chunkData, outputFile, rowNames = FALSE, overwrite = TRUE)
}

