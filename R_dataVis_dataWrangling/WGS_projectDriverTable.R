### This script is meant to join all the Michelson 1 and 2 datq
### We will also be filtering all SNVs and INDELs for only the 61 CLL genes of interest

## Load necessary packages -----------------------------------------------------
library(readxl)
library(openxlsx)
library(dplyr)

## Data wrangling --------------------------------------------------------------
# Save Somatic_SNV and Somatic_INDEL paths (files were stored across multiple sequencing deliveries)
dataDir <- "/path/to/Somatic_SNV/Annotation"

fileList <- list.files(dataDir, pattern = "\\.xls$", full.names = TRUE)

# print(fileList)

# Assign drivers variable with list of 61 relevant driver genes
# drivers <- read.csv("/path/to/61CLLgenes", header = FALSE)

# combinedData <- data.frame()

# Rename the column in the drivers table
# colnames(drivers) <- "GeneName"

# Loop to combine files from each directory
for (file in fileList) {
  # Read the file
  data <- read.delim(file, header = TRUE, sep = "\t")
  
  # Swap the column names at positions 83 and 84 *BATCH 2 ONLY*
  # colnames(data)[c(83, 84)] <- colnames(data)[c(84, 83)]
  
  # Filter to include only rows where GeneName is in the drivers table
  # filteredData <- data %>% filter(GeneName %in% drivers$GeneName)
  
  # Extract the sample ID from the file name (first 3-11 characters)
  sampleId <- substr(basename(file), 3, 11)
  
  # Add the sample ID as a new column filteredData WHEN USING FILTER
  filteredData <- data %>% mutate(sampleID = sampleId)
  
  # Remove noncoding and synonymous SNV
  filteredData <- filteredData %>% filter(Func %in% c("exonic", "UTR3", "UTR5", "splicing"))
  
  # Append the filtered table to combinedData
  combinedData <- bind_rows(combinedData, filteredData)
}

# Check the combined data
print(dim(combinedData))
head(combinedData)
tail(combinedData)

## Output results --------------------------------------------------------------
# Write <project>DriverSNV and <project>DriverINDEL
write.xlsx(combinedData, file = '/path/to/projectDriverSNV.xlsx')