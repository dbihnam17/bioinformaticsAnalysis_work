### This function will be used to filter tabular data on various conditions defined
### in individual functions, then save them as .xlsx files

### This function is designed in a way in which new filters can be added if needed
### Just add another argument for the user to supply the desired filtering function

### This was originally designed to just use excel files as inputs and outputs,
### but was modified to accept common alternate inputs

### For now, I will leave it only writing .xlsx outputs, but the script can easily
### be modified to include an extra user argument to specify output file type
### (Could be an optional flag...)

### >source("/path/to/filterTable.R")
### >filterTable(help) for usage and current filters

### NOTE: filePaths can accept individual file paths, directories, and glob patterns

## Load libraries --------------------------------------------------------------

library(readxl)
library(writexl)
library(tools)

## filterTable function --------------------------------------------------------

## This function handles the overall structure of this code
## It will call other filtering functions in this script
## For now there is only svFilt

filterTable <- function(filePaths, outDir, filtFunction) {
  
  # Handle confused users
  nArgs <- length(as.list(match.call())[-1])
  if ((nArgs != 3) || (length(filePaths) == 1 && tolower(filePaths) %in% c("help", "usage"))) {
    helpFunction()
    return(invisible(NULL))
  }
  
  # Handle wildcard paths or directory
  if (length(filePaths) == 1 && dir.exists(filePaths)) {
    filePaths <- list.files(filePaths, 
                            pattern = "\\.(xls|xlsx|csv|tsv|txt)$", 
                            full.names = TRUE, ignore.case = TRUE)
  } else if (length(filePaths) == 1 && grepl("[*?]", filePaths)) {
    filePaths <- Sys.glob(filePaths)
  }
  
  # Check files exist
  if (length(filePaths) == 0) {
    stop("No files found. Check your file paths or wildcard pattern.")
  }
  
  # Create output directory if needed
  if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = TRUE)
  }
  
  # Loop through files
  for (file in filePaths) {
    ext <- tolower(tools::file_ext(file))
    data <- NULL
    
    # Read based on file extension
    if (ext %in% c("xls", "xlsx")) {
      data <- tryCatch({
        readxl::read_excel(file)
      }, error = function(e) {
        message(sprintf("Warning: Failed to read %s as Excel, trying as TSV...", file))
        tryCatch(read.delim(file, stringsAsFactors = FALSE), error = function(e2) {
          stop(sprintf("Failed to read file %s as Excel or TSV.", file))
        })
      })
    
    # Handle .tsv/.txt files    
    } else if (ext %in% c("tsv", "txt")) {
      data <- tryCatch(read.delim(file, stringsAsFactors = FALSE), 
                       error = function(e) stop(sprintf("Failed to read TSV file %s", file)))
    # Handle .csv files
    } else if (ext == "csv") {
      data <- tryCatch(read.csv(file, stringsAsFactors = FALSE), 
                       error = function(e) stop(sprintf("Failed to read CSV file %s", file)))
    
    # Reject other file extensions    
    } else {
      stop(sprintf("Unsupported file extension '%s' for file %s", ext, file))
    }
    
    # Apply the filter function
    filtData <- filtFunction(data)
    
    # Truncate any overly long strings to avoid Excel errors
    filtData <- truncateLongStr(filtData)
    
    # Save as .xlsx always
    fileName <- tools::file_path_sans_ext(basename(file))
    outFile <- file.path(outDir, paste0(fileName, "_filtered.xlsx"))
    writexl::write_xlsx(filtData, outFile)
    
    message(sprintf("Filtered file saved to: %s", outFile))
  }
  
  invisible(NULL)
}

## truncateLongStr function ----------------------------------------------------

## This function will handle removing strings longer than 32,767 character
## to allow file to be saved as a .xlsx

truncateLongStr <- function(data, max_length = 32767) {
  data[] <- lapply(data, function(col) {
    if (is.character(col)) {
      # Truncate strings longer than max_length
      sapply(col, function(x) {
        if (nchar(x, type = "chars") > max_length) {
          substr(x, 1, max_length)
        } else {
          x
        }
      }, USE.NAMES = FALSE)
    } else {
      col
    }
  })
  return(data)
}

## filtSV function -------------------------------------------------------------

filtSV <- function(data) {
  # Filter for only translocations and inversions
  subset(data, SVType %in% c("INV", "BND"))
}

## Filters available -----------------------------------------------------------

functionDict <- list(
  filtSV = list(
    desc = "Function to return only translocations (BND) and inversions (INV) from delly SVType column in somatic structural variant files.",
    func = filtSV
  )
)

## helpFunction function -------------------------------------------------------

## This function handles invalid inputs and informs the user of proper usage

helpFunction <- function() {
  
  # Display generic function usage
  cat("\n--- filterTable() Usage ---\n\n")
  cat("filterTable(filePaths, outDir, filtFunction)\n\n")
  cat("Arguments:\n")
  cat("  - filePaths: Character vector of Excel file paths to process\n")
  cat("  - outDir: Directory path where filtered files will be saved\n")
  cat("  - filtFunction: Filtering function to apply on each file's data\n\n")
  
  # Describe file input flexibility
  cat("File input options:\n")
  cat("  - Accepts both .xls and .xlsx formats\n")
  cat("  - Can pass individual file paths: c(\"file1.xlsx\", \"file2.xls\")\n")
  cat("  - Can pass a directory containing Excel files\n")
  cat("  - Can use glob patterns: \"path/to/files/*.xls*\"\n\n")
  
  # Example usage
  cat("\nExample:\n")
  cat('  filterTable(\n')
  cat('    filePaths = c("file1.xlsx", "file2.xlsx"),\n')
  cat('    outDir = "filtered_output",\n')
  cat('    filtFunction = filtSV\n')
  cat('  )\n\n')
  
  # Print out available functions and their description
  cat("Available filtering functions:\n")
  for (fname in names(functionDict)) {
    desc <- functionDict[[fname]]$desc
    cat(sprintf("  - %s: %s\n", fname, desc))
  }
  
  # Exit cleanly
  invisible(NULL)
}
