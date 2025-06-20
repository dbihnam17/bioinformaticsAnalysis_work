### This script was used to create a graphical heatmap from an excel table containing
### counts of mutations per sample in genomic hotspots
### There was also a barchart on the right representing the prevalence of each
### mutation

### This is basically just a maftools oncoplot before I learned that you could
### use maftools instead...

### The imported csv was a curated matrix with the count (if any) of mutations present
### in that sample

### The data was collected by manually inspecting bam files in IGV and recording
### details about any mutations seen in the specific position

## Import libraries ------------------------------------------------------------
library(dplyr)
library(ComplexHeatmap)
library(circlize)

## Data wrangling --------------------------------------------------------------
# Import data
hotspotTable <- read.csv('/path/to/hotspotTableHeatmap.csv')

# Modify column names
colnames(hotspotTable) <- gsub("^X", "", colnames(hotspotTable))

# Ensure that the 'Gene' column is the rownames and separate the 'total' column
rownames(hotspotTable) <- hotspotTable$Gene # Set rownames to the 'Gene' column
hotspotTable$Gene <- NULL # Remove the 'Gene' column
totalCounts <- hotspotTable$total # Extract the 'total' column
hotspotTable$total <- NULL # Remove the 'total' column from the main dataframe

# Create groups based on the first four characters of the sample names
columnGroups <- substr(colnames(hotspotTable), 1, 2)

# Define discrete colors for each mutation count value
mutationLevels <- c(0, 1, 2)
discreteColors <- c("grey", "#3f88c5", "#ffba08")
names(discreteColors) <- mutationLevels

# Use a discrete color function to map each mutation count to its corresponding color
heatmapColors <- function(x) discreteColors[as.character(x)]

# Define gene grouping
geneGroups <- factor(1:nrow(hotspotTable))

## Plot the heatmap and barchart -----------------------------------------------
# Create the heatmap with row division
heatmap <- Heatmap(
  as.matrix(hotspotTable),
  name = "Mutations",
  col = heatmapColors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE, # Remove sample IDs from the bottom
  show_row_names = TRUE,
  column_split = columnGroups,
  row_split = geneGroups,
  row_title = "Genes",
  column_title = "Mutational Distribution of Key Genes Across Patients",
  row_names_side = "left",
  rect_gp = gpar(col = "white"),
  
  # Make titles and labels bold
  column_title_gp = gpar(fontface = "bold", fontsize = 12), # Column title (top)
  row_title_gp = gpar(fontface = "bold", fontsize = 12), # Row title ("Genes" label)
  row_names_gp = gpar(fontface = "bold", fontsize = 10), # Gene row names on the left
  column_names_gp = gpar(fontface = "bold", fontsize = 8), # Column names (patient IDs)
  
  # Make heatmap legend labels bold
  heatmap_legend_param = list(
    title = "Mutations",
    at = c(0, 1, 2), # Display legend values for mutation counts
    labels = c("0", "1", "2"), # Labels corresponding to the counts
    legend_height = unit(4, "cm"), # Height of the vertical legend
    direction = "vertical", # Set legend direction to vertical
    title_gp = gpar(fontface = "bold"), # Bold font for legend title
    labels_gp = gpar(fontface = "bold") # Bold font for legend labels
  ),
  
  # Set spacing between rows and columns
  column_gap = unit(4, "mm"), # Increase the horizontal gap between column groups
  row_gap = unit(3, "mm") # Increase the vertical gap between gene rows
)

# Create a bar plot for the 'total' column as a heatmap annotation
# Calculate percentages for the 'totalCounts' column
totalCounts_percentage <- (totalCounts / 16) * 100

total_barplot <- rowAnnotation(
  "Total Mutations" = anno_barplot(
    totalCounts_percentage, 
    border = FALSE, 
    gp = gpar(fill = "#3f88c5", col = "#3f88c5"), # Bar color
    axis_param = list(gp = gpar(fontface = "bold")), # Bold font for barplot axis labels
    axis = TRUE,
  ),
  width = unit(4, "cm"),
  
  # Set the label parameters for the row annotation (the title of the barplot)
  annotation_label = "Patients (%)",
  annotation_name_gp = gpar(fontface = "bold", fontsize = 10) # Make the label bold
)

# Combine and draw the heatmap with the annotation
# Draw the heatmap with the barplot and set the legend position to the left
draw(heatmap + total_barplot, 
     heatmap_legend_side = "left",
     annotation_legend_side = "left")
