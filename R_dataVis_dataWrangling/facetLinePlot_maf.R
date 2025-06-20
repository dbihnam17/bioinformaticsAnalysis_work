### This script is intended to use a tab-delimited txt file in to plot 
### mutation prevalence across a patient cohort

### Does not need to be in maf format, it is just the data format I had

## Necessary columns: |Hugo_Symbol|Timepoint_TMP|Tumor_Sample_Barcode|
## Optional columns: |Start_Position| (I used this to remove some NA rows in my maf file)

### This is being used to characterize the mutational landscape of CLL patients
### throughout their treatment course on BTKi drugs (ibrutinib, acalbrutinib, etc)

## Load libraries --------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

## Data wrangling --------------------------------------------------------------
# Read in data
maf <- read.table("/path/to/MAF.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Filter to samples with all 3 timepoints --Used in a previous version--
# filtered_maf <- maf %>%
#   filter(!is.na(Start_Position)) %>%
#   mutate(Sample_ID = str_sub(Tumor_Sample_Barcode, 1, 4)) %>%
#   group_by(Sample_ID) %>%
#   filter(all(c("Baseline", "On Treatment", "Progression") %in% Timepoint_TMP)) %>%
#   ungroup() %>%
#   select(-Sample_ID)  

# Filter to samples with at least 2 distinct timepoints
filtered_maf <- maf %>%
  filter(Timepoint_TMP != "Post Progression") %>% # Exclude post-progression timepoint
  filter(!is.na(Start_Position)) %>% # Remove NA rows
  mutate(Sample_ID = str_sub(Tumor_Sample_Barcode, 1, 4)) %>% # Extract the first 4 characters of Tumor_Sample_Barcode as Sample_ID
  group_by(Sample_ID) %>% # Group by Sample_ID
  filter(n_distinct(Timepoint_TMP) >= 2) %>% # Keep only samples with at least 2 distinct timepoints
  ungroup() %>%
  select(-Sample_ID) #Drop the Sample_ID column


# Count mutations per gene and timepoint
mutation_counts <- filtered_maf %>%
  filter(!is.na(Hugo_Symbol) & Hugo_Symbol != "") %>%
  group_by(Hugo_Symbol, Timepoint_TMP) %>%
  summarise(Mutation_Count = n(), .groups = "drop")

# Calculate total mutations per gene
gene_order <- mutation_counts %>%
  group_by(Hugo_Symbol) %>%
  summarise(Total_Mutations = sum(Mutation_Count)) %>%
  arrange(desc(Total_Mutations))

# Filter for top n genes
# For CLL BTKi study, PLCG2 should be shown as well as BTK to show drug resistance
top_genes <- gene_order %>%
  top_n(50, Total_Mutations) %>% # Get top n genes by total mutation count
  # bind_rows(filter(gene_order, Hugo_Symbol == "PLCG2")) %>%  #Add PLCG2 if missing
  distinct() %>% # Remove duplicates if PLCG2 was already in the top n
  pull(Hugo_Symbol) # Get gene names

# Filter data for the top n genes
filtered_top12 <- mutation_counts %>% filter(Hugo_Symbol %in% top_genes)

# Fill in empty timepoints for genes
complete_data <- filtered_top12 %>%
  complete(Hugo_Symbol, Timepoint_TMP = c("Baseline", "On Treatment", "Progression"), fill = list(Mutation_Count = 0))

## Visualization ---------------------------------------------------------------
# Plot faceted line graph for number of genes selected
ggplot(complete_data, aes(x = Timepoint_TMP, y = Mutation_Count, group = Hugo_Symbol, color = Hugo_Symbol)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  facet_wrap(~ factor(Hugo_Symbol, levels = top_genes), scales = "fixed") +
  theme_minimal() +
  labs(
    title = "Mutation Trends for Top 36 Genes",
    x = "Timepoint",
    y = "Mutation Count",
    color = "Gene") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(family = "Calibri"),
    legend.position = "none")


