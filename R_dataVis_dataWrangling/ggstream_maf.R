### This script is intended to use a tab-delimited txt file in to plot 
### tumor mutational burden across a patient cohort

### Does not need to be in maf format, it is just the data format I had

## Necessary columns: |Hugo_Symbol|Timepoint_TMP|
## Optional columns: |Start_Position| (I used this to remove some NA rows in my maf file)

### This is being used to characterize the mutational landscape of CLL patients
### throughout their treatment course on BTKi drugs (ibrutinib, acalbrutinib, etc)

## Load required libraries -----------------------------------------------------
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggstream)

## Data wrangling --------------------------------------------------------------

# Read in data (maf format not required, just what I had)
maf <- read.table("/path/to/MAF.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Remove NA rows
maf <- maf %>% filter(!is.na(Start_Position))

# View the top n genes
geneSummary <- maf %>%
  group_by(Hugo_Symbol) %>%
  summarise(Mutation_Count = n(), .groups = "drop") %>% 
  arrange(desc(Mutation_Count)) %>% 
  print(n = 20)

# Extract a vector of the top n genes (change to include more explicit genes)
# NOTE: Color palette was designed for 11 genes + other (12 total colors)
topGenes <- geneSummary %>% 
  head(11) %>% 
  pull(Hugo_Symbol)

# Add a new column to bin all other mutations into 'Other'
mafOther <- maf %>%
  mutate(Hugo_Symbol = ifelse(Hugo_Symbol %in% topGenes, Hugo_Symbol, "OTHER"))

# Filter out the "Other" category for the proportions calculation
proportionDF <- mafOther %>%
  group_by(Timepoint_TMP, Hugo_Symbol) %>%
  summarise(Mutation_Count = n(), .groups = "drop") %>%  # Count mutations per gene per timepoint
  group_by(Timepoint_TMP) %>%
  mutate(Timepoint_Total = sum(Mutation_Count),  # Compute total mutations per timepoint
         Gene_Proportion = Mutation_Count / Timepoint_Total) %>%  # Calculate proportions
  select(Hugo_Symbol, Timepoint_TMP, Gene_Proportion)

# Recalculate proportion sums
proportionSum <- proportionDF %>%
  group_by(Timepoint_TMP) %>%
  summarise(Total_Proportion = sum(Gene_Proportion))

# Ensure proportions add up to 1 for all timepoints
print(proportionSum)

# Make Timepoint_TMP a factor
# Put timepoints in the right order and remove post-progression
# NOTE: This code can be adjusted for more timepoints, post-progression is being
#       excluded due to low sample size, and only have TP53 mutations
proportionDF$Timepoint_TMP <- factor(proportionDF$Timepoint_TMP, 
                                         levels = c("Baseline", "On Treatment", "Progression"))

# Remove NA rows if present
proportionDF <- proportionDF %>% filter(!is.na(Timepoint_TMP))

# Convert Timepoint_TMP to numeric
proportionDF$Timepoint_Num <- as.numeric(proportionDF$Timepoint_TMP)

## Plotting --------------------------------------------------------------------

# Create color palette (Will need to be adjusted if using more than 12 groups)
custom_colors <- c(
  "#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", 
  "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", 
  "#5e4fa2", "navy")

# Manually adjust BTK and PLCG2 mutations so they visually start at TP2
# When just set to Timepoint_Num == 2, ggstream tries to make the data smooth
# And starts them at a low level around ~TP '1.7' which does not happen biologically
# Handle this on a per-case basis, the plot should match the data
proportionDF <- proportionDF %>%
  mutate(Timepoint_Num = ifelse(Hugo_Symbol %in% c("BTK", "PLCG2") & Timepoint_Num == 2.0, 2.75, Timepoint_Num))

# Reorder Hugo_Symbol to move the OTHER category to the bottom of the plot
proportionDF$Hugo_Symbol <- factor(proportionDF$Hugo_Symbol, 
                                       levels = c(topGenes, "OTHER"))

# Plot streamplot with OTHER
ggplot(proportionDF, aes(x = Timepoint_Num, y = Gene_Proportion, fill = Hugo_Symbol)) +
  geom_stream(type = "proportional", bw = 1.27) +
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(
    breaks = c(1, 2, 3), 
    labels = c("Baseline", "On-Treatment", "Progression")) +
  labs(title = "Tumor Mutational Burden of CLL Patients Throughout Treatment",
       x = "Timepoint",
       y = "TMB",
       fill = "Mutation") +
  theme_minimal() +
  theme(axis.text.x = element_text(),
        text = element_text(family = "Calibri"))

# Filter out OTHER
proportionDF_noOther <- proportionDF %>%
  filter(Hugo_Symbol != "OTHER")

# Plot streamplot without OTHER
ggplot(proportionDF_noOther, aes(x = Timepoint_Num, y = Gene_Proportion, fill = Hugo_Symbol)) +
  geom_stream(type = "proportional", bw = 1.27) +
  scale_fill_manual(values = custom_colors) + 
  scale_x_continuous(
    breaks = c(1, 2, 3), 
    labels = c("Baseline", "On-Treatment", "Progression")) +
  labs(title = "Tumor Mutational Burden of CLL Patients Throughout Treatment",
       x = "Timepoint",
       y = "TMB",
       fill = "Mutation") +
  theme_minimal() +
  theme(axis.text.x = element_text(),
        text = element_text(family = "Calibri"))
