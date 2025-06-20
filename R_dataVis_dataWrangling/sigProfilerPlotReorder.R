### This code is designed to re-order samples to show paired samples next to each other
### It will replot your de novo mutational signature data

## Import libraries ------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)

## Data wrangling --------------------------------------------------------------
# Import data
signatures <- read.table("/path/to/SBS96/All_Solutions/SBS96_5_Signatures/Activities/SBS96_S5_NMF_Activities.txt",
                         header = TRUE, sep = "\t")

# Import corresponding name file
sampleNames <- read.table("/path/to/sampleNames.tsv",
                          header = TRUE, sep = "\t")

# Replace names to match sampleNames table
signatures <- signatures %>%
  left_join(sampleNames, by = c("Samples" = "sourceID")) %>%
  mutate(ID = ifelse(!is.na(shortID), shortID, Samples)) %>%
  select(-shortID) %>% 
  mutate(Samples = ID) %>%
  select(-ID) %>% 
  mutate(ID_number = sub("B|P", "", Samples)) %>%
  arrange(ID_number, Samples)

# Modify Samples to a factor with levels in the desired order
signatures <- signatures %>%
  mutate(Samples = factor(Samples, levels = unique(Samples))) # Set custom order of Samples

# Reshape to long format
longSignatures <- signatures %>%
  pivot_longer(
    cols = starts_with("SBS"), # All columns representing signature counts
    names_to = "Signature", # New column for mutation signature types
    values_to = "Count" # New column for counts
  )

## Plotting --------------------------------------------------------------------
# Define a custom color palette
custom_colors <- c(
  "SBS96A" = "# 619b8a",
  "SBS96B" = "# a1c181",
  "SBS96C" = "# fcca46",
  "SBS96D" = "# fe7f2d",
  "SBS96E" = "# 233d4d")

# Plot signatures in stacked bar chart with custom colors
ggplot(longSignatures, aes(x = Samples, y = Count, fill = Signature)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Mutational Signatures Across Samples",
    x = "Samples",
    y = "Mutation Count"
  ) +
  scale_fill_manual(values = custom_colors) +  # Apply the custom color palette
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


