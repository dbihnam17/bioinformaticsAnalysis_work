### This code is designed to compare mutational contexts indicative of
### APOBEC and AID activity between samples

### It uses the Samples.txt file generated from running sigProfiler,
### basically just a SBS96 matrix

## Import packages -------------------------------------------------------------
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

## Data wrangling --------------------------------------------------------------
# Import data
counts <- read_table("/path/to/SBS96/Samples.txt")

normCounts <- counts
normCounts[, -1] <- sweep(counts[, -1], 2, colSums(counts[, -1]), FUN = "/")

labCounts <- normCounts %>%
  mutate(Mutation = case_when(
    MutationType %in% c('T[C>T]A', 'T[C>T]C', 'T[C>T]G', 'T[C>T]T', 'T[C>G]A', 'T[C>G]C', 'T[C>G]T', 'T[C>A]A') ~ "APOBEC",
    MutationType %in% c('G[C>T]T, G[C>T]C', 'G[C>T]A', 'G[C>G]T', 'A[T>A]A', 'T[T>C]T', 'A[T>A]T', 'T[T>A]T') ~ "AID",
    TRUE ~ 'Other'))

labCountsSummarized <- labCounts %>%
  select(-MutationType) %>% # Drop the MutationType column
  group_by(Mutation) %>%
  summarize(across(everything(), sum, na.rm = TRUE))

labCountsLong <- labCountsSummarized %>%
  pivot_longer(
    cols = -Mutation, # All columns except Mutation
    names_to = "Sample", # Name of the new column for sample names
    values_to = "Proportion") # Name of the new column for proportion values

## Plotting --------------------------------------------------------------------
# Set a color palette
myColors <- c("APOBEC" = "# ffba08", "AID" = "# d00000", "Other" = "# 1c3144")

# Plot bar graph
ggplot(labCountsLong, aes(x = Sample, y = Proportion, fill = Mutation)) +
  geom_bar(stat = "identity") +
  labs(title = "AA: APOBEC and AID Activity", x = "Sample", y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1.0))+
  scale_fill_manual(values = myColors)





