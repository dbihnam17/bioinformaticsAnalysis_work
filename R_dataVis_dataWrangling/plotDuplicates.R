### This script was intended to quickly visualize duplicate reads detected in
### flagstat reports between samples from different WGS library preps

## Load libraries --------------------------------------------------------------
library(ggplot2)
library(dplyr)

## Data wrangling --------------------------------------------------------------
# Create the data frame
data <- data.frame(
  Sample = rep(c("sample1", "sample2", "sample3", "sample3_reprocessed"), each = 2),
  Category = rep(c("Total Reads", "Duplicate Reads"), times = 4),
  Reads = c(1267328140, 181762913, 1195282284, 380066631, 685793732, 226166876, 1263253123, 370004053)
)

# Convert Sample to a factor to preserve the original order
data$Sample <- factor(data$Sample, levels = unique(data$Sample))

# Calculate the proportions
data <- data %>%
  group_by(Sample) %>%
  mutate(Proportion = Reads / sum(Reads))

## Plot the normalized stacked bar chart ---------------------------------------
ggplot(data, aes(x = Sample, y = Proportion, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Normalized Total and Duplicate Reads by Sample",
    x = "Sample",
    y = "Proportion",
    fill = "Category"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Total Reads" = "steelblue", "Duplicate Reads" = "firebrick")) +
  scale_y_continuous(labels = scales::percent_format())
