### This script is intended to plot the variance of each decomposed mutational signature
### across all samples from each group in a violin plot, then use a Wilcoxon test
### to determine if there is a significant difference in the prevalence of signatures
### between the groups

## Load in libraries -----------------------------------------------------------

library(tidyverse)
library(ggpubr)  # For adding p-values to ggplot

## Data wrangling --------------------------------------------------------------

# (Group 1)
g1 <- read.table("/path/to/group1/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt",
                 header = TRUE, row.names = 1)

# (Group 2)
g2 <- read.table("/path/to/group2/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt",
                   header = TRUE, row.names = 1)

# Add group labels
g1 <- g1 %>% 
  rownames_to_column("Sample") %>% 
  mutate(Group = "group1")
g2 <- g2 %>% 
  rownames_to_column("Sample") %>% 
  mutate(Group = "group2")

# Combine both groups
combined <- bind_rows(g1, g2)

# Convert data to long format
longData <- combined %>% 
  pivot_longer(cols = -c(Sample, Group), names_to = "Signature", values_to = "Count")

# Filter out signatures not present in either group
longDataFilt <- longData %>%
  group_by(Signature) %>%
  filter(sum(Count) > 0) %>%
  ungroup()

## Statistical testing ---------------------------------------------------------

# Perform Wilcoxon tests for each signature
p_values <- longDataFilt %>%
  group_by(Signature) %>%
  summarise(p_value = wilcox.test(Count ~ Group)$p.value)

# Join the p-values with the longDataFilt
longDataFilt <- left_join(longDataFilt, p_values, by = "Signature")

## Plotting --------------------------------------------------------------------

# Plot as violin plot with raw counts
ggplot(longDataFilt, aes(x = Group, y = Count, fill = Group)) +
  geom_violin(trim = FALSE, scale = "area", alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 0.8) +
  facet_wrap(~ Signature, scales = "free_y") +
  theme_minimal(base_size = 14) +
  theme(text = element_text(family = "Calibri")) +
  labs(title = "COSMIC Signature Mutational Contribution by Group") +
  geom_text(data = p_values, aes(x = 1.5, y = Inf, label = paste0("p = ", format(p_value, digits = 2))),
            inherit.aes = FALSE, size = 4, color = "black", hjust = 0.5, vjust = 1)
