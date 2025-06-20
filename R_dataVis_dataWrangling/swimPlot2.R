### This script is meant to plot a swimmer's plot based on the <project>
### patient data and available timepoints

## Due to differences in data format, swimmerPlot.R was not used for this
## Also, the PI that requested this plot had a different vision than the one I showed
## that resulted from swimmerPlot.R

## Load libraries --------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyverse)

## Data wrangling --------------------------------------------------------------
# Load in data
data <- read.csv('/path/to/example/swimPlot/swimPlotDataPrep.csv', header = TRUE)

# Make sure ID is treated as a factor
data$ID <- factor(data$ID)

# Reshape to long format
data_long <- data %>%
  pivot_longer(
    cols = c(baseline, on_treatment, progression, post_progression),
    names_to = "timepoint",
    values_to = "months"
  ) %>%
  filter(!is.na(months))

# Count the number of timepoints each patient has
timepoint_count <- data_long %>%
  group_by(ID) %>%
  summarise(n_timepoints = n_distinct(timepoint)) %>%
  arrange(desc(n_timepoints))

# Join the timepoint count back to data_long
data_long <- data_long %>%
  left_join(timepoint_count, by = "ID")

# Create a new factor for grouping patients by their exact set of timepoints
data_long <- data_long %>%
  group_by(ID) %>%
  mutate(timepoint_group = paste(sort(unique(timepoint)), collapse = "_")) %>%
  ungroup()

# Create a unique identifier based on n_timepoints and timepoint_group
data_long <- data_long %>%
  mutate(unique_order = paste(n_timepoints, timepoint_group, sep = "_"))

# Generate a new ordered factor for ID based on the unique_order
timepoint_order <- data_long %>%
  distinct(ID, unique_order, n_timepoints) %>%
  arrange(desc(n_timepoints), unique_order) %>%
  pull(ID)

# Reorder the ID factor based on the new timepoint_order
data_long$ID <- factor(data_long$ID, levels = timepoint_order)

# Reorder the 'timepoint' factor levels
data_long$timepoint <- factor(data_long$timepoint, levels = c("baseline", "on_treatment", "progression", "post_progression"))

# Manually adjust the months (years) value for specific rows to avoid visual overlap
data_long <- data_long %>%
  mutate(months = case_when(
    ID == 1234 & timepoint == "on_treatment" ~ months - .2,
    ID == 5678 & timepoint == "progression" ~ months - .2,
    TRUE ~ months
  ))

## Plot ------------------------------------------------------------------------
ggplot(data_long, aes(x = months, y = fct_rev(factor(ID)), group = ID)) +
  geom_vline(xintercept = seq(-10, 10, by = 5), linetype = "solid", color = "gray80", size = 1) +
  geom_vline(xintercept = seq(0, 0, by = 1), linetype = "solid", color = "grey30", size = 1) + # Darker center line at T=0
  geom_line(size = 4, color = "grey50", show.legend = FALSE) +
  geom_point(size = 3, aes(fill = timepoint), shape = 21, color = "black", stroke = 0.75) +
  scale_fill_manual(values = c(
    "baseline" = "dodgerblue", 
    "on_treatment" = "# 32CD32", 
    "progression" = "# FFA500", 
    "post_progression" = "# d00000"
  ), labels = c("Baseline", "On Treatment", "Progression", "Relapse")) +
  labs(
    x = "Time in Years (Treatment Start T = 0)",
    y = "Patient ID",
    fill = "Timepoint",
    title = "Example Project CLL Patient Treatment Timeline"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    text = element_text(size = 10, family = "Calibri"),
    legend.text = element_text(size = 10)
  ) +
  guides(fill = guide_legend(title = "Timepoint"))


