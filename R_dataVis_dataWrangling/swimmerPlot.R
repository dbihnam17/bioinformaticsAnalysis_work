### This script is used to create swimmer plots with patient treatment data
### Initially written for Castro/Fonseca lab MM CAR project

## Publication in review: Pre-infusion ferritin a predictive biomarker for CAR T-cell therapy in R/R multiple myeloma
## Authors: Eider Moreno-Cortes, Mariano Arribas, David Martinez, Juana C Figueroa Aguirre, Rohan Patel, Daniel T Bihnam, Udit Yadav, Erin Wiedmeier-Nutor, Leif Bergsagel, Saurabh Chhabra, Rafael Fonseca, Januario E. Castro and M.D.
## Frontiers in Hematology, section Gene Therapy, Cell Therapy and Hematology
## Submitted on: 25 Apr 2025

## Initial imports -------------------------------------------------------------
# Load necessary libraries
library(ggplot2)
# install.packages("devtools")
# devtools::install_github("biostatsPMH/swimplot", ref="main")
library(swimplot)
library(dplyr)
library(tidyr)
library(readxl)
library(openxlsx)
library(stringr)

## Data imports
# Import STATUS data
statusData <- read_excel('/path/to/patientSTATUS.xlsx')

# Import EVENTS data
eventData <- read_excel('/path/to/patientEVENTS.xlsx')

## Process STATUS data ---------------------------------------------------------

# Make data long
statusLong <- statusData %>%
  pivot_longer(
    cols = starts_with('response'),
    names_to = 'month',
    names_prefix = 'response',
    values_to = 'response'
  )

# Set responseEndMonth columns
# Mapping function
map_response_end_month <- function(month, totMonths) {
  if (month == 1) {
    end_month <- ifelse(totMonths >= 3, 3, totMonths)
  } else if (month == 3) {
    end_month <- ifelse(totMonths >= 6, 6, totMonths)
  } else if (month == 6) {
    end_month <- ifelse(totMonths >= 12, 12, totMonths)
  } else if (month == 12) {
    end_month <- ifelse(totMonths >= 24, 24, totMonths)
  } else if (month == 24) {
    end_month <- 'R'
  } else {
    end_month <- totMonths
  }
  return(end_month)
}

# Apply the mapping function to each row
statusLong$responseEndMonth <- mapply(map_response_end_month, statusLong$month, statusLong$totMonths)

# Filtering out unneccessary columns
statusLong <- statusLong %>%
  filter(!(relapsed == 0 & responseEndMonth == 'R'))

# Drop NA values
statusLong <- statusLong %>%
  drop_na(response)

# Fill 'R' with relpase month
for (i in 1:nrow(statusLong)) {
  if (statusLong$responseEndMonth[i] == 'R') {
    mrn_value <- statusLong$mrn[i]
    relapse_month <- eventData$relapseMonth[eventData$mrn == mrn_value]
    statusLong$responseEndMonth[i] <- relapse_month
  }
}

for (i in 1:nrow(statusLong)) {
  if (statusLong$month[i] == 'R') {
    mrn_value <- statusLong$mrn[i]
    relapse_month <- eventData$relapseMonth[eventData$mrn == mrn_value]
    statusLong$month[i] <- relapse_month
  }
}

# Convert columns to numerical values
statusLong$responseEndMonth <- as.numeric(statusLong$responseEndMonth)
statusLong$month <- as.numeric(statusLong$month)

# Export for manual QC/editing in Excel
# write.xlsx(statusLong, '/path/to/statusLong_modified.xlsx')

# Re-import manually edited file
# statusDataFinal <- read_excel('/path/to/statusLong_modifiedCUT.xlsx')

statusDataFinaldf <- as.data.frame(statusDataFinal)

## Process EVENT data ----------------------------------------------------------

eventDataFinal <- eventData %>%
  mutate(status = case_when(
    status == 1 ~ 'Deceased (related)',
    status == 2 ~ 'Deceased (unrelated)'
  ))

eventDataFinal <- eventDataFinal %>%
  mutate(relapsed = case_when(
    relapsed == 'Y' ~ 'Relapse',
    relapsed == 'N' ~ NA_character_
  ))

eventDataFinal <- eventDataFinal %>%
  mutate(alive = case_when(
    alive == 0 ~ 'Alive'
  ))

eventDatadf <- as.data.frame(eventDataFinal)

eventDatadfArrow <- eventDatadf %>% filter(finStatus != 'PD')

eventDatadfDecR <- eventDatadf %>% filter(status == 'Deceased (related)')
eventDatadfDecU <- eventDatadf %>% filter(status == 'Deceased (unrelated)')


## Plot status data ------------------------------------------------------------

swimPlot <- swimmer_plot(df = statusDataFinaldf,id='mrn',end='responseEndMonth',name_fill='response',
                         col='white',alpha=0.85,width=0.75)
swimPlot

swimPlotPoint <- swimPlot + swimmer_points(df_points=
                                             eventDatadfDecR,id='mrn',time='totMonths', name_shape =
                                             'status',size=2.15,fill='black',col='black')+swimmer_points(df_points=
                                                                                                             eventDatadfDecU,id='mrn',time='totMonths', name_shape =
                                                                                                             'status',size=2.15,fill='# A9A9A9',col='black')
swimPlotPoint

swimPlotPoint2 <- swimPlotPoint + swimmer_arrows(df_arrows = eventDatadfArrow,id='mrn',arrow_start='totMonths',
                                                 cont='alive',name_col='finStatus',show.legend=FALSE,
                                                 type = 'open',cex=0.5,length=0.06,arrow_positions = c(0.1,2))+scale_color_discrete(drop=FALSE)
swimPlotPoint2

swimPlotPoint3 <- swimPlotPoint2 + swimmer_points(df_points=
                                                    eventDatadf,id='mrn',time='relapseMonth', name_shape =
                                                    'relapsed',size=1.85,fill='red',col='black')

swimPlotPoint3

swimPlotFin <- swimPlotPoint3 +
  scale_y_continuous(minor_breaks = seq(0, 60, by = 4), breaks = seq(0, 60, by = 6)) +
  theme(panel.grid.major.x = element_line(colour = "grey95", size = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_manual(name = "Treatment Response", values = c("sCR" = "# 1f4d2e", "CR" = "# 008040", "VGPR" = '# 6bb359', "PR" = "# ffcc00", "SD" = "# ff8c19", "PD" = "# ff400d"), 
                    breaks = c("sCR", "CR", "VGPR", "PR", "SD", "PD")) +
  scale_color_manual(name = "Treatment Response", values = c("sCR" = "# 1f4d2e", "CR" = "# 008040", "VGPR" = '# 6bb359', "PR" = "# ffcc00", "SD" = "# ff8c19", "PD" = "# ff400d"), 
                     breaks = c("sCR", "CR", "VGPR", "PR", "SD", "PD")) +
  scale_shape_manual(name = "Patient Status", values = c("Deceased (related)" = 22, "Deceased (unrelated)" = 25, "Relapse" = 23), 
                     breaks = c("Deceased (related)", "Deceased (unrelated)", "Relapse"))+
  ylab('Months Since CAR-T Infusion')+
  xlab('')

swimPlotFin

swimPlotFin <- swimPlotPoint3 +
  scale_y_continuous(minor_breaks = seq(0, 60, by = 4), breaks = seq(0, 60, by = 6)) +
  theme(panel.grid.major.x = element_line(colour = "grey95", size = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_manual(name = "Treatment Response", 
                    values = c("sCR" = "# 1f4d2e", "CR" = "# 008040", "VGPR" = '# 6bb359', 
                               "PR" = "# ffcc00", "SD" = "# ff8c19", "PD" = "# ff400d"), 
                    breaks = c("sCR", "CR", "VGPR", "PR", "SD", "PD")) +
  scale_color_manual(name = "Treatment Response", 
                     values = c("sCR" = "# 1f4d2e", "CR" = "# 008040", "VGPR" = '# 6bb359', 
                                "PR" = "# ffcc00", "SD" = "# ff8c19", "PD" = "# ff400d"), 
                     breaks = c("sCR", "CR", "VGPR", "PR", "SD", "PD")) +
  scale_shape_manual(name = "Patient Status", 
                     values = c("Deceased (related)" = 22, "Deceased (unrelated)" = 25, 
                                "Relapse" = 23), 
                     breaks = c("Deceased (related)", "Deceased (unrelated)", "Relapse")) +
  ylab('Months Since CAR-T Infusion') +
  xlab('') +
  theme(legend.title = element_text(size = 12, family = "Calibri"),
        legend.text = element_text(size = 10, family = "Calibri"),
        legend.key.size = unit(1.5, "lines"),
        legend.spacing.y = unit(0.5, "lines"),
        legend.key.height = unit(1, "lines"),
        legend.key.width = unit(1.5, "lines"),
        text = element_text(size = 10, family = "Calibri"))

swimPlotFin







