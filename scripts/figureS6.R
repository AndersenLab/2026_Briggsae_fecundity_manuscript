

rm(list = ls())

library(ggplot2)
library(tidyr)
library(ggplot2)
library(dplyr)
library(ggrepel)
#install.packages("cowplot")
library(cowplot)
library(RColorBrewer)
#install.packages("readxl")
library(readxl)
library(dplyr)
#install.packages("readr")
library(readr)

data <- read.delim("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/311strains_avgBS_A3567_4_1_.tsv")
data_repwise <- read.delim("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/311strains_avgBS_A3567_4_1_repwise.tsv")

###Having only two colors:
### build color map with two colors:
# Define colors
divergent_col <- "#0072B2"
nondivergent_col <- "#E69F00"


fecundity_summary <- data_repwise %>%
  group_by(Strain, Div) %>%
  summarise(
    across(
      d1:d5,
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )


fecundity_long <- fecundity_summary %>%
  pivot_longer(
    cols = d1:d5,
    names_to = "Day",
    values_to = "Mean_fecundity"
  ) %>%
  mutate(
    Day = factor(
      Day,
      levels = paste0("d", 1:5),
      labels = as.character(1:5)
    )
  )

##Mean fecundity for each day##
library(dplyr)
library(tidyr)

df_daily_summary <- fecundity_summary %>%
  pivot_longer(
    cols = d1:d5,
    names_to = "Day",
    values_to = "Fecundity"
  ) %>%
  group_by(Day) %>%
  summarise(
    mean_fec = mean(Fecundity, na.rm = TRUE),
    sd_fec = sd(Fecundity, na.rm = TRUE),
    n = sum(!is.na(Fecundity)),
    .groups = "drop"
  )

df_daily_summary


p_fecundity <- ggplot(
  fecundity_long,
  aes(
    x = Day,
    y = Mean_fecundity
  )
) +
  
  # Boxplots
  geom_boxplot(
    fill = "grey90",
    color = "black",
    width = 0.6,
    outlier.shape = NA,
    linewidth = 0.4
  ) +
  
  # Individual strain means
  geom_jitter(
    aes(color = Div),
    width = 0.12,
    size = 0.5,
    alpha = 0.6
  ) +
  
  # Divergent / nondivergent colors
  scale_color_manual(
    values = c(
      "Div" = divergent_col,
      "NoDiv" = nondivergent_col
    )
  ) +
  
  labs(
    x = "Day",
    y = "Mean fecundity",
    color = "Divergence"
  ) +
  
  theme_classic(
    base_family = "Helvetica"
  ) +
  
  theme(
    
    axis.title = element_text(
      family = "Helvetica",
      size = 12,
      color = "black"
    ),
    
    axis.text = element_text(
      family = "Helvetica",
      size = 12,
      color = "black"
    ),
   
     legend.position = "none",
    
   # legend.title = element_text(
    #  family = "Helvetica",
    #  size = 10,
    #  color = "black"
   # ),
    
   # legend.text = element_text(
   #   family = "Helvetica",
   #   size = 10,
    #  color = "black"
   # ),
    
    strip.text = element_text(
      family = "Helvetica",
      size = 12,
      color = "black"
    ),
    
    strip.background = element_blank()
  )

p_fecundity

ggsave("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/plots/FigureS6.tif",
       plot = p_fecundity, width = 7.5, height = 2.5, units = "in", bg = "white")



