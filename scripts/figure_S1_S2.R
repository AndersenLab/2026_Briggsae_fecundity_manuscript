
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
library(lme4)
library(sommer)
library(MASS)
library(dplyr)
#install.packages("readr")
library(readr)
library(tidyverse)


df3_1 <- read.csv (file = "/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/Trial1_Precision_recall.csv",header = TRUE,sep = ',')

df4_1 <- read.csv (file = "/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/Trial2_Precision_recall.csv",header = TRUE,sep = ',')


####Figure S1 #####
##               ##
##               ##
##               ##
####Figure S1 #####


df_long_3 <- df3_1 %>%
  mutate(
    Day = as.numeric(
      sub(".*_D([0-9]+)_.*", "\\1", Assay_Day_Replicate)
    ),
    
    across(
      c(
        Precision_label,
        Recall_label,
        Precision_lawn,
        Recall_lawn,
        Precision_rest,
        Recall_rest,
        Precision_total,
        Recall_total
      ),
      ~ as.numeric(
        na_if(as.character(.x), "#DIV/0!")
      )
    )
  ) %>%
  pivot_longer(
    cols = c(
      Precision_label,
      Recall_label,
      Precision_lawn,
      Recall_lawn,
      Precision_rest,
      Recall_rest,
      Precision_total,
      Recall_total
    ),
    names_to = c("Metric", "Region"),
    names_sep = "_",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))


df_long_3$Region <- factor(
  df_long_3$Region,
  levels = c("label", "lawn", "rest", "total")
)


df_summary <- df_long_3 %>%
  group_by(Day, Region, Metric) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    .groups = "drop"
  )


p3 <- ggplot(
  df_summary,
  aes(
    x = Day,
    y = Mean,
    color = Metric,
    group = Metric,
    linetype = Metric
  )
) +
  geom_errorbar(
    aes(
      ymin = pmax(0, Mean - SD),
      ymax = pmin(1, Mean + SD)
    ),
    width = 0.12,
    linewidth = 0.5
  ) +
  geom_line(
    linewidth = 0.5
  ) +
  geom_point(
    size = 1.0
  ) +
  facet_wrap(~Region) +
  scale_color_manual(
    values = c(
      Precision = "black",
      Recall = "black"
    )
  ) +
  scale_linetype_manual(
    values = c(
      Precision = "solid",
      Recall = "dashed"
    )
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.4),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0))
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(
      family = "Helvetica",
      color = "black"
    ),
    axis.text = element_text(
      size = 11,
      family = "Helvetica",
      color = "black"
    ),
    axis.title = element_text(
      size = 12,
      family = "Helvetica",
      color = "black"
    ),
    strip.text = element_text(
      size = 12,
      family = "Helvetica",
      color = "black"
    )
  ) +
  labs(
    y = "Precision/Recall",
    x = "Day"
  )

p3

ggsave("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/plots/figureS1.tif",
       plot = p3, width = 7.5, height = 2.5, dpi = 350, units = "in", bg = "white")

####Figure S2 #####
##               ##
##               ##
##               ##
####Figure S2 #####

df_long_4 <- df4_1 %>%
  mutate(
    Day = as.numeric(
      sub(".*_D([0-9]+)_.*", "\\1", Assay_Day_Replicate)
    ),
    
    across(
      c(
        Trial_4_Presicion_0.45,
        Trial_4_Presicion_0.55,
        Trial_4_Presicion_0.65,
        Trial_4_Presicion_0.75,
        Trial_4_Recall_0.45,
        Trial_4_Recall_0.55,
        Trial_4_Recall_0.65,
        Trial_4_Recall_0.75
      ),
      ~ as.numeric(
        na_if(as.character(.x), "#DIV/0!")
      )
    )
  ) %>%
  pivot_longer(
    cols = c(
      Trial_4_Presicion_0.45,
      Trial_4_Recall_0.45,
      Trial_4_Presicion_0.55,
      Trial_4_Recall_0.55,
      Trial_4_Presicion_0.65,
      Trial_4_Recall_0.65,
      Trial_4_Presicion_0.75,
      Trial_4_Recall_0.75
    ),
    names_to = c("Metric", "Region"),
    names_pattern = "Trial_4_(Presicion|Recall)_(.*)",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = recode(
      Metric,
      Presicion = "Precision"
    )
  ) %>%
  filter(!is.na(Value))

df_summary_4 <- df_long_4 %>%
  group_by(Day, Region, Metric) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    .groups = "drop"
  )


p4 <- ggplot(
  df_summary_4,
  aes(
    x = Day,
    y = Mean,
    color = Metric,
    group = Metric,
    linetype = Metric
  )
) +
  geom_errorbar(
    aes(
      ymin = pmax(0, Mean - SD),
      ymax = pmin(1, Mean + SD)
    ),
    width = 0.12,
    linewidth = 0.5
  ) +
  geom_line(
    linewidth = 0.5
  ) +
  geom_point(
    size = 1.0
  ) +
  facet_wrap(~Region) +
  scale_color_manual(
    values = c(
      Precision = "black",
      Recall = "black"
    )
  ) +
  scale_linetype_manual(
    values = c(
      Precision = "solid",
      Recall = "dashed"
    )
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.4),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0))
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(
      family = "Helvetica",
      color = "black"
    ),
    axis.text = element_text(
      size = 11,
      family = "Helvetica",
      color = "black"
    ),
    axis.title = element_text(
      size = 12,
      family = "Helvetica",
      color = "black"
    ),
    strip.text = element_text(
      size = 12,
      family = "Helvetica",
      color = "black"
    )
  ) +
  labs(
    y = "Precision/Recall",
    x = "Day"
  )

p4


ggsave("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/plots/figureS2.tif",
       plot = p4, width = 7.5, height = 2.5, dpi = 350, units = "in", bg = "white")
