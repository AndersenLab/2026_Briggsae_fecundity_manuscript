
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

# ============================================================
# 1. LIFETIME FECUNDITY
# ============================================================

df_summary <- data_repwise %>%
  group_by(Strain) %>%
  summarise(
    Div = first(Div),
    mean_fec = mean(total_count, na.rm = TRUE),
    sem_fec = sd(total_count, na.rm = TRUE) /
      sqrt(sum(!is.na(total_count))),
    .groups = "drop"
  )

# Get the strain order used in the lifetime fecundity plot
strain_order <- df_summary %>%
  arrange(mean_fec) %>%
  pull(Strain)


# Apply the strain order
df_summary <- df_summary %>%
  mutate(
    Strain = factor(
      Strain,
      levels = strain_order
    )
  )


# ---------------------------------------------------------
# 1. Convert daily fecundity to cumulative fecundity
#    for each replicate
# ---------------------------------------------------------
data_daywise <- data_repwise %>%
  pivot_longer(
    cols = c(d1, d2, d3, d4, d5),
    names_to = "Day",
    values_to = "fecundity"
  ) %>%
  mutate(
    Day = factor(
      Day,
      levels = c("d1", "d2", "d3", "d4", "d5"),
      labels = c(
        "1",
        "2",
        "3",
        "4",
        "5"
      )
    ),
    # Use EXACTLY the same strain order as the lifetime plot
    Strain = factor(
      Strain,
      levels = strain_order
    )
  )


# ---------------------------------------------------------
# 1. Average daily fecundity across assays and replicates
#    for each strain, then calculate cumulative fecundity
# ---------------------------------------------------------

data_cumulative <- data_daywise %>%
  group_by(Strain, Day) %>%
  summarise(
    Div = first(Div),
    mean_fecundity = mean(fecundity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Strain, Day) %>%
  group_by(Strain) %>%
  mutate(
    cumulative_fec = cumsum(mean_fecundity)
  ) %>%
  ungroup()


# ---------------------------------------------------------
# 2. Summary: one value per strain per day
# ---------------------------------------------------------

df_summary_cumulative <- data_cumulative %>%
  dplyr::select(Strain, Div, Day, cumulative_fec)


# ---------------------------------------------------------
# 3. Colors
# ---------------------------------------------------------

divergent_col <- "#0072B2"
nondivergent_col <- "#E69F00"


# ---------------------------------------------------------
# 4. Plot cumulative fecundity
# ---------------------------------------------------------


p_cumulative <- ggplot(
  df_summary_cumulative,
  aes(
    x = Day,
    y = cumulative_fec,
    group = Strain,
    color = Div
  )
) +
  geom_line(
    linewidth = 0.25,
    alpha = 0.5
  ) +
  geom_point(
    size = 0.6,
    alpha = 0.5
  ) +
  scale_color_manual(
    values = c(
      "Div" = divergent_col,
      "NoDiv" = nondivergent_col
    )
  ) +
  labs(
    x = "Day",
    y = "Cumulative fecundity",
    color = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(
      family = "Helvetica",
      size = 12,
      color = "black"
    ),
    axis.text = element_text(
      family = "Helvetica",
      size = 12,
      color = "black"
    ),
    axis.title = element_text(
      family = "Helvetica",
      size = 12,
      color = "black"
    ),
    legend.text = element_text(
      family = "Helvetica",
      size = 12,
      color = "black"
    ),
    legend.position = "none"
  )
p_cumulative


ggsave("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/plots/FigureS7.tif",
       plot = p_cumulative, width = 7.5, height = 2.5, units = "in", dpi = 350, bg = "white")

