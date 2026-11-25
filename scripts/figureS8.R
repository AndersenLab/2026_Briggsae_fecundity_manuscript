
rm(list = ls())
# Load libraries
library(tidyverse)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(ggrepel)
#install.packages("ggbeeswarm")
library(ggbeeswarm)


# Load GWAS data
d.all1 <- fread(
  "/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/nemascan/20261115_262strains_repday_1_2//dayof_D50_bin/gwa/inbred/dayof_D50_bin_inbred.gwa" )


# Format data
d.all1 <- d.all1 %>%
  mutate(
    CHROM = CHR,
    POS = POS,
    marker = paste0(CHR, ":", POS),
    log10p = -log10(P)
  )

# Remove mitochondria if present
d.all1 <- d.all1 %>%
  filter(CHROM != "MtDNA")

# Compute Bonferroni threshold
BF <- -log10(
  0.05 / nrow(d.all1 %>% distinct(marker))
)

# Load independent test count from NemaScan
total_independent_tests_1 <- read.table(
  "/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/nemascan/20261115_262strains_repday_1_2/independent_tests.txt"
)

EIGEN <- -log10(
  0.05 / total_independent_tests_1[[1]]
)

# Assign significance
d.all12 <- d.all1 %>%
  mutate(
    sig = case_when(
      log10p > BF ~ "BF",
      log10p > EIGEN & log10p <= BF ~ "EIGEN",
      TRUE ~ "NONSIG"
    )
  )

# Significance colors
sig.colors <- c(
  "BF" = "red",
  "EIGEN" = "#EE4266",
  "NONSIG" = "black"
)

# Threshold dataframe
sig.df <- tibble(
  name = c("BF", "EIGEN"),
  sig = c(BF, EIGEN)
)

# Manhattan plot
manhattan_plot_1 <- ggplot(d.all12) +
  geom_point(
    aes(
      x = POS / 1e6,
      y = log10p,
      colour = sig,
      alpha = sig
    ),
    size = 0.7
  ) +
  scale_alpha_manual(
    values = c(
      "BF" = 1,
      "EIGEN" = 1,
      "NONSIG" = 0.25
    )
  ) +
  scale_colour_manual(
    values = sig.colors
  ) +
  geom_hline(
    data = sig.df,
    aes(
      yintercept = sig,
      linetype = name
    )
  ) +
  scale_linetype_manual(
    values = c(
      "BF" = 1,
      "EIGEN" = 3
    )
  ) +
  facet_grid(
    ~CHROM,
    scales = "free_x"
  ) +
  labs(
    x = "Genomic position (Mb)",
    y = expression(-log[10](italic(p)))
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    strip.text.x = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_line()
  )

manhattan_plot_1


ggsave("/vast/eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/plots/figureS8.tif",
       plot = manhattan_plot_1,  width = 7.5, height = 2.5, dpi = 350, units = "in", bg = "white")
