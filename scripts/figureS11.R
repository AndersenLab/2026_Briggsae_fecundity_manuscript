# ============================================================
# Five-day Manhattan plots: C. briggsae fecundity GWAS
# ============================================================
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
# Load packages
library(data.table)
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)


# ------------------------------------------------------------
# 1. File locations
# ------------------------------------------------------------

base_dir <- "/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/nemascan"

day_files <- list(
  
  day1 = list(
    gwa = file.path(base_dir, "20261008_262strains_day1/avg_count_d1/gwa/inbred/avg_count_d1_inbred.gwa"),
    independent_tests = file.path(base_dir, "20261008_262strains_day1/independent_tests.txt")
  ),
  
  day2 = list(
    gwa = file.path(base_dir, "20261008_262strains_day2/avg_count_d2/gwa/inbred/avg_count_d2_inbred.gwa"),
    independent_tests = file.path(base_dir, "20261008_262strains_day2/independent_tests.txt")
  ),
  
  day3 = list(
    gwa = file.path(base_dir, "20261008_262strains_day3/avg_count_d3/gwa/inbred/avg_count_d3_inbred.gwa"),
    independent_tests = file.path(base_dir, "20261008_262strains_day3/independent_tests.txt"),
    qtl = file.path(base_dir, "20261008_262strains_day3/avg_count_d3/gwa/inbred/avg_count_d3_inbred_qtl.tsv"),
    qtl_chromosome = "X"
  ),
  
  day4 = list(
    gwa = file.path(base_dir, "20261008_262strains_day4/avg_count_d4/gwa/inbred/avg_count_d4_inbred.gwa"),
    independent_tests = file.path(base_dir, "20261008_262strains_day4/independent_tests.txt"),
    qtl = file.path(base_dir, "20261008_262strains_day4/avg_count_d4/gwa/inbred/avg_count_d4_inbred_qtl.tsv"),
    qtl_chromosome = "IV"
  ),
  
  day5 = list(
    gwa = file.path(base_dir, "20261008_262strains_day5/avg_count_d5/gwa/inbred/avg_count_d5_inbred.gwa"),
    independent_tests = file.path(base_dir, "20261008_262strains_day5/independent_tests.txt"),
    qtl = file.path(base_dir, "20261008_262strains_day5/avg_count_d5/gwa/inbred/avg_count_d5_inbred_qtl.tsv"),
    qtl_chromosome = "IV"
  )
)


# ------------------------------------------------------------
# 2. Significance colours and settings
# ------------------------------------------------------------

sig.colors <- c(
  "BF" = "red",
  "EIGEN" = "#EE4266",
  "NONSIG" = "black"
)

sig.alpha <- c(
  "BF" = 1,
  "EIGEN" = 1,
  "NONSIG" = 0.25
)

sig.linetypes <- c(
  "BF" = "solid",
  "EIGEN" = "dashed"
)


# ------------------------------------------------------------
# 3. Function to load and prepare GWAS data
# ------------------------------------------------------------

prepare_gwas <- function(gwa_file, independent_tests_file) {
  
  # Load GWAS results
  dat <- fread(gwa_file) %>%
    mutate(
      CHROM = as.character(CHR),
      POS = as.numeric(POS),
      marker = paste0(CHROM, ":", POS),
      log10p = -log10(P)
    ) %>%
    filter(CHROM != "MtDNA")
  
  
  # Calculate Bonferroni threshold
  n_markers <- n_distinct(dat$marker)
  
  BF <- -log10(0.05 / n_markers)
  
  
  # Load independent test count
  independent_tests <- read.table(
    independent_tests_file,
    header = FALSE
  )
  
  EIGEN <- -log10(
    0.05 / independent_tests[[1]][1]
  )
  
  
  # Assign significance categories
  dat <- dat %>%
    mutate(
      sig = case_when(
        log10p > BF ~ "BF",
        log10p > EIGEN ~ "EIGEN",
        TRUE ~ "NONSIG"
      ),
      sig = factor(
        sig,
        levels = c("NONSIG", "EIGEN", "BF")
      )
    )
  
  
  # Threshold dataframe
  sig.df <- tibble(
    name = c("BF", "EIGEN"),
    threshold = c(BF, EIGEN)
  )
  
  
  return(
    list(
      data = dat,
      thresholds = sig.df,
      BF = BF,
      EIGEN = EIGEN
    )
  )
}


# ------------------------------------------------------------
# 4. Function to create Manhattan plots
# ------------------------------------------------------------

make_manhattan_plot <- function(prepared_data,
                                qtl_file = NULL,
                                qtl_chromosome = NULL,
                                day_label = NULL) {
  
  dat <- prepared_data$data
  sig.df <- prepared_data$thresholds
  
  p <- ggplot(dat)
  
  
  # Add QTL interval, if available
  if (!is.null(qtl_file)) {
    
    qtl.data <- read_tsv(
      qtl_file,
      show_col_types = FALSE
    ) %>%
      mutate(
        CHROM = as.character(CHROM)
      ) %>%
      filter(CHROM == qtl_chromosome) %>%
      mutate(
        xmin = as.numeric(startPOS) / 1e6,
        xmax = as.numeric(endPOS) / 1e6
      )
    
    
    p <- p +
      geom_rect(
        data = qtl.data,
        aes(
          xmin = xmin,
          xmax = xmax,
          ymin = -Inf,
          ymax = Inf
        ),
        fill = "rosybrown1",
        alpha = 0.25,
        inherit.aes = FALSE
      )
  }
  
  
  # Manhattan plot
  p <- p +
    
    geom_point(
      aes(
        x = POS / 1e6,
        y = log10p,
        colour = sig,
        alpha = sig
      ),
      size = 0.25
    ) +
    
    scale_alpha_manual(
      values = sig.alpha,
      drop = FALSE
    ) +
    
    scale_colour_manual(
      values = sig.colors,
      drop = FALSE
    ) +
    
    geom_hline(
      data = sig.df,
      aes(
        yintercept = threshold,
        linetype = name
      )
    ) +
    
    scale_linetype_manual(
      values = sig.linetypes
    ) +
    
    facet_grid(
      ~CHROM,
      scales = "free_x"
    ) +
    
    labs(
      title = day_label,
      x = "Genomic position (Mb)",
      y = expression(-log[10](italic(p)))
    ) +
    
    theme_bw(
      base_size = 9,
      base_family = "Helvetica"
    ) +
    
    theme(
      text = element_text(
        family = "Helvetica",
        colour = "black"
      ),
      
      plot.title = element_text(
        size = 9,
        face = "bold",
        hjust = 0
      ),
      
      legend.position = "none",
      
      panel.grid = element_blank(),
      
      strip.text.x = element_text(
        size = 8
      ),
      
      axis.text.x = element_blank(),
      
      axis.text.y = element_text(
        size = 7
      ),
      
      axis.title.x = element_text(
        size = 8
      ),
      
      axis.title.y = element_text(
        size = 8
      ),
      
      axis.ticks.x = element_line(),
      
      plot.margin = margin(
        1, 2, 1, 2,
        unit = "pt"
      )
    )
  
  return(p)
}


# ------------------------------------------------------------
# 5. Prepare data for all five days
# ------------------------------------------------------------

prepared <- lapply(
  day_files,
  function(files) {
    
    prepare_gwas(
      gwa_file = files$gwa,
      independent_tests_file = files$independent_tests
    )
    
  }
)


# ------------------------------------------------------------
# 6. Create individual Manhattan plots
# ------------------------------------------------------------

# Day 1
manhattan_plot_1 <- make_manhattan_plot(
  prepared_data = prepared$day1,
  day_label = "Day 1"
)


# Day 2
manhattan_plot_2 <- make_manhattan_plot(
  prepared_data = prepared$day2,
  day_label = "Day 2"
)


# Day 3: highlight chromosome X QTL
manhattan_plot_3 <- make_manhattan_plot(
  prepared_data = prepared$day3,
  #qtl_file = day_files$day3$qtl,
  #qtl_chromosome = day_files$day3$qtl_chromosome,
  day_label = "Day 3"
)


# Day 4: highlight chromosome IV QTL
manhattan_plot_4 <- make_manhattan_plot(
  prepared_data = prepared$day4,
 # qtl_file = day_files$day4$qtl,
 # qtl_chromosome = day_files$day4$qtl_chromosome,
  day_label = "Day 4"
)


# Day 5: highlight chromosome IV QTL
manhattan_plot_5 <- make_manhattan_plot(
  prepared_data = prepared$day5,
 # qtl_file = day_files$day5$qtl,
 # qtl_chromosome = day_files$day5$qtl_chromosome,
  day_label = "Day 5"
)


# View individual plots
manhattan_plot_1
manhattan_plot_2
manhattan_plot_3
manhattan_plot_4
manhattan_plot_5


# ------------------------------------------------------------
# 7. Make five vertically stacked strips
# ------------------------------------------------------------

make_strip <- function(p, show_x_title = FALSE) {
  
  p +
    
    labs(
      x = if (show_x_title) {
        "Genomic position (Mb)"
      } else {
        NULL
      },
      title = NULL
    ) +
    
    theme(
      
      axis.text.x = element_blank(),
      
      axis.title.x = element_text(
        size = 12,
        family = "Helvetica"
      ),
      
      axis.title.y = element_text(
        size = 12,
        family = "Helvetica"
      ),
      
      axis.text.y = element_text(
        size = 11
      ),
      
      strip.text.x = element_text(
        size = 11
      ),
      
      plot.margin = margin(
        1, 2, 1, 2,
        unit = "pt"
      )
    )
}


# Combine all five plots
combined_strips <- (
  
  make_strip(manhattan_plot_1) /
    
    make_strip(manhattan_plot_2) /
    
    make_strip(manhattan_plot_3) /
    
    make_strip(manhattan_plot_4) /
    
    make_strip(
      manhattan_plot_5,
      show_x_title = TRUE
    )
  
) +
  
  plot_layout(
    ncol = 1,
    heights = rep(1, 5)
  ) +
  
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "",
    tag_suffix = ""
  ) &
  
  theme(
    plot.tag = element_text(
      family = "Helvetica",
      face = "bold",
      size = 12,
      colour = "black"
    ),
    
    plot.tag.position = c(0, 1)
  )


# Display combined figure
combined_strips

ggsave(
  filename = "/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/plots/S11.tif",
  plot = combined_strips,
  width = 7.5,
  height = 9,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

