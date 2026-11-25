##────────────────────────────────────────────────
##  Brood size Manhattan plot Script - figure 3
##────────────────────────────────────────────────
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

#────────────────────────────────────────────────
# Panel A
#────────────────────────────────────────────────

# Load data
gwas_raw <- read_tsv(
  "/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/nemascan/20260903_1//avg_BS/gwa/inbred/avg_BS_inbred.gwa",
  col_types = cols()
)


gwas_raw <- gwas_raw %>% 
            mutate(
            CHROM = CHR,
            POS = POS,
            marker = paste0(CHR, ":", POS),
            log10p = -log10(P)
)



qtl_peaks <- read_tsv(
  "/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/nemascan/20260903_1/avg_BS/gwa/inbred/avg_BS_inbred_qtl.tsv"
)


qtl_chrX <- qtl_peaks %>%
  filter(CHROM == "X")

region_start <- qtl_chrX$startPOS / 1e6
region_end   <- qtl_chrX$endPOS  / 1e6


# Load brood size mapping
d.all <- fread("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/nemascan/20260903_1/avg_BS/gwa/inbred/avg_BS_inbred.gwa") 
d.all <- d.all %>% 
  mutate(
    CHROM = CHR,
    POS = POS,
    marker = paste0(CHR, ":", POS),
    log10p = -log10(P)
  )

# Remove mitochondria if present
d.all <- d.all %>%
  filter(CHROM != "MtDNA")

# Add log10p if not already present
if(!"log10p" %in% colnames(d.all)){
  d.all <- d.all %>%
    mutate(log10p = -log10(P))
}

# Compute thresholds
#────────────────────────────────────────
# Bonferroni (based on unique markers)
BF <- -log10(0.05 / nrow(d.all %>% distinct(marker)))

# Load independent test count from NemaScan
total_independent_tests <- read.table(
  "/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/nemascan//20260903_1/independent_tests.txt"
)
EIGEN <- -log10(0.05 / total_independent_tests[[1]])

# Assign significance
d.all2 <- d.all %>%
  mutate(sig = case_when(
    log10p > BF ~ "BF",
    log10p > EIGEN & log10p <= BF ~ "EIGEN",
    TRUE ~ "NONSIG"
  ))

sig.colors <- c("BF" = "red",
                "EIGEN" = "#EE4266",
                "NONSIG" = "black")

sig.df <- tibble(
  name = c("BF", "EIGEN"),
  sig  = c(BF, EIGEN)
)

# Manhattan plot
manhattan_plot <- ggplot(d.all2) +
  theme_bw() +
  # geom_rect(
  #   data = tibble(
  #     CHROM = "V",
  #     xmin = region_start,
  #     xmax = region_end,
  #     ymin = -Inf,
  #     ymax = Inf
  #   ),
  #   aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
  #   fill = "rosybrown1",
  #   alpha = 0.25,
  #   inherit.aes = FALSE
  # ) +
  geom_point(aes(
    x = POS / 1e6,
    y = log10p,
    colour = sig,
    alpha = sig
  ), size = 0.7) +
  scale_alpha_manual(values = c("BF" = 1,
                                "EIGEN" = 1,
                                "NONSIG" = 0.25)) +
  scale_colour_manual(values = sig.colors) +
  geom_hline(data = sig.df,
             aes(yintercept = sig, linetype = name)) +
  scale_linetype_manual(values = c("BF" = 1, "EIGEN" = 3)) +
  facet_grid(~CHROM, scales = "free_x") +
  labs(
    x = "Genomic position (Mb)",
    y = expression(-log[10](italic(p)))
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    strip.text.x = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_line()
  )
manhattan_plot


#────────────────────────────────────────────────
# Panel B P x G
#────────────────────────────────────────────────
#  load your trait data file
bs_trait_data <- read_tsv("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/262strains_avgBS_A3567_4_1_gwasmapping.tsv")
genotype_matrix <-read_tsv ("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/nemascan/20260903_1/genotype_matrix.tsv")
avg_bs_inbred_qtl <- read_tsv("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/nemascan/20260903_1/avg_BS/gwa/inbred/avg_BS_inbred_qtl.tsv")
peak_snp <- qtl_chrX$marker   
C_briggsae_relatedness_groups_Sheet1 <- read_csv("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/C.briggsae_relatednessgroups.csv")

###Rename stain in bs_trait_data
bs_trait_data <- bs_trait_data %>%
  rename(strain = Strain) %>%
  mutate(strain = recode(strain, "PB420" = "AF16"))

###Rename the column isotype to strain in the C_briggsae_relatedness_groups_Sheet1
C_briggsae_relatedness_groups_Sheet1 <- C_briggsae_relatedness_groups_Sheet1 %>%
  rename(strain = isotype)  %>%
  rename(relatedness_group = 'relatedness group') %>%
  rename(relatedness_group_color = 'relatedness group color')

#####For QTL on left arm#####:
filtered_df_qtl_X<- genotype_matrix %>%
  filter(CHROM == "X", POS == 16163640)

strain_df_qtl_X <- filtered_df_qtl_X %>%
  pivot_longer(
    cols = -(CHROM:ALT),
    names_to = "strain",
  )

combined_data_qtl_X <- merge(bs_trait_data, strain_df_qtl_X, by = "strain") %>%
  mutate(
    allele = factor(dplyr::recode(as.character(value), "-1" = "REF", "1" = "ALT"), levels = c("REF", "ALT")),
    value = avg_BS
  )

combined_data_qtl_X %>%
  count(allele)
##REF: 206 
##ALT: 56

##Combine datasets
master_data <- left_join(
  combined_data_qtl_X,
  C_briggsae_relatedness_groups_Sheet1,
  by = "strain"
) 

##plot
p_blackdots <- master_data %>%
  filter(!is.na(allele)) %>%
  ggplot(aes(x = allele, y = value)) +
  
  # Boxplot (draw first so points sit on top)
  geom_boxplot(
    fill = "gray90",
    color = "black",
    alpha = 0.8,
    width = 0.5,
    outlier.shape = NA   # avoids duplicate outlier dots
  ) +
  
  geom_jitter(
    fill = "black",
    color = "black",
    position = position_jitter(width = 0.12, height = 0),
    size = 0.45,
    shape = 21,
    alpha = 0.6
  ) +
  # Apply your colors
  #scale_fill_manual(values = COLOR_MAP_1) +
  
  # Labels and theme
  labs(x = "Genotype", y = "Lifetime fecundity") +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 12),
  ) +
  theme_classic()

p_blackdots

wilcox.test(
  value ~ allele,
  data = master_data %>% filter(!is.na(allele))
)

# p-value = 1.619e-06


#────────────────────────────────────────────────
# Panel C gene plot
#────────────────────────────────────────────────
fine_mapping <- read_tsv ("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/nemascan/20260903_1/avg_BS/fine_mapping/inbred/X_14531153_18159901/avg_BS_inbred_X_14531153_18159901.annotated.gwa")
#############Gene plot
plot_gene <- function(fine_mapping, qtl_row, method) {
  
  peak_pos <- as.numeric(qtl_row$peakPOS)
  
  gene_df <- fine_mapping %>%
    mutate(
      effect = case_when(
        CONSEQUENCE == "Intergenic" ~ "Intergenic",
        
        CONSEQUENCE %in% c(
          "Intronic",
          "3'-UTR",
          "5'-UTR",
          "Splice region"
        ) ~ "Genic Non-coding",
        
        CONSEQUENCE %in% c(
          "Missense",
          "Stop gained",
          "Stop lost",
          "Start lost",
          "Nonframeshift deletion"
        ) ~ "Coding-change",
        
        CONSEQUENCE == "Synonymous" ~ "Synonymous",
        
        CONSEQUENCE == "N/A" ~ "Intergenic",
        
        TRUE ~ "Intergenic"
      )
    ) %>%
    filter(!is.na(POS), !is.na(P))
  
  xs <- min(gene_df$POS, na.rm = TRUE)
  xe <- max(gene_df$POS, na.rm = TRUE)
  
  ggplot(
    gene_df,
    aes(
      x = POS / 1e6,
      y = -log10(P),
      color = effect
    )
  ) +
    geom_point(
      size = 0.6,
      alpha = 0.6
    ) +
    geom_vline(
      xintercept = peak_pos / 1e6,
      linetype = 3
    ) +
    scale_color_manual(
      values = c(
        "Intergenic" = "gray70",
        "Genic Non-coding" = "orange",
        "Coding-change" = "magenta3",
        "Synonymous" = "green3"
      ),
      breaks = c(
        "Intergenic",
        "Genic Non-coding",
        "Coding-change",
        "Synonymous"
      ),
      labels = c(
        "1", ## this was "Intergenic"
        "2", ## this was "Genic Non-coding"
        "3", ##this was "Coding-\nchange"
        "4" ##this was "Synonymous"
      ),
      name = "Consequence"
    ) +
    coord_cartesian(
      xlim = c(xs / 1e6, xe / 1e6)
    ) +
    labs(
      x = "Genomic position (Mb)",
      y = expression(-log[10](italic(p)))
    ) +
    theme_bw(base_size = 11) +
    theme(
      text = element_text(family = "Helvetica"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      legend.position = "top",
      panel.grid = element_blank()
    )
}

p_gene <- plot_gene(
  fine_mapping = fine_mapping,
  qtl_row = qtl_chrX,
  method = "inbred"
)

p_gene


# ─────────────────────────────────────────────
# Common publication theme
# ─────────────────────────────────────────────

theme_common <- theme_classic(
  base_size = 11,
  base_family = "Helvetica"
) +
  theme(
    text = element_text(
      family = "Helvetica",
      size = 11,
      color = "black"
    ),
    axis.title = element_text(
      family = "Helvetica",
      size = 12,
      color = "black"
    ),
    axis.text = element_text(
      family = "Helvetica",
      size = 11,
      color = "black"
    ),
    strip.text = element_text(
      family = "Helvetica",
      size = 11,
      color = "black"
    ),
    legend.title = element_text(
      family = "Helvetica",
      size = 12,
      color = "black"
    ),
    legend.text = element_text(
      family = "Helvetica",
      size = 11,
      color = "black"
    ),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_blank()
  )

# ─────────────────────────────────────────────
# Apply consistent formatting to each panel
# ─────────────────────────────────────────────

# Panel A: Manhattan plot
pA <- manhattan_plot +
  theme(
    text = element_text(
      family = "Helvetica",
      size = 12,
      color = "black"
    ),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    strip.text.x = element_text(size = 11),
    legend.position = "none",
    panel.grid = element_blank()
  )

pA

# Panel B: Genotype versus lifetime fecundity
pB <- p_blackdots +
  theme_common +
  theme(
    axis.text.x = element_text(
      family = "Helvetica",
      size = 11,
      color = "black"
    ),
    legend.position = "none"
  )
pB

# Panel C: Fine-mapping / functional consequence plot
pC <- p_gene +
  theme_common +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal"
  )
pC

# ─────────────────────────────────────────────
# Assemble the three panels
# ─────────────────────────────────────────────
library(ggplot2)
library(patchwork)

combined_figure <- (
  pA + pB + pC +
    plot_layout(
      design = "AAAA\nBCCC",
      heights = c(1, 1)
    ) +
    plot_annotation(
      tag_levels = "a"
    )
) &
  theme(
    plot.tag = element_text(
      family = "sans",
      face = "bold",
      size = 12,
      color = "black"
    ),
    plot.tag.position = "topleft"
  )

combined_figure


ggsave("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/plots/figure3.tif",
       plot = combined_figure, width = 7.5, height = 6, dpi = 350, units = "in", bg = "white")





