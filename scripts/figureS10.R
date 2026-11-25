
rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)
library(glue)
library(purrr)
#library(ggbeeswarm)
library(ggrepel)
library(dplyr)
#install.packages("emmeans")
library(emmeans)
library(patchwork)
#install.packages("cowplot")
library(cowplot)
# load relatedness groups 
C_briggsae_relatedness_groups_Sheet1 <- read_csv("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/C.briggsae_relatednessgroups.csv")

#  load your trait data file
bs_trait_data <- read_tsv("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/262strains_avgBS_A3567_4_1_gwasmapping.tsv")
genotype_matrix <-read_tsv ("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/nemascan/20260903_1/genotype_matrix.tsv")
avg_bs_inbred_qtl <- read_tsv("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/nemascan/20260903_1/avg_BS/gwa/inbred/avg_BS_inbred_qtl.tsv")

#load Briggsaw WI information to get substrate info
briggsae_WI_strain_info <- read.csv("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/Cbriggsae_WI_strain_info.csv")


###Rename the column isotype to strain in the C_briggsae_relatedness_groups_Sheet1
C_briggsae_relatedness_groups_Sheet1 <- C_briggsae_relatedness_groups_Sheet1 %>%
  rename(strain = isotype)  %>%
  rename(relatedness_group = 'relatedness group') %>%
  rename(relatedness_group_color = 'relatedness group color')


###Rename stain in bs_trait_data
bs_trait_data <- bs_trait_data %>%
  rename(strain = Strain) %>%
  mutate(strain = recode(strain, "PB420" = "AF16"))


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

master_data_subs <- 
  left_join(
  master_data,
  briggsae_WI_strain_info,
  by = "strain"
)



########### Phenotype-by-Genotype by relatedness##########
#
#. PANEL A
#
########## Phenotype-by-Genotype by relatedness###########
COLOR_MAP_2 <- with(
  master_data %>% distinct(relatedness_group, relatedness_group_color),
  setNames(relatedness_group_color, relatedness_group)
)

#Plot
p_relatedness <- master_data %>%
  filter(!is.na(allele)) %>%
  # Ensure allele is a factor for consistent x-axis ordering
  mutate(allele = factor(allele)) %>% 
  ggplot(aes(x = allele, y = value)) +
  
  geom_boxplot(
    fill = "white",      # White fill looks cleaner with jittered points
    color = "black",
    linewidth = 0.6,     # Slightly thicker lines
    width = 0.4,
    outlier.shape = NA,
    alpha = 0.5
  ) +
  
  geom_jitter(
    aes(
      fill = relatedness_group,
      color = relatedness_group
    ),
    position = position_jitter(width = 0.12, height = 0),
    size = 0.45,
    shape = 21,
    alpha = 0.6
  ) +
  
  scale_fill_manual(
    values = COLOR_MAP_2,
    name = "Relatedness"
  ) +
  
  scale_color_manual(
    values = COLOR_MAP_2,
    guide = "none"
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + # Adds breathing room at top
  
  labs(
    x = "Genotype", 
    y = "Lifetime Fecundity", # Title case
    fill = "Relatedness"      # Clean legend title
  ) +
  
  theme_classic() + # Removes gridlines and border for a clean look
  theme(
    # Text sizing for 1-column or 2-column prints
    axis.text.x = element_text(size = 11, color = "black", face = "italic"), # Italicized for alleles
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6)
  )

# Display plot
p_relatedness

#ANOVA
model <- lm(
  value ~ allele * relatedness_group,
  data = master_data
)

anova(model)


emm <- emmeans(
  model,
  ~ relatedness_group | allele
)

pairs(emm, adjust = "tukey")

########### Phenotype-by-Genotype by geography##########
#
# PANEL B
#
########## Phenotype-by-Genotype by geography###########
COLOR_MAP_3 <- with(
  master_data %>% distinct(geo, geo_color),
  setNames(geo_color, geo)
)

#plot
p_geography <- master_data %>%
  filter(!is.na(allele)) %>%
  # Ensure allele is a factor for consistent x-axis ordering
  mutate(allele = factor(allele)) %>% 
  ggplot(aes(x = allele, y = value)) +
  
  geom_boxplot(
    fill = "white",      
    color = "black",
    linewidth = 0.6,     
    width = 0.4,
    outlier.shape = NA,
    alpha = 0.5
  ) +
  
  geom_jitter(
    aes(
      fill = geo,
      color = geo
    ),
    position = position_jitter(width = 0.12, height = 0),
    size = 0.45,
    shape = 21,
    alpha = 0.6
  ) +
  
  scale_fill_manual(
    values = COLOR_MAP_3,
    name = "geo"
  ) +
  
  scale_color_manual(
    values = COLOR_MAP_3,
    guide = "none"
  ) +
  
  scale_fill_manual(values = COLOR_MAP_3) +
  # Adds a small buffer at the top of the y-axis for cleaner presentation
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + 
  
  labs(
    x = "Genotype", 
    y = "Lifetime fecundity",  # Title Case
    fill = "Geography"  # Clean legend title
  ) +
  
  theme_classic() + 
  theme(
    axis.text.x = element_text(size = 11, color = "black", face = "italic"), 
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6)
  )
p_geography



#ANOVA
model <- lm(
  value ~ allele * geo,
  data = master_data
)

anova(model)


emm <- emmeans(
  model,
  ~ geo | allele
)

pairs(emm, adjust = "tukey")



########### Phenotype-by-Genotype by substrate##########
#
# PANEL C
#
########## Phenotype-by-Genotype by substrate###########
##left

p_substrate <- master_data_subs %>%
  filter(!is.na(allele)) %>%
  # Ensure allele is a factor for consistent x-axis ordering
  mutate(allele = factor(allele)) %>% 
  ggplot(aes(x = allele, y = value)) +
  
  geom_boxplot(
    fill = "white",      
    color = "black",
    linewidth = 0.6,     
    width = 0.4,
    outlier.shape = NA,
    alpha = 0.5
  ) +
  
  geom_jitter(
    aes(
      fill = substrate,
      color = substrate
    ),
    position = position_jitter(width = 0.12, height = 0),
    size = 0.45,
    shape = 21,
    alpha = 0.6
  ) +
  
  guides(
    color = "none"
  ) +
  
  
  #scale_fill_manual(values = COLOR_MAP_3) +
  # Adds a small buffer at the top of the y-axis for cleaner presentation
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + 
  
  labs(
    x = "Genotype", 
    y = "Lifetime fecundity",  # Title Case
    fill = "Substrate"  # Clean legend title
  ) +
  
  theme_classic() + 
  theme(
    axis.text.x = element_text(size = 11, color = "black", face = "italic"), 
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6)
  )
p_substrate

#ANOVA
model <- lm(
  value ~ allele * substrate,
  data = master_data_subs
)

anova(model)


emm <- emmeans(
  model,
  ~ substrate | allele
)

pairs(emm, adjust = "tukey")


  
# ============================================================
# COMMON THEME
# ============================================================

common_theme <- theme_classic(base_family = "Helvetica") +
  theme(
    # All text
    text = element_text(
      family = "Helvetica",
      color = "black"
    ),
    
    # Axis text
    axis.text.x = element_text(
      family = "Helvetica",
      size = 11,
      color = "black",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      family = "Helvetica",
      size = 11,
      color = "black"
    ),
    
    # Axis titles
    axis.title.x = element_blank(),
    axis.title.y = element_text(
      family = "Helvetica",
      size = 12,
      face = "plain"
    ),
    
    # Small legend on the right
    legend.position = "right",
    legend.title = element_text(
      family = "Helvetica",
      size = 12,
      face = "bold"
    ),
    legend.text = element_text(
      family = "Helvetica",
      size = 11
    ),
    
    # Small legend keys
    legend.key.size = unit(0.2, "cm"),
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.2, "cm"),
    
    # Compact legend spacing
    legend.spacing.x = unit(0.05, "cm"),
    legend.spacing.y = unit(0.02, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    
    # Compact plot margins
    plot.margin = margin(
      t = 2,
      r = 2,
      b = 2,
      l = 2
    )
  )
# ============================================================
# APPLY COMMON THEME TO ALL FIVE EXISTING PLOTS
# ============================================================


p_relatedness <- p_relatedness + common_theme
p_geography <- p_geography + common_theme
p_substrate <- p_substrate + common_theme


# ============================================================
# FINAL FIGURE
# ============================================================

final_figure <- (
  wrap_plots(
    p_relatedness,
    p_geography,
    p_substrate,
    ncol = 1
  ) +
    plot_annotation(tag_levels = "a")
) &
  theme(
    plot.tag = element_text(
      family = "Helvetica",
      size = 12,
      face = "bold",
      color = "black"
    )
  )

final_figure


ggsave("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/plots/figureS10.tif",
       plot = final_figure,  width = 7.5, height = 9, dpi = 350, units = "in", bg = "white")

