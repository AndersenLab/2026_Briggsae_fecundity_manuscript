
rm(list = ls())

library(ggplot2)
library(tidyr)
library(ggplot2)
library(dplyr)
library(ggrepel)
#install.packages("cowplot")
library(cowplot)
library(RColorBrewer)
#install.packages("MASS")
library(MASS)
#install.packages("heritability")
library(heritability)
#install.packages("devtools")  
library(devtools)
#devtools::install_github('covaruber/sommer')
library(sommer)
library(ggplot2)


########Open the automation raw data
combined_df_1_2_7<-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/automation_raw_data_A_1_2_7.csv",header = TRUE,sep = ',')
combined_df_3_4_5_6<-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/automation_raw_data_A_3_4_5_6.csv",header = TRUE,sep = ',')

#######Combine the two automation dataframes
combined_df_automation <- bind_rows(combined_df_3_4_5_6, combined_df_1_2_7)

combined_df_analysis <- combined_df_automation %>%
  separate(fname,
           into = c("Assay", "Day", "Strain", "Replicate", "extra1", "extra2"),
           sep = "_",
           remove = FALSE)

combined_df_analysis <- combined_df_analysis %>% 
  rename(count = `Well.1`)

combined_df_analysis <- 
  dplyr::select(combined_df_analysis, Assay, Day, Strain, Replicate, count)

##Selecting only assays 1, 3, 4, 5, 6, 7:
combined_df_analysis <- combined_df_analysis %>%
  filter(Assay %in% c("Assay1","Assay3", "Assay4", "Assay5", "Assay6", "Assay7"))

##Renaming PB420 to CGC2
combined_df_analysis <- combined_df_analysis %>%
  mutate(Strain = recode(Strain, "PB420" = "CGC2"))

Sentinel_Colours <- c(
  "CGC2"   = "#53886C",
  "ECA2666" = "#466EB4",
  "ED3102"  = "#00A0E1",
  "JU3237"  = "#F7D42A",
  "NIC1667" = "#7C1D6F",
  "VX34"    = "#CC799D"
)

######only selecting sentinel strains
sentinel_Strains_automation <- combined_df_analysis %>%
  dplyr::filter(Strain %in% names(Sentinel_Colours))
##1135 rows

##Filtering the replicates that have counts for all 5 days
sentinel_Strains_automation_filtered_5days <- sentinel_Strains_automation %>%
  dplyr::group_by(Strain, Assay, Replicate) %>%
  dplyr::filter(dplyr::n_distinct(Day) == 5) %>%
  dplyr::ungroup()
##1025 rows


filtered_df <- sentinel_Strains_automation_filtered_5days %>%
  group_by(Assay, Strain, Replicate) %>%
  summarise(
    total_count = sum(count, na.rm = TRUE),
    .groups = "drop"
  )
##205 rows

filtered_df$Assay  <- droplevels(factor(filtered_df$Assay))
filtered_df$Strain <- droplevels(factor(filtered_df$Strain))

model_full_filtered <- aov(total_count ~ Assay + Strain, data = filtered_df)

anova_table_filtered <- summary(model_full_filtered)[[1]]
anova_table_filtered

anova_df_filtered <- as.data.frame(anova_table_filtered)

total_ss_filtered <- sum(anova_df_filtered$`Sum Sq`, na.rm = TRUE)

anova_df_filtered$percent_variance <-
  (anova_df_filtered$`Sum Sq` / total_ss_filtered) * 100

anova_df_filtered$percent_variance <-
  round(anova_df_filtered$percent_variance, 2)

anova_df_filtered

##leaving one assay out analysis
assays_filtered <- unique(filtered_df$Assay)

results_filtered <- data.frame(
  removed_assay = character(),
  F_value = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (a in assays_filtered) {
  
  data_subset_filtered <- filtered_df[filtered_df$Assay != a, ]
  
  data_subset_filtered$Assay  <- droplevels(data_subset_filtered$Assay)
  data_subset_filtered$Strain <- droplevels(data_subset_filtered$Strain)
  
  model_tmp_filtered <- aov(total_count ~ Assay + Strain,
                            data = data_subset_filtered)
  
  anova_tmp_filtered <- summary(model_tmp_filtered)[[1]]
  
  F_val <- anova_tmp_filtered["Assay","F value"]
  p_val <- anova_tmp_filtered["Assay","Pr(>F)"]
  
  results_filtered <- rbind(results_filtered,
                            data.frame(
                              removed_assay = a,
                              F_value = F_val,
                              p_value = p_val
                            ))
}

results_filtered


##To plot
results_filtered$removed_assay <- factor(
  results_filtered$removed_assay,
  levels = results_filtered$removed_assay
)

assay_colors_filtered <- brewer.pal(
  n = nlevels(results_filtered$removed_assay),
  "Set2"
)

names(assay_colors_filtered) <- levels(results_filtered$removed_assay)


p1_filtered <- ggplot(
  results_filtered,
  aes(
    x = removed_assay,
    y = F_value
  )
) +
  geom_col(
    width = 0.7,
    fill = "grey80",
    color = "black",
    linewidth = 0.4
  ) +
  geom_text(
    aes(label = round(F_value, 1)),
    vjust = -0.4,
    size = 4,
    family = "Helvetica",
    color = "black"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    x = "Assay removed",
    y = expression(italic(F) * "-value (assay effect)")
  ) +
  theme_classic(
    base_size = 14,
    base_family = "Helvetica"
  ) +
  theme(
    text = element_text(
      family = "Helvetica",
      color = "black",
      face = "plain"
    ),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.line = element_line(
      color = "black",
      linewidth = 0.4
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.4
    ),
    plot.margin = margin(8, 12, 8, 8)
  )

p1_filtered


ggsave("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/plots/figureS3.tif",
       plot = p1_filtered,  width = 7.5, height = 6, units = "in", dpi = 350, bg = "white")
