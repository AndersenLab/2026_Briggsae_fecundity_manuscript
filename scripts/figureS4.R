
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
#install.packages("pwr")
#library(pwr)
library(ggplot2)
library(patchwork)
library(ggplot2)
library(patchwork)


combined_df_1_2_7<-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/automation_raw_data_A_1_2_7.csv",header = TRUE,sep = ',')
combined_df_3_4_5_6<-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/automation_raw_data_A_3_4_5_6.csv",header = TRUE,sep = ',')

saved_replicates <- read.csv("/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/saved_replicates.csv", header = TRUE,sep = ',')



#############Step 1: Rename specific strains:#############
## JU2801_1 to JU2801 
## GXW002 to GXW22 (The Nemascan mapping pipeline does not recognize GXW0022. It needs GXW22.)
combined_df_3_4_5_6 <- combined_df_3_4_5_6 %>%
  mutate(fname = gsub("JU2801_1", "JU2801", fname)) %>%
  mutate(fname = gsub("GXW0022", "GXW22", fname))

combined_df_1_2_7 <- combined_df_1_2_7 %>%
  #mutate(fname = gsub("Assay1_Day5_JU2801_r_series0000frame00000_t_0.00_s.png", "Assay1_Day5_JU2801_r1_series0000frame00000_t_0.00_s.png", fname)) %>%
  mutate(fname = gsub("JU2801_1", "JU2801", fname)) %>%
  mutate(fname = gsub("GXW0022",  "GXW22",  fname)) %>%
  mutate(fname = gsub("JU2801_2", "JU2801", fname))


#### Bind the two dataframes together:
combined_df <- bind_rows(combined_df_3_4_5_6, combined_df_1_2_7)

##### Separating the column into Assay/Day/Strain/Replicate
combined_df_analysis <- combined_df %>%
  separate(fname,
           into = c("Assay", "Day", "Strain", "Replicate", "extra"),
           sep = "_",
           extra = "merge",
           remove = FALSE) %>%
  rename(count = `Well.1`) %>%
  dplyr::select(Assay, Day, Strain, Replicate, count)

#############Step 2: Rename specific rows:#############
## r2' to r2 (the r2 folder for A4 4 QG2883 was named "r2'" by mistake. This will fix the error)
## r9 to r8 (the r8 folder for A5 QG700 was named r9 by mistake. This will fix the error)
##r to r1(the r1 folder for A1 JU2801 was named "r" by mistake. This will fix the error)

combined_df_analysis <- combined_df_analysis %>%
  mutate(Replicate = recode(Replicate, "r2'" = "r2")) %>%
  mutate(Replicate = recode(Replicate, "r9" = "r8")) %>%
  mutate(Replicate = recode(Replicate, "r" = "r1")) %>%
  mutate(Day = as.numeric(gsub("Day", "", Day)))

##A3 NIC187 r3 is present in two rows 
##A3 BRC20294 is present in two rows
##remove one row of each
combined_df_analysis <- combined_df_analysis %>%
  group_by(Assay, Strain, Replicate, Day) %>%
  slice(1) %>%
  ungroup()


########################Making the master dataset########################
##Making a dataset for assays 1,3,4,5,6,7
########################Making the master dataset########################

df_filtered <- combined_df_analysis %>%
  filter(
    Assay %in% c("Assay1", "Assay3", "Assay4", "Assay5", "Assay6", "Assay7")
  )

##Total count per strain per replicate 
total_counts <- df_filtered %>%
  group_by(Assay, Strain, Replicate) %>%
  summarise(
    total_count = sum(count, na.rm = TRUE),   # total across all days
    n_days = n_distinct(Day),                  # how many days recorded
    .groups = "drop"
  )

unique_strains_step12 <- total_counts %>%
  distinct(Strain)
##366 unique strains


#############Step 3: Keeping replicates only have counts for all 5 days:#############
df_filtered_5days <- df_filtered %>%
  group_by(Assay, Strain, Replicate) %>%
  filter(n_distinct(Day) == 5) %>%
  ungroup()

df_filtered_5days <- df_filtered %>%
  group_by(Assay, Strain, Replicate) %>%
  filter(n_distinct(Day) == 5) %>%
  ungroup()


total_counts_5days <- df_filtered_5days %>%
  group_by(Assay, Strain, Replicate) %>%
  summarise(
    total_count = sum(count, na.rm = TRUE),   # total across all days
    n_days = n_distinct(Day),                  # how many days recorded
    .groups = "drop"
  )


unique_strains_step3 <- total_counts_5days %>%
  distinct(Strain)

##Average BS and SD of strains
avg_BS <- total_counts_5days %>%
  group_by(Strain, Assay) %>%
  summarise(
    avg_total_count = mean(total_count, na.rm = TRUE),
    sd_total_count = sd(total_count, na.rm = TRUE),
    n_replicates = sum(!is.na(total_count)),
    .groups = "drop"
  )


#############Step 4: Saved replicates analysis#############
##Starting from saving replicates that were stuck in the water droplets from A5.
##Want to analyse how the BS of the saved replicates differs from the average BS of the strain
##Keep the replicate only if the z score is between 1.5 and -1.5
##saved replicate analysis

####Total BS for the saved strains
saved_strain_bs <- saved_replicates %>%
  left_join(
    dplyr::select(total_counts_5days, Assay, Strain, Replicate, total_count),
    by = c("Assay", "Strain", "Replicate")
  )


###18 replicates of the 118 have NA meaning that those replicates were lost
##because they didn't have counts for all 5 days
####Remove these 18 replicates from further analysis
clean_df_saved_Strain_bs <- saved_strain_bs %>%
  filter(!is.na(total_count))
###100 replicates left

###Now I want to know how many unique strains are there
unique_strains <- clean_df_saved_Strain_bs %>%
  distinct(Strain)
##72 unique strains

####Join the average bs values for each strain to the clean_df
clean_df_saved_BS_with_avg <- clean_df_saved_Strain_bs %>%
  left_join(avg_BS, by = c("Strain", "Assay")) %>%
  mutate(diff = avg_total_count - total_count)

##Calculate Z score:
clean_df_saved_BS_with_avg <- clean_df_saved_BS_with_avg %>%
  mutate(
    replicate_z = (total_count - avg_total_count) / sd_total_count
  )




#############Step 5: BS>30 on day 5 analysis#############

df_day5_BSmorethan30 <- df_filtered_5days %>%
  group_by(Assay, Strain, Replicate) %>%
  filter(any(Day == "5" & count > 30)) %>%
  ungroup()

# total BS per replicate
replicate_totals_BS_morethan30_day5 <- df_day5_BSmorethan30 %>%
  group_by(Strain, Replicate, Assay) %>%
  summarise(total_count = sum(count, na.rm = TRUE), .groups = "drop")

####Join the average bs values for each strain to the clean_df
clean_df_BSmorethan30day5_with_avg <- replicate_totals_BS_morethan30_day5 %>%
  left_join(avg_BS, by = c("Strain", "Assay")) %>%
  mutate(diff = avg_total_count - total_count)

##Calculate Z score:
clean_df_BSmorethan30day5_with_avg <- clean_df_BSmorethan30day5_with_avg %>%
  mutate(
    replicate_z = (total_count - avg_total_count) / sd_total_count
  )




strain_levels <- sort(unique(c(
  as.character(clean_df_saved_BS_with_avg$Strain),
  as.character(clean_df_BSmorethan30day5_with_avg$Strain)
)))

clean_df_saved_BS_with_avg$Strain <- factor(
  clean_df_saved_BS_with_avg$Strain,
  levels = strain_levels
)

clean_df_BSmorethan30day5_with_avg$Strain <- factor(
  clean_df_BSmorethan30day5_with_avg$Strain,
  levels = strain_levels
)


# --------------------------------------------------
# Common y-axis
# --------------------------------------------------

y_max <- max(
  c(
    clean_df_saved_BS_with_avg$replicate_z,
    clean_df_BSmorethan30day5_with_avg$replicate_z
  ),
  na.rm = TRUE
)

# Set lower limit to -2.5
y_min <- -2.5

# Round upper limit outward
y_max <- ceiling(y_max)



# --------------------------------------------------
# Panel A: 
# --------------------------------------------------

pA <- ggplot(
  clean_df_saved_BS_with_avg,
  aes(x = Strain, y = replicate_z)
) +
  geom_jitter(
    width = 0,
    height = 0,
    alpha = 0.8,
    color = "black",
    size = 1.0
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = c(-1.5, 1.5),
    linewidth = 0.4
  ) +
  scale_y_continuous(
    breaks = pretty(c(y_min, y_max), n = 6)
  ) +
  coord_cartesian(
    ylim = c(y_min, y_max)
  ) +
  theme_classic(base_family = "Helvetica") +
  labs(
    x = "Strain",
    y = "z-score"
  ) +
  theme(
    text = element_text(family = "Helvetica", color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(family = "Helvetica", size = 12),
    axis.title.y = element_text(family = "Helvetica", size = 12),
    axis.text.y = element_text(
      family = "Helvetica",
      size = 11,
      color = "black"
    )
  )

pA


# --------------------------------------------------
# Panel B: 
# --------------------------------------------------

pB <- ggplot(
  clean_df_BSmorethan30day5_with_avg,
  aes(x = Strain, y = replicate_z)
) +
  geom_jitter(
    width = 0,
    height = 0,
    alpha = 0.8,
    color = "black",
    size = 1.0
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = c(-1.5, 1.5),
    linewidth = 0.4
  ) +
  scale_y_continuous(
    breaks = pretty(c(y_min, y_max), n = 6)
  ) +
  coord_cartesian(
    ylim = c(y_min, y_max)
  ) +
  theme_classic(base_family = "Helvetica") +
  labs(
    x = "Strain",
    y = "z-score"
  ) +
  theme(
    text = element_text(family = "Helvetica", color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(family = "Helvetica", size = 12),
    axis.title.y = element_text(family = "Helvetica", size = 12),
    axis.text.y = element_text(
      family = "Helvetica",
      size = 11,
      color = "black"
    ),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

pB


# --------------------------------------------------
# Combine panels
# --------------------------------------------------
combined_plot <- (
  wrap_plots(
    pA,
    pB,
    ncol = 2
  ) +
    plot_annotation(tag_levels = "a")
) &
  theme(
    text = element_text(family = "Helvetica", color = "black"),
    plot.tag = element_text(
      family = "Helvetica",
      face = "bold",
      size = 12,
      color = "black"
    )
  )

combined_plot


ggsave("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/plots/figureS4.tif",
       plot = combined_plot,  width = 7.5, height = 2.5, dpi = 350, units = "in", bg = "white")


