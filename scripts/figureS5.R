
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
#library(lme4)
#library(sommer)
#library(MASS)
library(dplyr)
#install.packages("readr")
library(readr)
#install.packages("pwr")
#library(pwr)#install.packages("purrr")   # run once if not installed
#library(purrr)
library(stringr)


combined_df_1_2_7<-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/automation_raw_data_A_1_2_7.csv",header = TRUE,sep = ',')
combined_df_3_4_5_6<-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/automation_raw_data_A_3_4_5_6.csv",header = TRUE,sep = ',')

###Open the list of 311 strains:
df_311_strains_list <- read_tsv("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/311strains_avgBS_A3567_4_1_.tsv")
#311 strains")

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



df_filtered <- combined_df_analysis %>%
  semi_join(df_311_strains_list, by = "Strain")

##renaming the day column
df_filtered <- df_filtered %>%
  mutate(Day = as.numeric(str_extract(Day, "\\d+"))) %>%
  filter(Day %in% 1:5)

df_filtered %>%
  distinct(Strain)

##removing assay 2
df_filtered_noassay2 <- df_filtered %>%
  filter(Assay != "Assay2")

###Remove the replicates that do not have counts for all 5 days
##You have to do this because you are pulling all replicates from
##all the 247 strains
##You want to now remove the replicates that do not have counts for all 5 days
df_filtered_5daysonly <- df_filtered_noassay2 %>%
  group_by(Strain, Replicate) %>%
  filter(n_distinct(Day) == 5) %>%
  ungroup()

df_filtered_5daysonly %>%
  distinct(Strain)
##311 strains

Divergent_Colours <- c(
  "QG1005"   = "#0072B2",
  "BRC20530" = "#0072B2",
  "ECA2670"  = "#0072B2",
  "QG4233"   = "#0072B2",
  "JU2801"   = "#0072B2",
  "QG4094"   = "#0072B2",
  "JU3203"   = "#0072B2",
  "QG4063"   = "#0072B2",
  "BRC20503" = "#0072B2",
  "ECA2618"  = "#0072B2",
  "JU1341"   = "#0072B2",
  "NIC19"    = "#0072B2",
  "QG4208"   = "#0072B2",
  "QR25"     = "#0072B2",
  "QG4095"   = "#0072B2",
  "HPT24"    = "#0072B2",
  "BRC20524" = "#0072B2",
  "ECA2621"  = "#0072B2",
  "JU3200"   = "#0072B2",
  "JU3206"   = "#0072B2",
  "HPT11"    = "#0072B2",
  "JU2767"   = "#0072B2",
  "ECA163"   = "#0072B2",
  "VSL2209"  = "#0072B2",
  "QG4131"   = "#0072B2",
  "HPT18"    = "#0072B2",
  "VSL2207"  = "#0072B2",
  "QG4214"   = "#0072B2",
  "JU3202"   = "#0072B2",
  "QG2893"   = "#0072B2",
  "QG2919"   = "#0072B2",
  "BRC20115" = "#0072B2",
  "NIC1607"  = "#0072B2",
  "QG4064"   = "#0072B2"
)


Others <- "#E69F00"


##Average replicates per day for each strain
df_summary_perday <- df_filtered_5daysonly %>%
  group_by(Strain, Day) %>%
  summarise(mean_BS = mean(count, na.rm = TRUE), .groups = "drop")

##Which day does each strain hav the highest BS
df_peak_day <- df_summary_perday %>%
  group_by(Strain) %>%
  slice_max(order_by = mean_BS, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Strain, Peak_Day = Day)

##join 
df_plot <- df_summary_perday %>%
  left_join(df_peak_day, by = "Strain")

df_plot <- df_plot %>%
  mutate(
    Strain_color = ifelse(
      Strain %in% names(Divergent_Colours),
      "#0072B2",
      "#E69F00"
    )
  )


##count # of strains that have highest BS on day 1 vs on day 2
df_peak_counts <- df_peak_day %>%
  count(Peak_Day) %>%
  arrange(Peak_Day)

df_peak_counts


p1 <- ggplot(df_plot, aes(x = Day, y = mean_BS, group = Strain)) +
  geom_line(
    aes(color = Strain_color),
    linewidth = 0.2,
    alpha = 0.4
  ) +
  geom_point(
    aes(color = Strain_color),
    size = 0.25,
    alpha = 0.6
  ) +
  facet_wrap(
    ~Peak_Day,
    nrow = 2,
    labeller = as_labeller(c(
      `1` = "Peak reproduction day 1",
      `2` = "Peak reproduction day 2"
    ))
  ) +
  scale_color_identity() +
  scale_y_continuous(limits = c(0, NA)) +
  theme_classic(base_family = "Helvetica") +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    strip.background = element_blank(),
    strip.text = element_text(size = 12)
  ) +
  labs(
    x = "Day",
    y = "Lifetime fecundity"
  )

p1

ggsave("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/plots/FigureS5.tif",
       plot = p1,  width = 7.5, height = 6, units = "in", dpi = 350,
      bg = "white")

