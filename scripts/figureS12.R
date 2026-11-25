##Using 4.4.1 R libs
##To compare counts of sentinel trains obtained by the BIP and students counting
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
#install.packages("readr")
library(readr)
library(dplyr)
library(tidyr)

#install.packages("purrr")  
library(purrr)
library(stringr)

########Open the automation raw data
combined_df_1_2_7<-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/automation_raw_data_A_1_2_7.csv",header = TRUE,sep = ',')
combined_df_3_4_5_6<-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/automation_raw_data_A_3_4_5_6.csv",header = TRUE,sep = ',')

#######Combine the two automation dataframes
combined_df_automation <- bind_rows(combined_df_3_4_5_6, combined_df_1_2_7)

#######Open the manual count data
#A1
assay_1_D1_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/Assay1_sentenialstrains_Day1.csv", header = TRUE, sep = ',')
assay_1_D2_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay1_sentenialstrains_Day2.csv", header = TRUE, sep = ',')
assay_1_D3_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/Assay1_sentenialstrains_Day3.csv", header = TRUE, sep = ',')
assay_1_D4_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/Assay1_sentenialstrains_Day4.csv", header = TRUE, sep = ',')
assay_1_D5_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/Assay1_sentenialstrains_Day5.csv", header = TRUE, sep = ',')


#A3
assay_3_D1_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/Assay3_sentinelstrains_Day1.csv", header = TRUE, sep = ',')
assay_3_D2_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay3_sentinelstrains_Day2.csv", header = TRUE, sep = ',')
assay_3_D3_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay3_sentinelstrains_Day3.csv", header = TRUE, sep = ',')
assay_3_D4_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay3_sentinelstrains_Day4.csv", header = TRUE, sep = ',')
assay_3_D5_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay3_sentinelstrains_Day5.csv", header = TRUE, sep = ',')


#A4
assay_4_D1_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/Assay4_sentenialstrains_Day1.csv", header = TRUE, sep = ',')
assay_4_D2_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay4_sentenialstrains_Day2.csv", header = TRUE, sep = ',')
assay_4_D3_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay4_sentenialstrains_Day3.csv", header = TRUE, sep = ',')
assay_4_D4_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay4_sentenialstrains_Day4.csv", header = TRUE, sep = ',')
assay_4_D5_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay4_sentenialstrains_Day5.csv", header = TRUE, sep = ',')


#A5
assay_5_D1_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay5_sentinelstrains_Day1.csv", header = TRUE, sep = ',')
assay_5_D2_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay5_sentinelstrains_Day2.csv", header = TRUE, sep = ',')
assay_5_D3_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay5_sentinelstrains_Day3.csv", header = TRUE, sep = ',')
assay_5_D4_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay5_sentinelstrains_Day4.csv", header = TRUE, sep = ',')
assay_5_D5_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay5_sentinelstrains_Day5.csv", header = TRUE, sep = ',')


#A6
assay_6_D1_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay6_sentinelstrains_Day1.csv", header = TRUE, sep = ',')
assay_6_D2_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay6_sentinelstrains_Day2.csv", header = TRUE, sep = ',')
assay_6_D3_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay6_sentinelstrains_Day3.csv", header = TRUE, sep = ',')
assay_6_D4_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay6_sentinelstrains_Day4.csv", header = TRUE, sep = ',')
assay_6_D5_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay6_sentinelstrains_Day5.csv", header = TRUE, sep = ',')

#A7
assay_7_D1_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay7_sentinelstrains_Day1.csv", header = TRUE, sep = ',')
assay_7_D2_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay7_sentinelstrains_Day2.csv", header = TRUE, sep = ',')
assay_7_D3_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay7_sentinelstrains_Day3.csv", header = TRUE, sep = ',')
assay_7_D4_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay7_sentinelstrains_Day4.csv", header = TRUE, sep = ',')
assay_7_D5_manual <-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025//2026_Fecundity_manuscript/raw_data/Assay7_sentinelstrains_Day5.csv", header = TRUE, sep = ',')


##Taking the average of the three students for each assay and each day
##Making a master sheet
assays <- c(1, 3, 4, 5, 6, 7)
days <- 1:5

master_manual <- do.call(rbind, lapply(assays, function(a) {
  do.call(rbind, lapply(days, function(d) {
    
    df_name <- paste0("assay_", a, "_D", d, "_manual")
    df <- get(df_name)
    
    df %>%
      dplyr::select(
        sample,
        total_offspring_1,
        total_offspring_1.1,
        total_offspring_1.2
      ) %>%
      dplyr::mutate(
        dplyr::across(
          starts_with("total_offspring"),
          ~ na_if(tolower(trimws(as.character(.))), "n/a")
        ),
        dplyr::across(
          starts_with("total_offspring"),
          as.numeric
        )
      ) %>%
      tidyr::pivot_longer(
        cols = starts_with("total_offspring"),
        names_to = "rep_measure",
        values_to = "Count"
      ) %>%
      dplyr::filter(!is.na(Count)) %>%
      tidyr::separate(
        sample,
        into = c("Strain", "Assay", "Day", "Replicate"),
        sep = "_"
      ) %>%
      dplyr::mutate(
        Day = stringr::str_remove(Day, "^D"),
        Replicate = stringr::str_remove(Replicate, "^R"),
        Assay = paste0(
          "Assay",
          stringr::str_remove(Assay, "^[A-Za-z]+")
        )
      ) %>%
      dplyr::group_by(Strain, Assay, Day, Replicate) %>%
      dplyr::summarise(
        Count = mean(Count),
        .groups = "drop"
      )
  }))
}))

master_manual

master_manual <- master_manual %>%
  mutate(Strain = recode(Strain, "PB420" = "CGC2"))

master_manual <- master_manual %>%
  dplyr::mutate(Day = as.numeric(Day)) %>%
  dplyr::mutate(Replicate = as.numeric(Replicate))
##1135 rows

##Filtering the replicates that have counts for all 5 days
master_manual_filtered_5days <- master_manual %>%
  dplyr::group_by(Strain, Assay, Replicate) %>%
  dplyr::filter(dplyr::n_distinct(Day) == 5) %>%
  dplyr::ungroup()
##1025 rows

Sentinel_Colours <- c(
  "CGC2"   = "#53886C",
  "ECA2666" = "#466EB4",
  "ED3102"  = "#00A0E1",
  "JU3237"  = "#F7D42A",
  "NIC1667" = "#7C1D6F",
  "VX34"    = "#CC799D"
)


####Now moving to automation data to get the same dataframe
###separating the column
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
  filter(Assay %in% c("Assay1", "Assay3", "Assay4" , "Assay5", "Assay6", "Assay7"))


combined_df_analysis <- combined_df_analysis %>%
  mutate(Strain = recode(Strain, "PB420" = "CGC2"))

##Selecting only sentinel strains
sentinel_Strains_automation <- combined_df_analysis %>%
  dplyr::filter(Strain %in% names(Sentinel_Colours))
##1135 rows

sentinel_Strains_automation <- sentinel_Strains_automation %>%
  dplyr::mutate(Replicate = as.numeric(stringr::str_remove(Replicate, "^r")))

sentinel_Strains_automation <- sentinel_Strains_automation %>%
  dplyr::rename(Count = count)

##Filtering the replicates that have counts for all 5 days
sentinel_Strains_automation_filtered_5days <- sentinel_Strains_automation %>%
  dplyr::group_by(Strain, Assay, Replicate) %>%
  dplyr::filter(dplyr::n_distinct(Day) == 5) %>%
  dplyr::ungroup()
##1025 rows

##master_manual_filtered_5days and sentinel_Strains_automation_filtered_5days 
##both have 1025 rows. which is good

##Making same column names for the automation and manual counts
manual <- master_manual_filtered_5days %>%
  rename(
    Assay = Assay,                # manual had lowercase 'assay'
    Strain = Strain,
    Replicate = Replicate,
    count_manual = Count
  )

auto <- sentinel_Strains_automation_filtered_5days %>%
  rename(
    Assay = Assay,                # adjust if needed
    Strain = Strain,
    Replicate = Replicate,
    count_automation = Count
    #n_days = n_days               # keep if you need it; otherwise drop
  )
auto <- auto %>%
  mutate(Day = as.numeric(str_remove(Day, "^Day")))

manual <- manual %>%
  mutate(
    Assay = as.character(Assay),
    Strain = as.character(Strain),
    Replicate = as.character(Replicate),
    Day = as.numeric(Day)
  )

auto <- auto %>%
  mutate(
    Assay = as.character(Assay),
    Strain = as.character(Strain),
    Replicate = as.character(Replicate),
    Day = as.numeric(Day)
  )

##Combining the auto and manual dataframes
combined <- left_join(manual, auto,
                      by = c("Assay", "Strain", "Replicate", "Day"))



###First count total fecundity across 5 days:
##Total count per strain per replicate for auto
total_counts_automation <- sentinel_Strains_automation_filtered_5days %>%
  group_by(Assay, Strain, Replicate) %>%
  summarise(
    total_count = sum(Count, na.rm = TRUE),   # total across all days
    n_days = n_distinct(Day),                  # how many days recorded
    .groups = "drop"
  )

##Total count per strain per replicate for manual
total_counts_manual <- master_manual_filtered_5days %>%
  group_by(Assay, Strain, Replicate) %>%
  summarise(
    total_count = sum(Count, na.rm = TRUE),   # total across all days
    n_days = n_distinct(Day),                  # how many days recorded
    .groups = "drop"
  )


##Making same column names for the automation and manual counts
auto_total <- total_counts_automation %>%
  rename(
    Assay = Assay,                # adjust if needed
    Strain = Strain,
    Replicate = Replicate,
    count_automation = total_count
    #n_days = n_days               # keep if you need it; otherwise drop
  )

manual_total <- total_counts_manual %>%
  rename(
    Assay = Assay,                # manual had lowercase 'assay'
    Strain = Strain,
    Replicate = Replicate,
    count_manual = total_count
  )

##Combine the auto and manual dataframes
  combined_total <- left_join(manual_total, auto_total,
                      by = c("Assay", "Strain", "Replicate", "n_days"))

#pivot for the plot
plot_combined_total <- combined_total %>%
  pivot_longer(cols = c(count_manual, count_automation),
               names_to = "Method",
               values_to = "Total_Count") %>%
  mutate(
    Shape = ifelse(Method == "count_manual", 16, 17), # 16 = circle, 17 = triangle
    Method = factor(Method, levels = c("count_manual", "count_automation")),
    Replicate = factor(Replicate)
  )

#####To make a correlation plot:
#First get mean for each strain in each assay

df_mean <- combined_total %>%
  group_by(Strain, Assay) %>%
  summarise(
    mean_manual = mean(count_manual, na.rm = TRUE),
    mean_auto = mean(count_automation, na.rm = TRUE),
    .groups = "drop"
  )


##Know the min and max nubmer of replicates
replicate_summary <- combined_total %>%
  group_by(Strain, Assay) %>%
  summarise(
    n_replicates = n_distinct(Replicate),
    .groups = "drop"
  )

#Assign strain colors
Sentinel_Colours <- c(
  "CGC2"   = "#53886C",
  "ECA2666" = "#466EB4",
  "ED3102"  = "#00A0E1",
  "JU3237"  = "#F7D42A",
  "NIC1667" = "#7C1D6F",
  "VX34"    = "#CC799D"
)

# Fit the regression
fit <- lm(mean_auto ~ mean_manual, data = df_mean)

# Extract slope, intercept, and R-squared
slope <- coef(fit)["mean_manual"]
intercept <- coef(fit)["(Intercept)"]
r2 <- summary(fit)$r.squared

# Create label text
#label_text <- paste0("y = ", round(slope, 2), "x + ", round(intercept,1),
# "\nR² = ", round(r2, 2))

##Label only for R²
label_text <- paste0("R² = ", round(r2, 2))

# Determine axis limits with buffer
axis_range <- range(c(df_mean$mean_manual, df_mean$mean_auto), na.rm = TRUE)
buffer <- 0.05 * diff(axis_range)
axis_limits <- c(axis_range[1] - buffer, axis_range[2] + buffer)

####Bland altman plot
# Calculate mean and difference per replicate
bland_altman_df <- combined_total %>%
  mutate(
    mean_count = (count_manual + count_automation) / 2,
    diff_count = count_automation - count_manual
  )

# Calculate mean difference and 95% limits of agreement
mean_diff <- mean(bland_altman_df$diff_count, na.rm = TRUE)
sd_diff <- sd(bland_altman_df$diff_count, na.rm = TRUE)
loa_upper <- mean_diff + 1.96 * sd_diff
loa_lower <- mean_diff - 1.96 * sd_diff


#### Correlation plot

p11 <- ggplot(
  df_mean,
  aes(
    x = mean_manual,
    y = mean_auto,
    color = Strain,
    shape = Assay
  )
) +
  geom_point(
    size = 1.5,
    alpha = 0.8
  ) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    se = FALSE,
    color = "black",
    linewidth = 0.5
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.5
  ) +
  scale_color_manual(
    values = Sentinel_Colours
  ) +
  coord_fixed(
    xlim = axis_limits,
    ylim = axis_limits
  ) +
  labs(
    x = "Manual Count",
    y = "Automation Count",
    color = "Strain",
    shape = "Assay"
  ) +
  annotate(
    "text",
    x = axis_limits[1] + 0.05 * diff(axis_limits),
    y = axis_limits[2] - 0.05 * diff(axis_limits),
    label = label_text,
    hjust = 0,
    vjust = 1,
    size = 4
  ) +
  theme_classic(
    base_family = "Helvetica"
  ) +
  theme(
    text = element_text(
      family = "Helvetica",
      size = 12
    ),
    axis.title = element_text(
      family = "Helvetica",
      size = 12
    ),
    axis.text = element_text(
      family = "Helvetica",
      size = 11
    ),
    legend.title = element_text(
      family = "Helvetica",
      size = 12
    ),
    legend.text = element_text(
      family = "Helvetica",
      size = 11
    )
  )
p11

#### Bland–Altman plot

bland_altman_df <- combined_total %>%
  mutate(
    mean_count = (count_manual + count_automation) / 2,
    diff_count = count_automation - count_manual
  )

mean_diff <- mean(
  bland_altman_df$diff_count,
  na.rm = TRUE
)

sd_diff <- sd(
  bland_altman_df$diff_count,
  na.rm = TRUE
)

loa_upper <- mean_diff + 1.96 * sd_diff
loa_lower <- mean_diff - 1.96 * sd_diff


p_ba <- ggplot(
  bland_altman_df,
  aes(
    x = mean_count,
    y = diff_count,
    color = Strain,
    shape = Assay
  )
) +
  geom_point(
    size = 1.5,
    alpha = 0.8
  ) +
  geom_hline(
    yintercept = mean_diff,
    color = "red",
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = loa_upper,
    color = "blue",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = loa_lower,
    color = "blue",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_color_manual(
    values = Sentinel_Colours
  ) +
  labs(
    x = "Mean of Manual and Automation Counts",
    y = "Difference (Automation - Manual)",
    color = "Strain",
    shape = "Assay"
  ) +
  theme_classic(
    base_family = "Helvetica"
  ) +
  theme(
    text = element_text(
      family = "Helvetica",
      size = 12
    ),
    axis.title = element_text(
      family = "Helvetica",
      size = 12
    ),
    axis.text = element_text(
      family = "Helvetica",
      size = 11
    ),
    legend.title = element_text(
      family = "Helvetica",
      size = 12
    ),
    legend.text = element_text(
      family = "Helvetica",
      size = 11
    ),
    legend.position = "none"
  )

p_ba

library(patchwork)

combined_plot <- (p11 / free(p_ba)) +
  plot_annotation(
    tag_levels = "a"
  ) &
  theme(
    plot.tag = element_text(
      family = "Helvetica",
      size = 12,
      face = "bold"
    )
  )

combined_plot

ggplot2::ggsave(
  filename = "/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/plots/figureS12.tif",
  plot = combined_plot,
  width = 7.5,      # in inches
  height = 7.5,     # in inches
  dpi = 350       # high resolution
)

