##Using 4.4.1 R libs
##Script used to generate supp table 1 for the briggsae BS manuscript

##NOTE: Assays 1, 2, 3 are the same set of strains
##NOTE: Assays 4 and 7 are the same set of strains

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


########Open file for assays 1, 2, and 7
combined_df_1_2_7<-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/automation_raw_data_A_1_2_7.csv",header = TRUE,sep = ',')
combined_df_3_4_5_6<-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/automation_raw_data_A_3_4_5_6.csv",header = TRUE,sep = ',')

##Open file that has information about saved replicates:
saved_replicates <- read.csv("/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/saved_replicates.csv", header = TRUE,sep = ',')

##Open the file that has information about the problematic strains 
##The file also has info about substrate, geography, etc. 
briggsae_WI<-read.csv(file = "/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/Cbriggsae_WI_strain_info.csv",header = TRUE,sep = ',')

briggsae_WI <- briggsae_WI %>%
  rename(Strain = strain)

briggsae_problematicstrains <- briggsae_WI %>%
  filter(issues == TRUE)
##68 strains

briggsae_problematicstrains <- briggsae_problematicstrains %>%
  rename(Strain = strain)

##Open the file having Ryan's best mapping panel list of strains
Ryan_25_outliers <- read.table("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/cb_allout25_strains.txt", header = FALSE, sep = "\t")

##Rename GXW0022 to GXW22
Ryan_25_outliers <- Ryan_25_outliers %>%
  mutate(V1 = recode(V1, "GXW0022" = "GXW22"))

# Open file that has relatedness groups information
C_briggsae_relatedness_groups_Sheet1 <- read_csv("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/C.briggsae_relatednessgroups.csv")

C_briggsae_relatedness_groups_Sheet1 <- C_briggsae_relatedness_groups_Sheet1 %>%
  rename(Strain = isotype)


Sentinel_Colours <- c(
  "PB420"   = "#53886C",
  "ECA2666" = "#466EB4",
  "ED3102"  = "#00A0E1",
  "JU3237"  = "#F7D42A",
  "NIC1667" = "#7C1D6F",
  "VX34"    = "#CC799D"
)

Divergent_Colours <- c(
  "QG1005"   = "#E24A33",
  "BRC20530" = "#E24A33",
  "ECA2670"  = "#E24A33",
  "QG4233"   = "#E24A33",
  "JU2801"   = "#E24A33",
  "QG4094"   = "#E24A33",
  "JU3203"   = "#E24A33",
  "QG4063"   = "#E24A33",
  "BRC20503" = "#E24A33",
  "ECA2618"  = "#E24A33",
  "JU1341"   = "#E24A33",
  "NIC19"    = "#E24A33",
  "QG4208"   = "#E24A33",
  "QR25"     = "#E24A33",
  "QG4095"   = "#E24A33",
  "HPT24"    = "#E24A33",
  "BRC20524" = "#E24A33",
  "ECA2621"  = "#E24A33",
  "JU3200"   = "#E24A33",
  "JU3206"   = "#E24A33",
  "HPT11"    = "#E24A33",
  "JU2767"   = "#E24A33",
  "ECA163"   = "#E24A33",
  "VSL2209"  = "#E24A33",
  "QG4131"   = "#E24A33",
  "HPT18"    = "#E24A33",
  "VSL2207"  = "#E24A33",
  "QG4214"   = "#E24A33",
  "JU3202"   = "#E24A33",
  "QG2893"   = "#E24A33",
  "QG2919"   = "#E24A33",
  "BRC20115" = "#E24A33",
  "NIC1607"  = "#E24A33",
  "QG4064"   = "#E24A33"
)


Others <- "gray70"
  
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

##Heritability##
herit_df_step12 <- total_counts %>%
  dplyr::rename(
    strain = Strain,
    phenotype = total_count
  )

herit_df_step12$strain <- as.factor(herit_df_step12$strain)
herit_df_step12$Assay  <- as.factor(herit_df_step12$Assay)


##Fit model only for strain:
model_step12 <- mmer(
  phenotype ~ 1,
  random = ~ strain,
  data = herit_df_step12
)

### Extract variance components
vc_step12 <- summary(model_step12)$varcomp

### View model summary (optional check)
summary(model_step12)

Vg_step12 <- vc_step12["strain.phenotype-phenotype", "VarComp"]
Ve_step12 <- vc_step12["units.phenotype-phenotype", "VarComp"]

H2_step12 <- Vg_step12 / (Vg_step12 + Ve_step12)
H2_step12
##0.265

##Fit model for strain + assay
model_1 <- mmer(
  phenotype ~ 1,
  random = ~ strain + Assay,
  data = herit_df_step12
)

##Modeling:
##Phenotype=overall mean+Strain+Assay+Residual
##The model estimates how much of the variation in phenotype comes from:
#Strain: differences among strains
#Assay: differences among assays/experimental batches
#Residual: everything else not explained by strain or assay
#VP=VG+VAssay+VE

### View model summary (optional check)
summary(model_1)

### Extract variance components
vc_1 <- summary(model_1)$varcomp

### Variance components
Vg_1 <- vc_1["strain.phenotype-phenotype", "VarComp"]   # genetic variance
Va_1 <- vc_1["Assay.phenotype-phenotype", "VarComp"]    # assay/batch variance
Ve_1 <- vc_1["units.phenotype-phenotype", "VarComp"]    # residual variance

### Broad-sense heritability
H2 <- Vg_1 / (Vg_1 + Va_1 + Ve_1)
H2
#0.260
##26.0% of the total phenotype variance is attributed
##to differences among strains, with assay and residual variation included in the denominator.


#############Step 3: Keeping replicates only have counts for all 5 days:#############
df_filtered_5days <- df_filtered %>%
  group_by(Assay, Strain, Replicate) %>%
  filter(n_distinct(Day) == 5) %>%
  ungroup()

df_filtered_5days <- df_filtered %>%
  group_by(Assay, Strain, Replicate) %>%
  filter(n_distinct(Day) == 5) %>%
  ungroup()

##1 
total_counts_5days <- df_filtered_5days %>%
  group_by(Assay, Strain, Replicate) %>%
  summarise(
    total_count = sum(count, na.rm = TRUE),   # total across all days
    n_days = n_distinct(Day),                  # how many days recorded
    .groups = "drop"
  )


unique_strains_step3 <- total_counts_5day %>%
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


##Heritability##
herit_df_step3 <- total_counts_5day %>%
  dplyr::rename(
    strain = Strain,
    phenotype = total_count
  )

herit_df_step3$strain <- as.factor(herit_df_step3$strain)
herit_df_step3$Assay  <- as.factor(herit_df_step3$Assay)


##Fit model only for strain:
model_step3 <- mmer(
  phenotype ~ 1,
  random = ~ strain,
  data = herit_df_step3
)

### Extract variance components
vc_step3 <- summary(model_step3)$varcomp

### View model summary (optional check)
summary(model_step3)

Vg_step3 <- vc_step3["strain.phenotype-phenotype", "VarComp"]
Ve_step3 <- vc_step3["units.phenotype-phenotype", "VarComp"]

H2_step3 <- Vg_step3 / (Vg_step3 + Ve_step3)
H2_step3
##0.317

##Fit model for strain + assay
model_3_1 <- mmer(
  phenotype ~ 1,
  random = ~ strain + Assay,
  data = herit_df_step3
)

##Modeling:
##Phenotype=overall mean+Strain+Assay+Residual
##The model estimates how much of the variation in phenotype comes from:
#Strain: differences among strains
#Assay: differences among assays/experimental batches
#Residual: everything else not explained by strain or assay
#VP=VG+VAssay+VE

### View model summary (optional check)
summary(model_3_1)

### Extract variance components
vc_3_1 <- summary(model_3_1)$varcomp

### Variance components
Vg_3_1 <- vc_3_1["strain.phenotype-phenotype", "VarComp"]   # genetic variance
Va_3_1 <- vc_3_1["Assay.phenotype-phenotype", "VarComp"]    # assay/batch variance
Ve_3_1 <- vc_3_1["units.phenotype-phenotype", "VarComp"]    # residual variance

### Broad-sense heritability
H2 <- Vg_3_1 / (Vg_3_1 + Va_3_1 + Ve_3_1)
H2
#0.313
##31.3% of the total phenotype variance is attributed
##to differences among strains, with assay and residual variation included in the denominator.


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


####Now to get rid of the replicates that have a Z score of 
####Greater than 1.5 or less than -1.5

kept_saved_BS_zreplicate  <- clean_df_saved_BS_with_avg %>%
  filter(replicate_z >= -1.5 & replicate_z <= 1.5)

##Kept 84 replicates
unique_strains_kept <- kept_saved_BS_zreplicate %>%
  distinct(Strain)
##63 unique strains

removed_saved_BS_zreplicate  <- clean_df_saved_BS_with_avg %>%
  filter(replicate_z > 1.5 | replicate_z < -1.5)
##Removed 16 replicates

unique_strains_removed <- removed_saved_BS_zreplicate %>%
  distinct(Strain)

##15 unique strains
 

####Now remove the above replicates from calculation of the average BS for these strains
cleaned_total_counts_saved_analysis <- total_counts_5days %>%
  anti_join(removed_saved_BS_zreplicate, by = c("Assay", "Strain", "Replicate"))


###Removed 2399 -  2384 = 15 lines, (Why not 16?? - this is because one strain:
###A5 JU2868  r2 was saved on day 1 and day 2. so it will be only one line in the total_counts
###because it is at replicate-wise, not day-wise structure)
unique_strains_step4 <- cleaned_total_counts_saved_analysis %>%
  distinct(Strain)
###366 unique strains


##Heritability##
herit_df_step4 <- cleaned_total_counts_saved_analysis %>%
  dplyr::rename(
    strain = Strain,
    phenotype = total_count
  )

herit_df_step4$strain <- as.factor(herit_df_step4$strain)
herit_df_step4$Assay  <- as.factor(herit_df_step4$Assay)


##Fit model only for strain:
model_step4 <- mmer(
  phenotype ~ 1,
  random = ~ strain,
  data = herit_df_step4
)

### Extract variance components
vc_step4 <- summary(model_step4)$varcomp

### View model summary (optional check)
summary(model_step4)

Vg_step4 <- vc_step4["strain.phenotype-phenotype", "VarComp"]
Ve_step4 <- vc_step4["units.phenotype-phenotype", "VarComp"]

H2_step4 <- Vg_step4 / (Vg_step4 + Ve_step4)
H2_step4
##0.325

##Fit model for strain + assay ##This part is not working. it is crashing. so I am not running it
#model_4_1 <- mmer(
 # phenotype ~ 1,
 # random = ~ strain + Assay,
#  data = herit_df_step4
#)

##Modeling:
##Phenotype=overall mean+Strain+Assay+Residual
##The model estimates how much of the variation in phenotype comes from:
#Strain: differences among strains
#Assay: differences among assays/experimental batches
#Residual: everything else not explained by strain or assay
#VP=VG+VAssay+VE

### View model summary (optional check)
#summary(model_4_1)

### Extract variance components
#vc_4_1 <- summary(model_4_1)$varcomp

### Variance components
#Vg_4_1 <- vc_4_1["strain.phenotype-phenotype", "VarComp"]   # genetic variance
#Va_4_1 <- vc_4_1["Assay.phenotype-phenotype", "VarComp"]    # assay/batch variance
#Ve_4_1 <- vc_4_1["units.phenotype-phenotype", "VarComp"]    # residual variance

### Broad-sense heritability
#H2 <- Vg_4_1 / (Vg_4_1 + Va_4_1 + Ve_4_1)
#H2
#
##xx.xxx% of the total phenotype variance is attributed
##to differences among strains, with assay and residual variation included in the denominator.


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


####Now to get rid of the replicates that have a Z score of 
####Greater than 1.5 or less than -1.5
kept_BS_morethan30day5_zreplicate  <- clean_df_BSmorethan30day5_with_avg %>%
  filter(replicate_z >= -1.5 & replicate_z <= 1.5)

##Kept 3 replicate
unique_strains_kept_BSmorehtna30day5 <- kept_BS_morethan30day5_zreplicate %>%
  distinct(Strain)
##3 unique strains

removed_BSmorehtna30day5  <-clean_df_BSmorethan30day5_with_avg %>%
  filter(replicate_z > 1.5 | replicate_z < -1.5)
##Removed 2 replicates

unique_strains_removed_BSmorehtna30day5 <- removed_BSmorehtna30day5 %>%
  distinct(Strain)
##2 unique strains

####Now remove the above replicates from calculation of the average BS for these strains
cleaned_total_counts_saved_BSmorehtan30_analysis <- cleaned_total_counts_saved_analysis %>%
  anti_join(removed_BSmorehtna30day5, by = c("Strain", "Replicate", "total_count"))

###Removed 3847 - 2 = 3845 lines
unique_strains_clean_df_saved_BSmorehtan30 <- cleaned_total_counts_saved_BSmorehtan30_analysis %>%
  distinct(Strain)
###366 unique strains

##Add a column of number of replicates per strain to each row
cleaned_total_counts_saved_BSmorehtan30_analysis <- cleaned_total_counts_saved_BSmorehtan30_analysis %>%
  group_by(Strain, Assay) %>%
  mutate(n_replicates = n_distinct(Replicate)) %>%
  ungroup()


rep_dist <- cleaned_total_counts_saved_BSmorehtan30_analysis %>%
  group_by(Strain, Assay) %>%
  summarise(n_replicates = n_distinct(Replicate), .groups = "drop")

replicate_counts <- rep_dist %>%
  group_by(n_replicates, Assay) %>%
  summarise(n_strains = n(), .groups = "drop")



##Heritability##
herit_df_step5 <- cleaned_total_counts_saved_BSmorehtan30_analysis %>%
  dplyr::rename(
    strain = Strain,
    phenotype = total_count
  )

herit_df_step5$strain <- as.factor(herit_df_step5$strain)
herit_df_step5$Assay  <- as.factor(herit_df_step5$Assay)


##Fit model only for strain:
model_step5 <- mmer(
  phenotype ~ 1,
  random = ~ strain,
  data = herit_df_step5
)

### Extract variance components
vc_step5 <- summary(model_step5)$varcomp

### View model summary (optional check)
summary(model_step5)

Vg_step5 <- vc_step5["strain.phenotype-phenotype", "VarComp"]
Ve_step5 <- vc_step5["units.phenotype-phenotype", "VarComp"]

H2_step5 <- Vg_step5 / (Vg_step5 + Ve_step5)
H2_step5
##0.325

##Fit model for strain + assay ##This step is not running 
#model_5_1 <- mmer(
 # phenotype ~ 1,
  #random = ~ strain + Assay,
 # data = herit_df_step5
#)

##Modeling:
##Phenotype=overall mean+Strain+Assay+Residual
##The model estimates how much of the variation in phenotype comes from:
#Strain: differences among strains
#Assay: differences among assays/experimental batches
#Residual: everything else not explained by strain or assay
#VP=VG+VAssay+VE

### View model summary (optional check)
#summary(model_5_1)

### Extract variance components
#vc_5_1 <- summary(model_5_1)$varcomp

### Variance components
#Vg_5_1 <- vc_5_1["strain.phenotype-phenotype", "VarComp"]   # genetic variance
#Va_5_1 <- vc_5_1["Assay.phenotype-phenotype", "VarComp"]    # assay/batch variance
#Ve_5_1 <- vc_5_1["units.phenotype-phenotype", "VarComp"]    # residual variance

### Broad-sense heritability
#H2 <- Vg_5_1 / (Vg_5_1 + Va_5_1 + Ve_5_1)
#H2
##
##xx.xxx% of the total phenotype variance is attributed
##to differences among strains, with assay and residual variation included in the denominator.

#############Step 6: Only keep strains that have 5 or more replicates#############

df_filtered_processed_data1 <- cleaned_total_counts_saved_BSmorehtan30_analysis %>%
  filter(n_replicates >= 5)

unique_strains_step6 <- df_filtered_processed_data1 %>%
  distinct(Strain)
##321 unique strains


df_filtered_daywise <- df_filtered_5days %>%
  semi_join(df_filtered_processed_data1, by = c("Assay", "Strain", "Replicate"))

##Gut checks to make sure the correct strains, assays, replicates, etc are being pulled
#Gut check1 for total # of strains in the processed_Dataset1 which is 2328
df_filtered_daywise %>%
  summarise(n_strains = n_distinct(Strain))
##321 strains, which is correct

##
df_filtered_daywise %>%
  group_by(Assay) %>%
  summarise(n_rows = n(), .groups = "drop")

##
df_filtered_daywise %>%
  group_by(Assay) %>%
  summarise(n_strains = n_distinct(Strain), .groups = "drop")

#Gut check2 for total #of rows for assay,strain, replicate combo which
#should be equivalent to #of rows in df_filtered_processed_data1
#which is 2966
df_filtered_daywise %>%
  distinct(Assay, Strain, Replicate) %>%
  nrow()
##2966

##Total count per strain per replicate 
total_counts_filtered <- df_filtered_daywise %>%
  group_by(Assay, Strain, Replicate) %>%
  summarise(
    total_count = sum(count, na.rm = TRUE),   # total across all days
    n_days = n_distinct(Day),                  # how many days recorded
    .groups = "drop"
  )

##Heritability##
herit_df_step6 <- total_counts_filtered %>%
  dplyr::rename(
    strain = Strain,
    phenotype = total_count
  )

herit_df_step6$strain <- as.factor(herit_df_step6$strain)
herit_df_step6$Assay  <- as.factor(herit_df_step6$Assay)


##Fit model only for strain:
model_step6 <- mmer(
  phenotype ~ 1,
  random = ~ strain,
  data = herit_df_step6
)

### Extract variance components
vc_step6 <- summary(model_step6)$varcomp

### View model summary (optional check)
summary(model_step6)

Vg_step6 <- vc_step6["strain.phenotype-phenotype", "VarComp"]
Ve_step6 <- vc_step6["units.phenotype-phenotype", "VarComp"]

H2_step6 <- Vg_step6 / (Vg_step6 + Ve_step6)
H2_step6
##0.335

##Average BS and SD of strains
avg_BS_filtered <- total_counts_filtered %>%
  group_by(Strain) %>%
  summarise(
    avg_total_count = mean(total_count, na.rm = TRUE),
    sd_total_count = sd(total_count, na.rm = TRUE),
    n_replicates = sum(!is.na(total_count)),
    .groups = "drop"
  )

###NowPivot Days into columns - will go at the end after all the cleaning
df_wide <- df_filtered_daywise %>%
  dplyr::mutate(day = paste0("d", Day)) %>%
  dplyr::select(-Day) %>%
  tidyr::pivot_wider(
    names_from = day,
    values_from = count,
    values_fill = NA
  )
##This also has 2966 rows, which is exactly what it should have

##Add columns for:
##1. for total_count for each row,
##2. how many days do you have counts for that particular row
##4. if problematic strain or not
##5. if the strain is divergent or not
##6. overlapping with Ryan's best mapping panel set

divergent_strains <- names(Divergent_Colours)


df_wide <- df_wide %>%
  mutate(
    n_days = rowSums(!is.na(across(starts_with("d")))),
    total_count = rowSums(across(starts_with("d")), na.rm = TRUE),
    Div = ifelse(Strain %in% names(Divergent_Colours), "Div", "NoDiv")
  )

df_wide <- df_wide %>%
  mutate(
    problematic = ifelse(
      Strain %in% briggsae_problematicstrains$Strain,
      "Yes",
      "No"
    )
  )

df_wide <- df_wide %>%
  mutate(
    Ryan_mapping_panel = ifelse(
      Strain %in% Ryan_25_outliers$V1,
      "Yes",
      "No"
    )
  )


df_strain_avg <- df_wide %>%
  group_by(Strain, Div) %>%
  summarise(
    avg_BS = mean(total_count, na.rm = TRUE),
    n_reps = n(),
    .groups = "drop"
  )



##having the final processed_Data_1 with all rows:
final_processed_data_1 <- df_wide %>%
  left_join(
    df_filtered_processed_data1 %>%
      dplyr::select(Assay, Strain, Replicate, total_count, n_replicates),
    by = c("Assay", "Strain", "Replicate", "total_count")
  )

##final_processed_data_1 is processed_data_1


#############Step 7: Making the final dataset:#############
####All strains from assay 3, 5, 6, 7
####Strains from assay 1 that are not in assay 3
####Strains from assay 4 that are not in assay 7

unique_assay1_not_assay3 <- final_processed_data_1 %>%
  filter(Assay == "Assay1") %>%
  distinct(Strain) %>%
  filter(!Strain %in% (
    final_processed_data_1 %>%
      filter(Assay == "Assay3") %>%
      pull(Strain)
  ))

##3 strains unique to assay 1

unique_assay4_not_assay7 <- final_processed_data_1 %>%
  filter(Assay == "Assay4") %>%
  distinct(Strain) %>%
  filter(!Strain %in% (
    final_processed_data_1 %>%
      filter(Assay == "Assay7") %>%
      pull(Strain)
  ))
##8 strains unique to assay 4

final_processed_data_2 <- final_processed_data_1 %>%
  filter(
    Assay %in% c("Assay3", "Assay5", "Assay6", "Assay7") |
      
      # Assay 1: only strains NOT in Assay 3
      (Assay == "Assay1" &
         !Strain %in% (final_processed_data_1 %>%
                         filter(Assay == "Assay3") %>%
                         pull(Strain))) |
      
      # Assay 4: only strains NOT in Assay 7
      (Assay == "Assay4" &
         !Strain %in% (final_processed_data_1 %>%
                         filter(Assay == "Assay7") %>%
                         pull(Strain)))
  )


##
final_processed_data_2 %>%
  group_by(Assay) %>%
  summarise(n_rows = n(), .groups = "drop")

##
final_processed_data_2 %>%
  group_by(Assay) %>%
  summarise(n_strains = n_distinct(Strain), .groups = "drop")

unique_strains_step7 <- final_processed_data_2 %>%
  distinct(Strain)
##321 unique strains


##Add columns for latitude, longitude, substrate 
final_processed_data_2 <- final_processed_data_2 %>%
  left_join(
    briggsae_WI %>%
      dplyr::select(Strain, latitude, longitude, substrate, substrate_temp),
    by = "Strain"
  )

##Add columns for relatedness group, geography, abslat, etc
C_briggsae_relatedness_groups_Sheet1 <- C_briggsae_relatedness_groups_Sheet1 %>%
  dplyr::rename(
    relatedness_group = `relatedness group`,
    relatedness_group_color = `relatedness group color`
  )

final_processed_data_2 <- final_processed_data_2 %>%
  left_join(
    C_briggsae_relatedness_groups_Sheet1 %>%
      dplyr::select(Strain, relatedness_group, lat, long, geo, geo_color, abslat, relatedness_group_color),
    by = "Strain"
  )


##This is supplementary table S1
write_tsv(final_processed_data_2, "/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/tableS1.tsv")


##heritability
herit_df_step7 <- final_processed_data_2 %>%
  dplyr::rename(
    strain = Strain,
    phenotype = total_count
  )

herit_df_step7$strain <- as.factor(herit_df_step7$strain)
herit_df_step7$Assay  <- as.factor(herit_df_step7$Assay)


##Fit model only for strain:
model_step7 <- mmer(
  phenotype ~ 1,
  random = ~ strain,
  data = herit_df_step7
)

### Extract variance components
vc_step7 <- summary(model_step7)$varcomp

### View model summary (optional check)
summary(model_step7)

Vg_step7 <- vc_step7["strain.phenotype-phenotype", "VarComp"]
Ve_step7 <- vc_step7["units.phenotype-phenotype", "VarComp"]

H2_step7 <- Vg_step7 / (Vg_step7 + Ve_step7)
H2_step7
##0.342



#############Step 8: Removing problematic strains:#############
df_filtered_removeproblematic <- final_processed_data_2 %>%
  filter(problematic != "Yes")


df_filtered_removeproblematic %>%
  summarise(n_strains = n_distinct(Strain))
##311 unique strains


df_filtered_removeproblematic %>%
  group_by(Assay) %>%
  summarise(n_rows = n(), .groups = "drop")

##
df_filtered_removeproblematic %>%
  group_by(Assay) %>%
  summarise(n_strains = n_distinct(Strain), .groups = "drop")

##Add columns for latitude, longitude, substrate 
df_filtered_removeproblematic <- df_filtered_removeproblematic %>%
  left_join(
    briggsae_WI %>%
      dplyr::select(Strain, latitude, longitude, substrate, substrate_temp),
    by = "Strain"
  )

##Add columns for relatedness group, geography, abslat, etc
C_briggsae_relatedness_groups_Sheet1 <- C_briggsae_relatedness_groups_Sheet1 %>%
  dplyr::rename(
    relatedness_group = `relatedness group`,
    relatedness_group_color = `relatedness group color`
  )

df_filtered_removeproblematic <- df_filtered_removeproblematic %>%
  left_join(
    C_briggsae_relatedness_groups_Sheet1 %>%
      dplyr::select(Strain, relatedness_group, lat, long, geo, geo_color, abslat, relatedness_group_color),
    by = "Strain"
  )

df_strain_avg_step8 <- df_filtered_removeproblematic %>%
  group_by(Strain) %>%
  summarise(
    avg_BS = mean(total_count, na.rm = TRUE),
    n_reps = n(),
    .groups = "drop"
  )


##Heritability##
herit_df_step8 <- df_filtered_removeproblematic %>%
  dplyr::rename(
    strain = Strain,
    phenotype = total_count
  )

herit_df_step8$strain <- as.factor(herit_df_step8$strain)
herit_df_step8$Assay  <- as.factor(herit_df_step8$Assay)


##Fit model only for strain:
model_step8 <- mmer(
  phenotype ~ 1,
  random = ~ strain,
  data = herit_df_step8
)

### Extract variance components
vc_step8 <- summary(model_step8)$varcomp

### View model summary (optional check)
summary(model_step8)

Vg_step8 <- vc_step8["strain.phenotype-phenotype", "VarComp"]
Ve_step8 <- vc_step8["units.phenotype-phenotype", "VarComp"]

H2_step8 <- Vg_step8 / (Vg_step8 + Ve_step8)
H2_step8
##0.337


write_tsv(df_strain_avg_step8, "/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/311strains_avgBS_A3567_4_1_.tsv")
write_tsv(df_filtered_removeproblematic, "/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/311strains_avgBS_A3567_4_1_repwise.tsv")

#############Step 9: Removing divergent strains:#############
df_filtered_removeproblematic_removedivergent <- df_filtered_removeproblematic %>%
  filter(Div != "Div")


df_filtered_removeproblematic_removedivergent %>%
  summarise(n_strains = n_distinct(Strain))
##283 unique strains


df_filtered_removeproblematic_removedivergent %>%
  group_by(Assay) %>%
  summarise(n_rows = n(), .groups = "drop")

##
df_filtered_removeproblematic_removedivergent %>%
  group_by(Assay) %>%
  summarise(n_strains = n_distinct(Strain), .groups = "drop")


##Heritability##
herit_df_step9 <- df_filtered_removeproblematic_removedivergent %>%
  dplyr::rename(
    strain = Strain,
    phenotype = total_count
  )

herit_df_step9$strain <- as.factor(herit_df_step9$strain)
herit_df_step9$Assay  <- as.factor(herit_df_step9$Assay)


##Fit model only for strain:
model_step9 <- mmer(
  phenotype ~ 1,
  random = ~ strain,
  data = herit_df_step9
)

### Extract variance components
vc_step9 <- summary(model_step9)$varcomp

### View model summary (optional check)
summary(model_step9)

Vg_step9 <- vc_step9["strain.phenotype-phenotype", "VarComp"]
Ve_step9 <- vc_step9["units.phenotype-phenotype", "VarComp"]

H2_step9 <- Vg_step9 / (Vg_step9 + Ve_step9)
H2_step9
##0.332


#############Step 10: Removing strains not overlapping with Ryan's best mapping panel:#############
df_filtered_step10 <- df_filtered_removeproblematic_removedivergent %>%
  filter(Ryan_mapping_panel != "No")


df_filtered_step10 %>%
  summarise(n_strains = n_distinct(Strain))
##256 unique strains


df_filtered_step10 %>%
  group_by(Assay) %>%
  summarise(n_rows = n(), .groups = "drop")

##
df_filtered_step10 %>%
  group_by(Assay) %>%
  summarise(n_strains = n_distinct(Strain), .groups = "drop")


##Heritability##
herit_df_step10 <- df_filtered_step10 %>%
  dplyr::rename(
    strain = Strain,
    phenotype = total_count
  )

herit_df_step10$strain <- as.factor(herit_df_step10$strain)
herit_df_step10$Assay  <- as.factor(herit_df_step10$Assay)


##Fit model only for strain:
model_step10 <- mmer(
  phenotype ~ 1,
  random = ~ strain,
  data = herit_df_step10
)

### Extract variance components
vc_step10 <- summary(model_step10)$varcomp

### View model summary (optional check)
summary(model_step10)

Vg_step10 <- vc_step10["strain.phenotype-phenotype", "VarComp"]
Ve_step10 <- vc_step10["units.phenotype-phenotype", "VarComp"]

H2_step10 <- Vg_step10 / (Vg_step10 + Ve_step10)
H2_step10
##0.327


df_strain_avg_step10 <- df_filtered_step10 %>%
  group_by(Strain) %>%
  summarise(
    avg_BS = mean(total_count, na.rm = TRUE),
    n_reps = n(),
    .groups = "drop"
  )
##262 strains
df_strain_avg_step10 <- df_strain_avg_step10 %>%
  dplyr::select(-n_reps)


df_avg_daywise <- df_filtered_step10 %>%
  group_by(Strain) %>%
  summarise(
    d1_avg = mean(d1, na.rm = TRUE),
    d2_avg = mean(d2, na.rm = TRUE),
    d3_avg = mean(d3, na.rm = TRUE),
    d4_avg = mean(d4, na.rm = TRUE),
    d5_avg = mean(d5, na.rm = TRUE),
    .groups = "drop"
  )

df_d1 <- df_avg_daywise %>%
  group_by(Strain) %>%
  summarise(
    avg_count_d1 = mean(d1_avg, na.rm = TRUE),
    .groups = "drop"
  )

df_d2 <- df_avg_daywise %>%
  group_by(Strain) %>%
  summarise(
    avg_count_d2 = mean(d2_avg, na.rm = TRUE),
    .groups = "drop"
  )

df_d3 <- df_avg_daywise %>%
  group_by(Strain) %>%
  summarise(
    avg_count_d3 = mean(d3_avg, na.rm = TRUE),
    .groups = "drop"
  )

df_d4 <- df_avg_daywise %>%
  group_by(Strain) %>%
  summarise(
    avg_count_d4 = mean(d4_avg, na.rm = TRUE),
    .groups = "drop"
  )

df_d5 <- df_avg_daywise %>%
  group_by(Strain) %>%
  summarise(
    avg_count_d5 = mean(d5_avg, na.rm = TRUE),
    .groups = "drop"
  )

##This is the final_dataset that will be used for GWAS mapping
write_tsv(df_strain_avg_step10, "/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/262strains_avgBS_A3567_4_1_gwasmapping.tsv")

write_tsv(df_d1, "/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/day1_262strains_GWASmapping.tsv")
write_tsv(df_d2, "/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/day2_262strains_GWASmapping.tsv")
write_tsv(df_d3, "/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/day3_262strains_GWASmapping.tsv")
write_tsv(df_d4, "/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/day4_262strains_GWASmapping.tsv")
write_tsv(df_d5, "/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/day5_262strains_GWASmapping.tsv")

