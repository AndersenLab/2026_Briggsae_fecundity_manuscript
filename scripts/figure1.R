

rm(list = ls())

#install.packages("readlxl")
library(readxl)
library(dplyr)
library(ggplot2)
#install.packages("ggthemes")
library(ggthemes)
#library(maps)
library(ggpubr)
library(readr)
#install.packages("rnaturalearth") 
#install.packages("sf")
#install.packages("rnaturalearthdata")
library(rnaturalearth)
library(sf)
library(rnaturalearthdata)
library(cowplot)

world <- ne_countries(scale = "medium", returnclass = "sf")


##Open the briggsae file
all_data_C.briggsae_Strains <- read.csv("/home/njhaver4/vast-eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/Cbriggsae_WI_strain_info.csv")
all_data_C.briggsae_Strains$latitude<-ifelse(all_data_C.briggsae_Strains$latitude == "NA",NA,all_data_C.briggsae_Strains$latitude)
  
###Open the Briggsae gwas data files:
data_1 <- read_tsv("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/311strains_avgBS_A3567_4_1_.tsv")
data_repwise <- read.delim("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/311strains_avgBS_A3567_4_1_repwise.tsv")


##########Panel a #########
##                       ##
##                       ##
##                       ##
##                       ##
##########Panel a #########


data <- data_1 %>%
  select(-n_reps)

C.briggsae_Strains <- all_data_C.briggsae_Strains %>%
  #dplyr::filter(strain == isotype)  %>% 
  dplyr::select(strain,latitude,longitude)%>%
  dplyr::filter(!is.na(latitude) & !is.na(longitude)) %>%
  dplyr::rename("lat"="latitude", "long"="longitude") 
colnames(C.briggsae_Strains)

C.briggsae_Strains$lat<-as.numeric(C.briggsae_Strains$lat)
C.briggsae_Strains$long<-as.numeric(C.briggsae_Strains$long)

#

############# Geographic distribution of 311 strains########
data <- data %>%
  dplyr::rename("strain"="Strain") 

data_lat_long <- C.briggsae_Strains %>%
  filter(strain %in% data$strain)

#Plot the 311 strains

plot_world_311 <- ggplot() +
    geom_sf(
    data = world,
    fill = "#F7F7F7",
    color = "black",
    linewidth = 0.2
  ) +
  
  geom_point(
    data = data_lat_long,
    aes(x = long, y = lat),
    shape = 16,
    color = "#FF0000",
    size = 0.9,
    alpha = 0.8
  ) +
  
  coord_sf(
    xlim = c(-180, 191),
    ylim = c(-58, 84),
    expand = FALSE
  ) +
  
  theme_void() +
  
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.4, "cm"),
    #panel.border = element_rect(
     # color = "black",
     # fill = NA,
     # linewidth = 0.8
 # )
  )

plot_world_311

plot_world_311_a <- plot_world_311 + annotate(
  "text",
  x = -175,
  y = 80,
  label = "a",
  fontface = "bold",
  size = 5,
  hjust = 0
)


##########Panel a #########
##                       ##
##                       ##
##                       ##
##                       ##
##########Panel a #########

###Having only two colors:
### build color map with two colors:
# Define colors
divergent_col <- "#0072B2"
nondivergent_col <- "#E69F00"

df_summary <- data_repwise %>%
  group_by(Strain) %>%
  summarise(
    Div = first(Div),
    mean_fec = mean(total_count, na.rm = TRUE),
    sd_fec = sd(total_count, na.rm = TRUE),
    .groups = "drop"
  )

strain_order <- df_summary %>%
  arrange(mean_fec) %>%
  pull(Strain)


df_summary <- df_summary %>%
  mutate(
    Strain = factor(
      Strain,
      levels = strain_order
    )
  )

p <- ggplot(
  df_summary,
  aes(
    x = Strain,
    y = mean_fec,
    fill = Div
  )
) +
  geom_col(
    width = 0.85,
    color = "white"
  ) +
  geom_errorbar(
    aes(
      ymin = mean_fec - sd_fec,
      ymax = mean_fec + sd_fec
    ),
    width = 0.25,
    linewidth = 0.4,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Div" = divergent_col,
      "NoDiv" = nondivergent_col
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.10))
  ) +
  scale_x_discrete(
    expand = c(0, 0)
  ) +
  theme_classic(
    base_family = "Helvetica"
  ) +
  theme(
    text = element_text(
      family = "Helvetica"
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  ) +
  labs(
    x = "Strain",
    y = "Lifetime fecundity"
  )

p



final_plot <- cowplot::plot_grid(
  plot_world_311, p,
  align = "v",
  nrow = 2,
  labels = c("a", "b"),
  label_size = 12
)

final_plot

ggsave("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/plots/figure1.tif",
       plot = final_plot, width = 7.5, height = 5, dpi = 350, units = "in", bg = "white")



