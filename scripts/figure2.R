

# Making figure 2:
##Panel A: box plot of relatedness for 311 strains
##Panel B: box plot of geography for 311 strains
##Panel C: box plot of substrate for 311 strains
##Panel D: scatter plot of abs latitude for 311 strains
##Panel E: scatter plot of abs longitude for 311 strains
##Panel F: scatter plot of abs temperature for 311 strains


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
library(ggplot2)
library(patchwork)
library(grid)

data <- read.delim("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/311strains_avgBS_A3567_4_1_.tsv")
data_repwise <- read.delim("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/311strains_avgBS_A3567_4_1_repwise.tsv")


#load Briggsaw WI information to get substrate info
briggsae_WI_strain_info <- read.csv("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/raw_data/Cbriggsae_WI_strain_info.csv")


###Rename stain in briggsae_WI_strain_info
briggsae_WI_strain_info <- briggsae_WI_strain_info #%>%
  #mutate(strain = recode(strain, "CGC2" = "PB420"))

###Rename stain in bs_trait_data
data <- data %>%
  rename(strain = Strain) #%>%
 # mutate(strain = recode(strain, "CGC2" = "PB420"))

data_repwise <- data_repwise %>%
  rename(strain = Strain) #%>%
 # mutate(strain = recode(strain, "CGC2" = "PB420"))

##Combine datasets
master_data <- left_join(
  data_repwise,
  briggsae_WI_strain_info,
  by = "strain"
) 

master_data_1 <- left_join(
  data,
  briggsae_WI_strain_info,
  by = "strain"
) 

# Summarize lifetime fecundity by strain
df_summary <- data_repwise %>%
  group_by(strain) %>%
  summarise(
    Relatedness_group = first(relatedness_group),
    Relatedness_group_color = first(relatedness_group_color),
    substrate = first(substrate),
    geo = first(geo),
    geo_color = first(geo_color),
    mean_fec = mean(total_count, na.rm = TRUE),
    sd_fec = sd(total_count, na.rm = TRUE),
    .groups = "drop"
  )

# Order strains by mean lifetime fecundity
strain_order <- df_summary %>%
  arrange(mean_fec) %>%
  pull(strain)

# Apply strain order
df_summary <- df_summary %>%
  mutate(
    Strain = factor(
      strain,
      levels = strain_order
    )
  )

##Add relatedness group info for PB420
df_summary <- df_summary %>%
  mutate(
    Relatedness_group = if_else(
      strain == "PB420",
      "Tropical",
      Relatedness_group
    ),
    Relatedness_group_color = if_else(
      strain == "PB420",
      "#ff0000",
      Relatedness_group_color
    )
  )

##Add geography info for PB420
df_summary <- df_summary %>%
  mutate(
    geo = if_else(
      strain == "PB420",
      "Asia",
      geo
    ),
    geo_color = if_else(
      strain == "PB420",
      "#399EB8",
      geo_color
    )
  )

# ============================================================
# Panel A: bar plot with relatednessgroup
# ============================================================


# Create a named color vector from the summary data
relatedness_colors <- df_summary %>%
  distinct(Relatedness_group, Relatedness_group_color) %>%
  tibble::deframe()

p_a <- ggplot(
  df_summary,
  aes(
    x = factor(strain, levels = strain_order),
    y = mean_fec,
    fill = Relatedness_group
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
    width = 0.15, #.25
    linewidth = 0.25, #.4
    color = "black"
  ) +
  scale_fill_manual(
    values = relatedness_colors
  ) +
  
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.10))
  ) +
  #scale_y_continuous(
   # expand = expansion(mult = c(0, 0.10))
 # ) +
  scale_x_discrete(
    expand = c(0, 0)
  ) +
  theme_classic(
    base_family = "Helvetica"
  ) +
  theme(
    text = element_text(
      family = "Helvetica",
      color = "black"
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "Strain",
    y = "Fecundity"
  )

p_a

# One-way ANOVA
model <- aov(mean_fec ~ Relatedness_group, data = df_summary)

# ANOVA results
summary(model)
# ============================================================
# Panel B: bar plot with geography
# ============================================================
geography_colors <- df_summary %>%
  distinct(geo, geo_color) %>%
  tibble::deframe()

p_b <- ggplot(
  df_summary,
  aes(
    x = factor(strain, levels = strain_order),
    y = mean_fec,
    fill = geo
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
    width = 0.15, #.25
    linewidth = 0.25, #.4
    color = "black"
  ) +
  scale_fill_manual(
    values = geography_colors
  ) +
  scale_y_continuous(
    limits = c(0, NA),
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
      family = "Helvetica",
      color = "black"
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "Strain",
    y = "Fecundity"
  )

p_b

# One-way ANOVA
model <- aov(mean_fec ~ geo, data = df_summary)

# ANOVA results
summary(model)

# ============================================================
# Panel C: bar plot with substrate
# ============================================================
df_summary <- df_summary %>%
  mutate(
    substrate_num = factor(
      match(substrate, sort(unique(substrate)))
    )
  )


p_c <- ggplot(
  df_summary,
  aes(
    x = factor(strain, levels = strain_order),
    y = mean_fec,
    fill = substrate_num
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
    width = 0.15, #.25
    linewidth = 0.25, #.4
    color = "black"
  ) +
  #scale_fill_manual(
    #values = geography_colors
  #) +
  scale_y_continuous(
    limits = c(0, NA),
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
      family = "Helvetica",
      color = "black"
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "Strain",
    y = "Fecundity"
  )

p_c

# One-way ANOVA
model <- aov(mean_fec ~ substrate, data = df_summary)

# ANOVA results
summary(model)

# ============================================================
# Panel D: scatter plot with abs latitude
# ============================================================
# Prepare data using absolute latitude
cor_data_lat <- master_data_1 %>%
  filter(!is.na(latitude), !is.na(avg_BS)) %>%
  mutate(abs_latitude = abs(latitude))

# Spearman correlation using absolute latitude
spearman_lat <- cor.test(
  cor_data_lat$abs_latitude,
  cor_data_lat$avg_BS,
  method = "spearman",
  exact = FALSE
)

# Spearman's rho
rho_lat <- unname(spearman_lat$estimate)

# P-value
p_lat <- spearman_lat$p.value

cat("Spearman's rho =", round(rho_lat, 4), "\n")
cat("P-value =", format.pval(p_lat, digits = 4), "\n")


# Linear regression model
model_lat <- lm(avg_BS ~ abs_latitude, data = cor_data_lat)

# Extract R-squared
r2 <- summary(model_lat)$r.squared

# Create plot

# ============================================================
# Panel D: scatter plot with absolute latitude
# ============================================================


# Create plot
p_d <- ggplot(
  cor_data_lat,
  aes(x = abs_latitude, y = avg_BS)
) +
  geom_point(
    size = 0.5,
    shape = 21,
    fill = "black",
    color = "black",
    alpha = 0.7
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    color = "black",
    linewidth = 0.25
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    #label = sprintf("R² = %.3f", r2),
    label = sprintf("Spearman's \u03C1 = %.3f", rho_lat),
    hjust = 1.1,
    vjust = 1.5,
    family = "Helvetica",
    size = 3,
    color = "black"
  ) +
  labs(
    x = "Absolute latitude",
    y = "Fecundity"
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(
    text = element_text(
      family = "Helvetica",
      color = "black"
    ),
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    axis.title = element_text(
      size = 12,
      color = "black"
    ),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    aspect.ratio = 1,
    plot.margin = margin(2, 2, 2, 2)
  )

# Display plot
p_d



# ============================================================
# Panel E: scatter plot with longitude
# ============================================================

# Prepare data for longitude analysis
cor_data_lon <- master_data_1 %>%
  filter(!is.na(longitude), !is.na(avg_BS))


# Spearman correlation 
spearman_long <- cor.test(
  cor_data_lon$longitude,
  cor_data_lon$avg_BS,
  method = "spearman",
  exact = FALSE
)

# Spearman's rho
rho_lon <- unname(spearman_long$estimate)

# P-value
p_lon <- spearman_long$p.value

cat("Spearman's rho =", round(rho_lat, 4), "\n")
cat("P-value =", format.pval(p_lat, digits = 4), "\n")

# Linear regression model
model_long <- lm(avg_BS ~ longitude, data = cor_data_lon)

# Extract R-squared
r2_long <- summary(model_long)$r.squared

# Create plot
p_e <- ggplot(
  cor_data_lon,
  aes(x = longitude, y = avg_BS)
) +
  
  # Data points
  geom_point(
    size = 0.5,
    shape = 21,
    fill = "black",
    color = "black",
    alpha = 0.7
  ) +
  
  # y = x reference line (dashed grey)
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.25
  ) +
  
  # Best-fit linear regression line (solid black)
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    color = "black",
    linewidth = 0.25
  ) +
  
  # R-squared annotation
  annotate(
    "text",
    x = Inf,
    y = Inf,
    #label = sprintf("R² = %.3f", r2_long),
    label = sprintf("Spearman's \u03C1 = %.3f", abs(rho_lon)),
    hjust = 1.1,
    vjust = 1.5,
    family = "Helvetica",
    size = 3,
    color = "black"
  ) +
  
  # Axis labels
  labs(
    x = "Longitude",
    y = "Fecundity"
  ) +
  
  # Publication-style theme
  theme_classic(base_family = "Helvetica") +
  theme(
    text = element_text(
      family = "Helvetica",
      color = "black"
    ),
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    axis.title = element_text(
      size = 12,
      color = "black"
    ),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    aspect.ratio = 1,
    plot.margin = margin(2, 2, 2, 2)
  )

# Display plot
p_e

# ============================================================
# Panel F: scatter plot with temperature
# ============================================================

# Prepare data for substrate temperature analysis
cor_data_temp <- master_data_1 %>%
  filter(!is.na(substrate_temp), !is.na(avg_BS))

# Spearman correlation 
spearman_temp <- cor.test(
  cor_data_temp$substrate_temp,
  cor_data_temp$avg_BS,
  method = "spearman",
  exact = FALSE
)

# Spearman's rho
rho_temp <- unname(spearman_temp$estimate)

# P-value
p_temp <- spearman_temp$p.value

cat("Spearman's rho =", round(rho_temp, 4), "\n")
cat("P-value =", format.pval(p_temp, digits = 4), "\n")

# Linear regression model
model_temp <- lm(avg_BS ~ substrate_temp, data = cor_data_temp)

# Extract R-squared
r2_temp <- summary(model_temp)$r.squared

# Create plot
p_temp <- ggplot(
  cor_data_temp,
  aes(x = substrate_temp, y = avg_BS)
) +
  
  # Data points
  geom_point(
    size = 0.5,
    shape = 21,
    fill = "black",
    color = "black",
    alpha = 0.7
  ) +
  
  # y = x reference line (dashed grey)
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.25
  ) +
  
  # Best-fit linear regression line (solid black)
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    color = "black",
    linewidth = 0.25
  ) +
  
  # R-squared annotation
  annotate(
    "text",
    x = Inf,
    y = Inf,
    #label = sprintf("R² = %.3f", r2_temp),
    label = sprintf("Spearman's \u03C1 = %.3f", abs(rho_temp)),
    hjust = 1.1,
    vjust = 1.5,
    family = "Helvetica",
    size = 3,
    color = "black"
  ) +
  
  # Axis labels
  labs(
    x = "Substrate temperature",
    y = "Fecundity"
  ) +
  
  # Publication-style theme
  theme_classic(base_family = "Helvetica") +
  theme(
    text = element_text(
      family = "Helvetica",
      color = "black"
    ),
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    axis.title = element_text(
      size = 12,
      color = "black"
    ),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    aspect.ratio = 1,
    plot.margin = margin(2, 2, 2, 2)
  )

# Display plot
p_temp



# ============================================================
# Compact theme for the combined figure
# ============================================================

small_theme <- theme(
  text = element_text(
    family = "Helvetica",
    color = "black"
  ),
  
  # Axis text and titles
  axis.text.x = element_text(size = 11, color = "black"),
  axis.text.y = element_text(size = 11, color = "black"),
  axis.title.x = element_text(size = 11, color = "black"),
  axis.title.y = element_text(size = 11, color = "black"),
  
  # Smaller, more compact legends
  legend.text = element_text(size = 11, color = "black"),
  legend.title = element_blank(),
  legend.key.size = unit(0.25, "cm"),
  legend.key.width = unit(0.25, "cm"),
  legend.key.height = unit(0.25, "cm"),
  legend.spacing.y = unit(0.02, "cm"),
  legend.spacing.x = unit(0.05, "cm"),
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(0, 0, 0, 0),
  legend.box.spacing = unit(0.05, "cm"),
  
  # Reduce whitespace within each panel
  plot.margin = margin(2, 2, 2, 2),
  
  # Thinner axes and ticks
  axis.line = element_line(linewidth = 0.4),
  axis.ticks = element_line(linewidth = 0.4)
)

# ============================================================
# Apply compact formatting to each panel
# ============================================================

# Panels A–C: remove x-axis labels and ticks
p_a2 <- p_a + small_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
    #axis.title.x = element_blank()
  )  +
  labs(
    x = "Strain"
  )

p_b2 <- p_b + small_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
 #   axis.title.x = element_blank()
  )  +
  labs(
    x = "Strain"
  )

p_c2 <- p_c +
  small_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    x = "Strain"
  )

# Remove fixed square aspect ratio from scatter plots
p_d2 <- p_d + small_theme +
  theme(aspect.ratio = NULL)

p_e2 <- p_e + small_theme +
  theme(aspect.ratio = NULL)

p_f2 <- p_temp + small_theme +
  theme(aspect.ratio = NULL)

# ============================================================
# Make the legend appear on the bottom of each panel
# ============================================================

# Panel A: Relatedness group
p_a2 <- p_a2 +
  guides(
    fill = guide_legend(nrow = 2, byrow = TRUE),
    color = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(
      family = "helvetica",
      size = 11,
      color = "black"
    ),
    legend.title = element_blank(),
    legend.key.size = unit(0.2, "cm"),
    legend.spacing.x = unit(0.08, "cm"),
    legend.spacing.y = unit(0.02, "cm"),
    legend.box = "horizontal",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )

# Panel B: Geography
p_b2 <- p_b2 +
  guides(
    fill = guide_legend(nrow = 2, byrow = TRUE),
    color = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(
      family = "helvetica",
      size = 11,
      color = "black"
    ),
    legend.title = element_blank(),
    legend.key.size = unit(0.2, "cm"),
    legend.spacing.x = unit(0.08, "cm"),
    legend.spacing.y = unit(0.02, "cm"),
    legend.box = "horizontal",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )

# Panel C: Substrate
p_c2 <- p_c2 +
  guides(
    fill = guide_legend(nrow = 2, byrow = TRUE),
    color = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(
      family = "helvetica",
      size = 11,
      color = "black"
    ),
    legend.title = element_blank(),
    legend.key.size = unit(0.2, "cm"),
    legend.spacing.x = unit(0.08, "cm"),
    legend.spacing.y = unit(0.02, "cm"),
    legend.box = "horizontal",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )


# ============================================================
# Combine all six panels
# ============================================================

# Combine all six panels
combined_plot <- wrap_plots(
  A = p_a2,
  B = p_b2,
  C = p_c2,
  D = p_d2,
  E = p_e2,
  F = p_f2,
  design = "
AA
BB
CC
DE
F.
",
  guides = "keep"
) +
  plot_annotation(
    tag_levels = "a"
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

# Display figure
combined_plot
ggsave("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/plots/figure2.tif",
       plot = combined_plot, width = 7.5, height = 9, dpi = 350, units = "in", bg = "white")

