#!/usr/bin/env Rscript
rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
  library(ggbeeswarm)
  library(ggrepel)
  library(glue)
  library(ggnewscale)
})

# ============================================================
# USER SETTINGS
# ============================================================
gwas_raw <- read_tsv("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/nemascan/20260903_1/avg_BS/gwa/inbred/avg_BS_inbred.gwa",
  col_types = cols()
)

mapping_file <- read_tsv ("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/262strains_avgBS_A3567_4_1_gwasmapping.tsv")
genotype_matrix_file <- read_tsv("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/Fecundity_allstrains_A3_A4_A5_A6/Publication/Nemascan/20260903_1/genotype_matrix.tsv")
independent_tests_file <- read.table("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/nemascan/20260903_1/independent_tests.txt")
qtl <- read_tsv("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/nemascan/20260903_1/avg_BS/gwa/inbred/avg_BS_inbred_qtl.tsv")
fine_mapping <- read_tsv ("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/processed_data/nemascan/20260903_1/avg_BS/fine_mapping/inbred/X_14531153_18159901/avg_BS_inbred_X_14531153_18159901.annotated.gwa")


# ============================================================
# HELPERS
# ============================================================

safe_read_tsv <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  x <- tryCatch(
    fread(path, data.table = FALSE, sep = "\t", header = TRUE, fill = TRUE),
    error = function(e) NULL
  )
  if (is.null(x) || nrow(x) == 0) return(NULL)
  names(x) <- trimws(names(x))
  x
}

safe_read_csv <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  x <- tryCatch(
    fread(path, data.table = FALSE, sep = ",", header = TRUE, fill = TRUE),
    error = function(e) NULL
  )
  if (is.null(x) || nrow(x) == 0) return(NULL)
  names(x) <- trimws(names(x))
  x
}



read_genotype_matrix <- function(path) {
  x <- safe_read_tsv(path)
  if (is.null(x)) return(NULL)
  names(x) <- trimws(names(x))
  if (!all(c("CHROM", "POS", "REF", "ALT") %in% names(x))) return(NULL)
  x %>%
    mutate(
      CHROM = as.character(CHROM),
      POS = suppressWarnings(as.numeric(POS))
    )
}

standardize_gwa <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  names(df) <- trimws(names(df))
  
  # harmonize chromosome column
  if (!("CHR" %in% names(df))) {
    if ("CHROM" %in% names(df)) df <- dplyr::rename(df, CHR = CHROM)
    if ("Chr"   %in% names(df)) df <- dplyr::rename(df, CHR = Chr)
    if ("chr"   %in% names(df)) df <- dplyr::rename(df, CHR = chr)
  }
  
  # harmonize position column
  if (!("POS" %in% names(df))) {
    if ("BP" %in% names(df)) df <- dplyr::rename(df, POS = BP)
    if ("bp" %in% names(df)) df <- dplyr::rename(df, POS = bp)
    if ("pos" %in% names(df)) df <- dplyr::rename(df, POS = pos)
  }
  
  # harmonize p-value column
  if (!("P" %in% names(df))) {
    if ("p" %in% names(df)) df <- dplyr::rename(df, P = p)
    if ("PVAL" %in% names(df)) df <- dplyr::rename(df, P = PVAL)
    if ("pval" %in% names(df)) df <- dplyr::rename(df, P = pval)
  }
  
  if (!all(c("CHR", "POS", "P") %in% names(df))) return(NULL)
  
  df %>%
    mutate(
      CHR = as.character(CHR),
      POS = suppressWarnings(as.numeric(POS)),
      P   = suppressWarnings(as.numeric(P))
    ) %>%
    filter(!is.na(CHR), !is.na(POS), !is.na(P), P > 0)
}

find_trait_dirs <- function(root_dir) {
  dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[!basename(dirs) %in% c("Reports", "summary_figures")]
  dirs
}

find_gwa_file <- function(trait_dir, method) {
  exact <- file.path(trait_dir, "gwa", method, paste0(basename(trait_dir), "_", method, ".gwa"))
  if (file.exists(exact)) return(exact)
  
  alt1 <- Sys.glob(file.path(trait_dir, "gwa", method, paste0("*_", method, ".gwa")))
  if (length(alt1) > 0) return(alt1[[1]])
  
  alt2 <- Sys.glob(file.path(trait_dir, "gwa", method, "*.gwa"))
  if (length(alt2) > 0) return(alt2[[1]])
  
  NULL
}

find_qtl_file <- function(trait_dir, method) {
  candidates <- c(
    file.path(trait_dir, "gwa", method, paste0(basename(trait_dir), "_", method, "_qtl.tsv")),
    file.path(trait_dir, "gwa", method, paste0("*_", method, "_qtl.tsv")),
    file.path(trait_dir, "gwa", method, "*.qtl.tsv")
  )
  for (pat in candidates) {
    files <- Sys.glob(pat)
    if (length(files) > 0) return(files[[1]])
  }
  return(NULL)
}

read_qtl_and_sort <- function(qtl_file, method) {
  qtl <- safe_read_tsv(qtl_file)
  if (is.null(qtl) || nrow(qtl) == 0) return(NULL)
  
  qtl <- as_tibble(qtl)
  # keep only rows with marker / peakPOS
  qtl <- qtl %>%
    mutate(
      method = method,
      CHROM = as.character(CHROM),
      peakPOS = as.numeric(peakPOS),
      startPOS = as.numeric(startPOS),
      endPOS = as.numeric(endPOS)
    ) %>%
    arrange(
      factor(CHROM, levels = roman_chr_levels),
      peakPOS
    )
  qtl
}

find_finemap_folder <- function(trait_dir, method, qtl_row) {
  if (!all(c("CHROM", "startPOS", "endPOS") %in% names(qtl_row))) return(NULL)
  region_name <- paste0(qtl_row$CHROM, "_", qtl_row$startPOS, "_", qtl_row$endPOS)
  folder <- file.path(trait_dir, "fine_mapping", method, region_name)
  if (!dir.exists(folder)) return(NULL)
  folder
}

find_annotated_gwa <- function(finemap_folder) {
  if (is.null(finemap_folder) || !dir.exists(finemap_folder)) return(NULL)
  files <- Sys.glob(file.path(finemap_folder, "*.annotated.gwa"))
  if (length(files) == 0) return(NULL)
  files[[1]]
}

roman_chr_levels <- c("I", "II", "III", "IV", "V", "X", "MtDNA")

# ============================================================
# BLOCK COLOR PALETTE
# ============================================================

# Define colors for blocks
block_colors <- c(
  "N2" = "#000000",
  "block_1" = "#E41A1C",
  "block_2" = "#377EB8",
  "block_3" = "#4DAF4A",
  "block_4" = "#984EA3",
  "block_5" = "#FF7F00",
  "block_6" = "#FFFF33",
  "block_7" = "#A65628",
  "block_8" = "#F781BF",
  "block_9" = "#999999",
  "block_10" = "#66C2A5",
  "block_11" = "#FC8D62",
  "block_12" = "#8DA0CB",
  "unknown" = "#CCCCCC"
)



# ============================================================
# LD 
# ===========================================================


plot_ld <- function(annot_df, method = "inbred") {
  
  if (is.null(annot_df) || nrow(annot_df) == 0 ||
      !all(c("POS", "P", "LD", "SNP") %in% names(annot_df))) {
    return(ggplot() + theme_void() + ggtitle("Missing required data"))
  }
  
  # Convert columns to numeric if necessary
  annot_df$POS <- as.numeric(as.character(annot_df$POS))
  annot_df$P   <- as.numeric(as.character(annot_df$P))
  annot_df$LD  <- as.numeric(as.character(annot_df$LD))
  
  # Identify the most significant SNP
  valid_rows <- which(!is.na(annot_df$P) & annot_df$P > 0)
  
  if (length(valid_rows) == 0) {
    return(ggplot() + theme_void() + ggtitle("No valid P-values"))
  }
  
  peak_row <- annot_df[valid_rows[which.min(annot_df$P[valid_rows])], ]
  peak_pos <- peak_row$POS[1]
  
  ggplot(annot_df, aes(x = POS / 1e6, y = -log10(P), fill = LD)) +
    geom_point(shape = 23, size = 2, na.rm = TRUE) +
    geom_point(
      data = peak_row,
      shape = 23,
      size = 2,
      fill = "red"
    ) +
    geom_vline(
      xintercept = peak_pos / 1e6,
      linetype = 3,
      color = "black"
    ) +
    scale_fill_distiller(
      palette = "YlGnBu",
      direction = 1,
      name = expression(R^2)
    ) +
    labs(
      title = paste0(toupper(method), " LD: ", peak_row$SNP[1]),
      x = "Genomic position (Mb)",
      y = expression(-log[10](italic(p)))
    ) +
    theme_bw(base_size = 12, base_family = "Helvetica") +
    theme(
      text = element_text(family = "Helvetica", color = "black"),
      plot.title = element_blank(),
        #family = "Helvetica",
       # hjust = 0.5,
        #face = "bold"
      #),
      axis.title = element_text(family = "Helvetica"),
      axis.text = element_text(family = "Helvetica"),
      legend.title = element_text(family = "Helvetica"),
      legend.text = element_text(family = "Helvetica"),
      panel.grid = element_blank()
    )
}

p_ld <- plot_ld(
  annot_df = fine_mapping,
  method = "inbred"
)

p_ld

ggsave("/vast/eande106/projects/Nikita/briggsae_GWAS_2025/2026_Fecundity_manuscript/plots/figureS9.tif",
       plot = p_ld, width = 7.5, height = 2.5, dpi = 350, units = "in", bg = "white")
