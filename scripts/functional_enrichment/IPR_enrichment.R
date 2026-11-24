library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
# library(enrichplot)
library(stringr)
library(ape)

# ======================================================================================================================================================================================== #
# QTL: V:14531153-16163640
# ======================================================================================================================================================================================== #
gffCB <- ape::read.gff("../../processed_data/genomic_resources/c_briggsae.QX1410_20250929.csq.gff")

allQX <- gffCB %>%
  dplyr::filter(type=="mRNA") %>%
  dplyr::select(attributes) %>%
  dplyr::filter(grepl("biotype=protein_coding", attributes)) %>%
  tidyr::separate(attributes, into=c("tran","rest"), sep = ";Parent=gene:") %>%
  dplyr::mutate(tran = gsub("ID=", "", tran)) %>%
  tidyr::separate(rest, into = c("QX1410", "blah"), sep = ";sequence_name") %>%
  dplyr::select(tran,QX1410)

ipr <- readr::read_tsv("../../processed_data/functional_enrichment/QX1410_IPR_annotations.tsv", col_names = c("tran", "MD5_digest", "seq_length", "app", "signature_accession", "signature_description", "start", "end", "score", "status", "date", "IPR_accession","IPR_description","GO", "pathways")) %>%
  dplyr::left_join(allQX, by = 'tran') %>%
  dplyr::select(-tran) %>%
  dplyr::select(QX1410,IPR_accession,IPR_description,GO) %>%
  dplyr::filter(IPR_description != "-")

write.table(ipr %>% dplyr::rename(`QX1410 transcript` = QX1410), "../../tables/supplementary/", quote = F, col.names = T, row.name = F, sep = "\t")

#==============================================================================================================================================================================================================================#

# INTERPROSCAN enrichment test - entire pangenome as background

#==============================================================================================================================================================================================================================#
# Define universe & QTL (annotated-only universe) 
univ_genes <- ipr_background_genes
hdr_genes  <- intersect(ws_hdr_priv_genes, ipr_background_genes) 

N <- length(univ_genes)
n <- length(hdr_genes)

# Counts per IPR (k = in universe, x = in QTL subset) 
k_tbl <- all_ipr_background %>%
  dplyr::count(IPR_accession, name = "k")

x_tbl <- all_ipr_background %>%
  dplyr::filter(tran %in% ws_hdr_priv_genes) %>%
  dplyr::count(IPR_accession, name = "x")

desc_tbl <- all_ipr_background %>%
  dplyr::distinct(IPR_accession, IPR_description)

# Hypergeometric enrichment (one-sided)
ipr_enrichment <- k_tbl %>%
  dplyr::left_join(x_tbl, by = "IPR_accession") %>%
  dplyr::mutate(x = tidyr::replace_na(x, 0L)) %>%
  dplyr::mutate(
    pval = stats::phyper(q = x - 1, m = k, n = N - k, k = n, lower.tail = FALSE),
    expected = (n * k) / N, # if IPR genes are randomly distributed, you’d expect this many HDR genes to carry the IPR.
    enrich_ratio = dplyr::if_else(expected > 0, x / expected, NA_real_), # (x HDR with IPR / k background with IPR) / (n HDR genes / N background genes)
    # odds ratio with Haldane–Anscombe correction (adding 0.5 to each cell to avoid infinities)
    OR = {
      a <- x + 0.5                                # HDR & has IPR
      b <- (n - x) + 0.5                          # HDR & no IPR
      c <- (k - x) + 0.5                          # non-HDR & has IPR
      d <- (N - n - (k - x)) + 0.5                # non-HDR & no IPR
      (a / b) / (c / d)
    },
    FDR_p.adjust = stats::p.adjust(pval, method = "BH")
  ) %>%
  dplyr::left_join(desc_tbl, by = "IPR_accession") %>%
  dplyr::mutate(N = N, n = n) %>%
  dplyr::select(IPR_accession, IPR_description, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust) %>%
  dplyr::arrange(FDR_p.adjust, dplyr::desc(enrich_ratio))

ipr_sig <- ipr_enrichment %>%
  dplyr::filter(FDR_p.adjust < 0.05)

ipr_sig %>% dplyr::slice_head(n = 20)

ipr_sig_gene_collapsed <- all_ipr_background %>%
  dplyr::filter(IPR_accession %in% ipr_sig$IPR_accession) %>%
  dplyr::group_by(IPR_accession) %>%
  dplyr::summarise(
    IPR_description = dplyr::first(stats::na.omit(IPR_description)),
    n_genes_HDR = dplyr::n_distinct(tran[tran %in% ws_hdr_priv_genes]),
    genes_HDR   = paste(sort(unique(tran[tran %in% ws_hdr_priv_genes])), collapse = ", "),
    n_genes_all = dplyr::n_distinct(tran),
    genes_all   = paste(sort(unique(tran)), collapse = ", "),
    .groups = "drop") %>%
  dplyr::left_join(ipr_sig %>% dplyr::select(IPR_accession, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust), by = "IPR_accession") %>%
  dplyr::mutate(Region = "Private pangenome")

data_plt_priv <- ipr_sig_gene_collapsed %>% 
  dplyr::filter(IPR_accession != "IPR008164") %>% # Excluding annotation of "Repeat of unknown function XGLTT"
  dplyr::arrange(FDR_p.adjust) %>% 
  dplyr::filter(IPR_accession != "") %>% 
  dplyr::slice_head(n = 20) %>% dplyr::arrange(desc(FDR_p.adjust)) %>% 
  dplyr::mutate(plotpoint = dplyr::row_number())

data_plt_priv_plt <- ipr_sig_gene_collapsed %>% 
  dplyr::filter(IPR_accession != "IPR008164") %>% # Excluding annotation of "Repeat of unknown function XGLTT"
  dplyr::arrange(FDR_p.adjust) %>%
  dplyr::filter(IPR_accession != "") %>%
  dplyr::slice_head(n = 30) %>% dplyr::arrange(n_genes_HDR) %>%
  dplyr::mutate(plotpoint = dplyr::row_number()) %>%
  dplyr::rename(`Gene set` = Region) %>% 
  dplyr::mutate(IPR_description = ifelse(IPR_description == "BTB/POZ domain-containing adapter for CUL3-mediated RhoA degradation protein 1-3", 
                                         "BTB/POZ domain-containing adapter for CUL3-mediated\nRhoA degradation protein 1-3", IPR_description))


# Plotting enrichmeent for genes in QTL
priv_HDR_enrich <- ggplot(data_plt_priv_plt) +
  geom_point(aes(x = n_genes_HDR, y = plotpoint, size = enrich_ratio, fill = -log10(FDR_p.adjust), shape = `Gene set`)) +
  scale_y_continuous(breaks = data_plt_priv_plt$plotpoint, labels = data_plt_priv_plt$IPR_description, name = "", expand = c(0.02,0.02)) +
  scale_shape_manual(values =  c("Private pangenome" = 24)) +
  scale_fill_gradient(low = "yellow", high = "red", breaks = c(round(min(-log10(data_plt_priv_plt$FDR_p.adjust))), round((max(-log10(data_plt_priv_plt$FDR_p.adjust)) + min(-log10(data_plt_priv_plt$FDR_p.adjust))) / 2), round(max(-log10(data_plt_priv_plt$FDR_p.adjust))))) +
  scale_size_continuous(range = c(1, 4), name = "Fold enrichment", breaks = pretty(data_plt_priv_plt$enrich_ratio, n = 4)) +
  coord_cartesian(xlim = c(0, 200)) +
  theme(axis.text.x = element_text(size=9, color='black'),
        axis.text.y = element_text(size=6.5, color='black'),
        axis.title = element_text(size=9, color='black'),
        plot.title = element_blank(),
        legend.title = element_text(size = 9, color='black', hjust = 1),
        legend.text = element_text(size = 9, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.7, 0.25),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        text = element_text(family="Helvetica"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(b = 5, t = 5, r = 5, l = 1, unit = "pt")) +
  guides(
    fill = guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 1),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE),
    shape = guide_legend(nrow = 3, order = 3, title.position = "top", title.hjust = 1, override.aes = list(size = 3))) +
  labs(title = "Enriched IPR terms for wild strain genes in HDRs",  x = "Gene count", size = "Fold enrichment", fill = expression(-log[10]~"(corrected p-value)"))
priv_HDR_enrich

# Save the plot
# ggsave("../../figures/supplementary/private_HDR_enrichment.png", priv_HDR_enrich, width = 7.5, height = 6, dpi = 600)
