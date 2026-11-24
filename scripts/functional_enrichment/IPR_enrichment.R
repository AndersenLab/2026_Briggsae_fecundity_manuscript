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
gffCB <- ape::read.gff("../../processed_data/genomic_resources/c_briggsae.QX1410_20250929.csq.gff") %>% dplyr::filter(seqid != "MtDNA")

allQX <- gffCB %>%
  dplyr::filter(type=="mRNA") %>%
  dplyr::select(chrom = seqid, start, end, attributes) %>%
  dplyr::filter(grepl("biotype=protein_coding", attributes)) %>%
  tidyr::separate(attributes, into=c("tran","rest"), sep = ";Parent=gene:") %>%
  dplyr::mutate(tran = gsub("ID=", "", tran)) %>%
  tidyr::separate(rest, into = c("QX1410", "blah"), sep = ";sequence_name") %>%
  dplyr::select(tran,QX1410, chrom, start, end)

ipr <- readr::read_tsv("../../processed_data/functional_enrichment/QX1410_IPR_annotations.tsv", col_names = c("tran", "MD5_digest", "seq_length", "app", "signature_accession", "signature_description", "start", "end", "score", "status", "date", "IPR_accession","IPR_description","GO", "pathways")) %>%
  dplyr::left_join(allQX %>% dplyr::select(-chrom, -start, -end), by = 'tran') %>%
  dplyr::select(-tran) %>%
  dplyr::select(QX1410,IPR_accession,IPR_description,GO) %>%
  dplyr::filter(IPR_description != "-") %>%
  dplyr::distinct()

# write.table(ipr %>% dplyr::rename(`QX1410 transcript` = QX1410, `IPR accession` = IPR_accession, `IPR description` = IPR_description, `Gene ontology ID` = GO), "../../tables/supplementary/IPR_annotations_QX1410.tsv", quote = F, col.names = T, row.name = F, sep = "\t")

#==============================================================================================================================================================================================================================#

# INTERPROSCAN FUNCTIONAL ENRICHMENT

#==============================================================================================================================================================================================================================#
ipr_background_genes <- ipr %>% dplyr::distinct(QX1410) %>% dplyr::pull() # 15,289 / 22,196
  
# Pulling QX1410 genes that fall inside the QTL
qtl_genes <- allQX %>% dplyr::filter(chrom == "X" & start <= 16163640 & end >= 14531153) %>% 
  dplyr::distinct(QX1410) %>% dplyr::pull() # 311 genes in QTL 

# Define universe & QTL (annotated-only universe) 
univ_genes <- ipr_background_genes
qtl_genes_ipr  <- intersect(qtl_genes, ipr_background_genes) # 121 genes in the QTL have an IPR annotation

N <- length(univ_genes)
n <- length(qtl_genes_ipr)

# Counts per IPR (k = in universe, x = in QTL subset) 
k_tbl <- ipr %>%
  dplyr::count(IPR_accession, name = "k")

x_tbl <- ipr %>%
  dplyr::filter(QX1410 %in% qtl_genes) %>%
  dplyr::count(IPR_accession, name = "x")

desc_tbl <- ipr %>%
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

ipr_sig <- ipr_enrichment %>% # NONE PASS FDR-ADJUSTED P-VALUE CUTOFF!!!!!!
  dplyr::filter(pval < 0.05)

ipr_sig_gene_collapsed <- ipr %>%
  dplyr::filter(IPR_accession %in% ipr_sig$IPR_accession) %>%
  dplyr::group_by(IPR_accession) %>%
  dplyr::summarise(
    IPR_description = dplyr::first(stats::na.omit(IPR_description)),
    n_genes_QTL = dplyr::n_distinct(QX1410[QX1410 %in% qtl_genes]),
    genes_QTL   = paste(sort(unique(QX1410[QX1410 %in% qtl_genes])), collapse = ", "),
    n_genes_all = dplyr::n_distinct(QX1410),
    genes_all   = paste(sort(unique(QX1410)), collapse = ", "),
    .groups = "drop") %>%
  dplyr::left_join(ipr_sig %>% dplyr::select(n, N, `IPR accession` = IPR_accession, Expected = expected, `Fold enrichment` = enrich_ratio, `Odds ratio` = OR, `p-value` = pval, `FDR-adjusted p-value` = FDR_p.adjust), by = c("IPR_accession" = "IPR accession")) %>%
  dplyr::arrange(`p-value`) %>%
  dplyr::rename(`Number of genes in QTL with specific IPR annotation` = n_genes_QTL, `Genes in QTL with specific IPR annotation` = genes_QTL, 
                `Number of genes with specific IPR annotation in entire genome` = n_genes_all, `Genes in the genome with specific IPR annotation` = genes_all,
                `All genes in QTL with any IPR annotation` = n, `All genes in the entire genome with any IPR annotation` = N)

# Save a table
write.table(ipr_sig_gene_collapsed, "../../tables/supplementary/IPR_enrichment_QTL.tsv", row.names = F, col.names = T, sep = "\t")
