# Processed Data

This directory contains processed datasets, GWAS results, and genomic resources used in the *Caenorhabditis briggsae* fecundity manuscript.

## Directory Structure

```text
processed_data/
│
├── 262strains_avgBS_A3567_4_1_gwasmapping.tsv
├── 311strains_avgBS_A3567_4_1_repwise.tsv
├── 311strains_avgBS_A3567_4_1_.tsv
│
├── day1_262strains_GWASmapping.tsv
├── day2_262strains_GWASmapping.tsv
├── day3_262strains_GWASmapping.tsv
├── day4_262strains_GWASmapping.tsv
├── day5_262strains_GWASmapping.tsv
│

│
├── functional_enrichment/
│   └── QX1410_IPR_annotations.tsv
│
├── genomic_resources/
│   ├── c_briggsae.QX1410_20250929.csq.gff
│   └── c_briggsae.QX1410_20250929.csq.longest.nomito.prot.fa
│
├── nemascan/
│   ├── 20260903_1/
│   ├── 20261008_262strains_day1/
│   ├── 20261008_262strains_day2/
│   ├── 20261008_262strains_day3/
│   ├── 20261008_262strains_day4/
│   ├── 20261008_262strains_day5/
│   └── 20261115_262strains_repday_1_2/
│
└── README.md
```

### Processed  Datasets

These files contain processed fecundity measurements and datasets used for GWAS mapping and related analyses.

- `311strains_avgBS_A3567_4_1_.tsv`
  - Processed average fecundity dataset for 311 strains.

- `311strains_avgBS_A3567_4_1_repwise.tsv`
  - Replicate-wise fecundity measurements for 311 strains.

- `262strains_avgBS_A3567_4_1_gwasmapping.tsv`
  - Processed average fecundity dataset for  262 strains used in  GWAS mapping.

### Daily GWAS Mapping Datasets

These files contain the datasets used for GWAS mapping of fecundity measurements across the five reproductive days in the 262-strain panel.

- `day1_262strains_GWASmapping.tsv`
- `day2_262strains_GWASmapping.tsv`
- `day3_262strains_GWASmapping.tsv`
- `day4_262strains_GWASmapping.tsv`
- `day5_262strains_GWASmapping.tsv`

### NemaScan GWAS Results

The `nemascan/` directory contains GWAS results generated using NemaScan for fecundity-associated traits in the 262-strain mapping panel.

Each analysis directory contains relevant GWAS outputs and supporting files, including genotype matrices, heritability estimates, independent-test information, strain quality-control information, software version information, and NemaScan reports.

### Analysis directories

- `20260903_1/`
  - GWAS analysis of average fecundity of 262 strains.

- `20261008_262strains_day1/`
  - GWAS analysis of day 1 fecundity of 262 strains.

- `20261008_262strains_day2/`
  - GWAS analysis of day 2 fecundity of 262 strains.

- `20261008_262strains_day3/`
  - GWAS analysis of day 3 fecundity of 262 strains.

- `20261008_262strains_day4/`
  - GWAS analysis of day 4 fecundity of 262 strains.

- `20261008_262strains_day5/`
  - GWAS analysis of day 5 fecundity of 262 strains.

- `20261115_262strains_repday_1_2/`
  - GWAS analysis of peak reproductive day of 262 strains.

### Genomic Resources

The `genomic_resources/` directory contains genomic annotation and protein sequence resources for the *C. briggsae* QX1410 reference genome.

- `c_briggsae.QX1410_20250929.csq.gff`
  - GFF file containing genomic feature and consequence annotation information.

- `c_briggsae.QX1410_20250929.csq.longest.nomito.prot.fa`
  - FASTA file containing the longest protein isoform per gene, excluding mitochondrial proteins.

### Functional Enrichment

The `functional_enrichment/` directory contains protein functional annotation data used for functional characterization and enrichment analyses.

- `QX1410_IPR_annotations.tsv`
  - InterPro annotation results containing protein domain and functional annotation information for the QX1410 reference genome.

