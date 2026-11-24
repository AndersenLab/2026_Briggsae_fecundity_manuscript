#!/bin/bash

IN="../../processed_data/genomic_resources"
wkdir="../../processed_data/functional_enrichment"

./interproscan.sh \
	--formats TSV \
	--input $IN/c_briggsae.QX1410_20250929.csq.longest.nomito.prot.fa \
	--goterms \
	--cpu 48 \
	--iprlookup \
	--disable-precalc \
	--output-file-base $wkdir/QX_IPR_allApps_20251006 \
	--tempdir $TMP
