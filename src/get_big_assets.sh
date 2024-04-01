#!/bin/bash

# Get GENCODE v43 annotations
wget -qO- https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz \
    | gunzip -c - \
    | gtf2bed - \
    | grep -w gene - \
    > gencode.v43.genes.bed

# Get Human Genome version hg38
wget -q https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# Get RMVar data
wget -q https://rmvar.renlab.org/download/RMVar_Human_Clinvar_info.txt
wget -q https://rmvar.renlab.org/download/RMVar_Human_basic_info_m6A.txt

# Get latest ClinVar VCF file
wget -q https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
