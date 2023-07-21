#!/bin/bash

wget -qO- https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz \
    | gunzip -c - \
    | gtf2bed - \
    | grep -w gene - \
    > gencode.v43.genes.bed
