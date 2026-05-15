# mdRAD_v2
# RAD-like method for profiling DNA methylation

This RAD protocol generates reads initiated by methylated CpG sites (about half of all methylated sites in genome), and makes it possible to identfy exacly which CpG site was responsible. The data can be analyzed on per-region basis, by simply counting reads mapping to the region, or on per-CpG basis. The per-region analysis is more robust for initial data expoloration, while per-CpG analysis helps pinpoint the most differentilly methylated sites within a candidate region.

## Overview

Briefly, DNA is digested by two restriction enzymes: the methylation-sensitive MspJI recognizing mCNNR sites, and methylation-insensitive MseI recognizing TTAA sites. Then, adaptors are ligated allowing for selective PCR amplification of MspJI – MseI fragments. These fragments are then sequenced from the MspJI side. 

The lab protocol is `mdRAD_protocol_v2_mse.docx`, the oligos required (including some of the standard Illumina barcodes) are listed in `mdRAD_oligos_v2_mse.xlsx`.

## Analysis 

- trimming and deduplicaton: `trim_mdRAD.pl`
- mapping to reference genome: `mdRAD_mapping_dataExtraction.sh` (note: this script is for Lonestar6 system of Texas Advanced Computing Center; cannibalize the `bowtie2` command for your system. The script also contains the data extraction loop for CpG analysis; we also have a separate script for that)
- extraction of per-CpG counts: `extract_mdrad.sh` (this script was created by Christopher Peterson)
- alternatively, for per-gene analysis, count reads mapping to your regions (most commonly, genes) using `featureCount`, https://doi.org/10.1093/bioinformatics/btt656, and analyze like RNA-seq data.

