# mdRAD_v2
## RAD-like method for profiling DNA methylation

This RAD protocol generates reads initiated by methylated CpG sites (about half of all methylated sites in genome), and makes it possible to identfy exacly which CpG site was responsible. The data can be analyzed on per-region basis, by simply counting reads mapping to the region, or on per-CpG basis. The per-region analysis is more robust for initial data expoloration, while per-CpG analysis helps pinpoint the most differentilly methylated sites within a candidate region.

# Overview

Briefly, DNA is digested by two restriction enzymes: the methylation-sensitive MspJI recognizing mCNNR sites, and methylation-insensitive MseI recognizing TTAA sites. Then, adaptors are ligated allowing for selective PCR amplification of MspJI – MseI fragments. These fragments are then sequenced from the MspJI side. 

# Analysis 

- trimming and deduplicaton: `trim_mdRAD.pl`
- 

