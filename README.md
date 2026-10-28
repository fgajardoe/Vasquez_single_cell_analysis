# Vasquez expression and colocalization analyses

Analyses for the research work by Daniel Vasquez et al.

## Overview

This repository contains two bulk-data analyses: differential expression of public microarray studies, and ChIP-seq colocalization.

## Analyses

### Differential expression (microarray)

`Differential_expression/run.R` reads precomputed differential expression tables and draws volcano plots (`logFC` vs `-log(adj.P.Val)`). Genes labeled when they pass an adjusted p-value cutoff of 0.05 are TET1, TET2, TET3, OGA, OGT, KMT2A–D, KMT2F, KMT2G, SETD1A, and SETD1B. Vertical guides are drawn at log2 fold-change ±2.

Studies used in `Differential_expression/run.sh`:

- Haferlach et al., 2010 — GSE13204, T-ALL vs healthy
- Homminga et al., 2011 — GSE26713
- Cramer-Morales et al., 2013 — GSE48558

`Differential_expression/merge.R` combines the three volcano plots into one panel.

`Differential_expression/boxplot_Dani/` uses GSE13204 (GPL570). Intensities are background-corrected (`normexp`) and normalized between arrays, then plotted as boxplots of TET1, KMT2A, and OGT for T-ALL versus non-leukemia and healthy bone marrow.

### ChIP-seq colocalization

`Colocalization_analysis/run.R` works with mouse embryonic stem cell ChIP-seq:

- Ogt, Vella et al., 2013 (GSE39154), filtered against IgG
- Tet1, Williams et al., 2011 (GSE24843), filtered against IgG

Peaks are restricted to chromosome 10, filtered by score, reduced, and counted against mm9 promoters (`TxDb.Mmusculus.UCSC.mm9.knownGene`).

## Citation

If you use this repository, please cite the work by Daniel Vasquez et al.

## Authors

Felipe Gajardo  
Daniel Vásquez
