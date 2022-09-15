## The Cancer Genome Atlas (TCGA) RNA-seq Differential Expression

# Project Description
This project simplifies the process of pulling RNA-sequencing data from the open source data repository The Cancer Genome Atlas (TCGA). It then simplifies the process of performing both EdgeR and DESEQ2 differential expression. In developing a simplified workflow, this repo aims to encourage running multiple different differential expression pipelines to target exaggerated false positives As shown by Yumei Li, et al. (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02648-4)

# future developments
1. Addition of Simple Wilcoxon rank-sum test (Proven best FDR control by Yumei Li, et al.) (September 2022)
2. Graphical Visualization of all Differential Expression analysis (October 2022)
3. Comparative Statistics for each Differential Expression analysis (November 2022)

# Dependancies
TCGAbiolinks (Bioconductor Package) https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html
apeglm (Bioconductor Package) https://bioconductor.org/packages/release/bioc/html/apeglm.html
biomaRt (Bioconductor Package) https://bioconductor.org/packages/release/bioc/html/biomaRt.html
DESeq2 (Bioconductor Package) https://bioconductor.org/packages/release/bioc/html/DESeq2.html
edgeR (Bioconductor Package) https://bioconductor.org/packages/release/bioc/html/edgeR.html
org.Hs.eg.db (Bioconductor Package) https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html
