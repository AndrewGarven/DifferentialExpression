# The Cancer Genome Atlas (TCGA) RNA-seq Differential Expression

## Project Description
This project simplifies the process of pulling RNA-sequencing data from the open source data repository The Cancer Genome Atlas (TCGA). It then simplifies the process of performing both EdgeR and DESEQ2 differential expression. In developing a simplified workflow, this repo aims to encourage running multiple differential expression pipelines to target exaggerated false positives As shown by Yumei Li, et al. (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02648-4)

## future developments
1. Addition of Simple Wilcoxon rank-sum test (Proven best FDR control by Yumei Li, et al.) (September 2022) (ADDED: SEPT 15, 2022)
2. Graphical Visualization of all Differential Expression analysis (October 2022)
3. Comparative Statistics for each Differential Expression analysis (November 2022)

## Dependancies

### Differential Expression
1. <ins>TCGAbiolinks</ins> (Bioconductor Package) https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html
2. <ins>apeglm</ins> (Bioconductor Package) https://bioconductor.org/packages/release/bioc/html/apeglm.html
4. <ins>DESeq2</ins> (Bioconductor Package) https://bioconductor.org/packages/release/bioc/html/DESeq2.html
5. <ins>edgeR</ins> (Bioconductor Package) https://bioconductor.org/packages/release/bioc/html/edgeR.html
6. <ins>org.Hs.eg.db</ins> (Bioconductor Package) https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html

### Plotting
7. <ins>GGally</ins> (CRAN package) https://cran.r-project.org/web/packages/GGally/index.html
8. <ins>svglite</ins> (CRAN package) https://cran.r-project.org/web/packages/svglite/index.html

 
