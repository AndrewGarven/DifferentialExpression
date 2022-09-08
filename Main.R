source("TCGAdownload.R")
source("DESEQDifferentialExpression.R")

data <- TCGADownload('TCGA-BLCA','Transcriptome Profiling', 'RNA-Seq', 'STAR - Counts', 'Gene Expression Quantification')

rna <- as.data.frame(SummarizedExperiment::assay(data))
clinical <- data.frame(data@colData)

DESEQDifferentialExpression(rna, clinical, 'definition', 'Solid_Tissue_Normal', 'OutputFIle.csv')