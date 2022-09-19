source("TCGAdownload.R")
source("DESEQDifferentialExpression.R")
source("EDGERDifferentialExpression.R")
source('WILCOXONDifferentialExpression.R')
source('ParallelCoordinatePlot.R')

# import TCGA dataset
data <- TCGADownload('TCGA-BLCA','Transcriptome Profiling', 'RNA-Seq', 'STAR - Counts', 'Gene Expression Quantification')

# seperate dataset into RNA-counts and clinical information
rna <- as.data.frame(SummarizedExperiment::assay(data))
clinical <- data.frame(data@colData)

write.csv(rna, file='allRNA.csv')

# run DESeq2 differential expression - details in 'DESEQDifferentialExpression.R' results saved to 'OutputFileName'
DESEQDifferentialExpression(rna, clinical, 'definition', 'Solid_Tissue_Normal', 'OutputFIle.csv')

# run edgeR differential expression - details in 'EDGERDifferentialExpression.R' results saved to 'OutputFileName'
EDGERDifferentialExpression(rna, clinical, 'definition', 'Solid_Tissue_Normal', 'edgeR_DE.csv')

# run wilcoxon differential expression - details in 'WILCOXONDifferentialExpression.R' results saved to 'OutputFileName'
WILCOXONDifferentialExpression(rna, clinical, 'definition', 'Solid_Tissue_Normal', 'WilcoxonDE.csv')

# generate a parallel Coordinate Plot from top differentially expressed genes - details in 'ParallelCoordinatePlot.R' results saved to 'OutputFileName'
ParallelCoordinatePlot(rna, clinical, 'definition', 'Solid Tissue Normal', 19, topRNA, 'PCP.svg')
