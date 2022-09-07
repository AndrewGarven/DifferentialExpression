BiocManager::install("TCGAbiolinks")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
BiocManager::install("biomaRt")

library('apeglm')
library('edgeR')
library("TCGAbiolinks")
library("DESeq2")
library("biomaRt")

TCGADownload <- function(project, dataCategory, experimentalStratergy, WorkflowType, dataType) {
  # INPUT
  # project = string of TCGA project type (ex. 'TCGA-BLCA')
  # dataCategory = string of TCGA data.category (ex. 'Transcriptome Profiling')
  # expreimentalStratergy = string of TCGA experimental.stratergy (ex. 'RNA-Seq')
  # WorkflowType = string of TCGA workflow.type (ex 'STAR - Counts')
  # dataType = string of TCGA data.type (ex. 'Gene Expression Quantification')
  
  # OUTPUT
  # large data object containing both expression and clinical information
  
  # search using GDC API for given parameters
  TCGAquery = GDCquery(
    project = project,
    data.category = dataCategory,
    experimental.strategy = experimentalStratergy,
    workflow.type = WorkflowType,
    data.type = dataType)
  
  # download data elucidated in search query 
  GDCdownload(query = TCGAquery)
  
  # prepare downloaded data into R object
  data <- GDCprepare(query = TCGAquery, save = FALSE)
  
  return(data)
}

DESEQDifferentialExpression <- function(rna, clinical, clinicalFactor, referenceFactor, OutputFileName) {
  #INPUT
  # rna = data.frame object of transcriptional expression
  # clinical = data.frame object containing clinical information for subjects in 'rna'
  # clinicalFactor = string of 'clinical' column name used as factor for differential expression
  # referenceFactor = string of 'clinicalFactor' baseline comparator (ex. 'Healthy Tissue')
  
  if (clinicalFactor %in% colnames(clinical)) {
    index = match(clinicalFactor,colnames(clinical))
    clinical$clinicalFactor <- clinical[,index]
    clinical$clinicalFactor <-  gsub(" ", "_", clinical$clinicalFactor)
    clinical$clinicalFactor <- as.factor(clinical$clinicalFactor)
    clinical$clinicalFactor <- relevel(clinical$clinicalFactor, ref = referenceFactor)
    dds <- DESeqDataSetFromMatrix(countData = rna,
                                  colData = clinical,
                                  design = ~ clinicalFactor)
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    dds <- DESeq(dds)
    
    res <- results(dds, alpha = 0.05,  altHypothesis = "greaterAbs", lfcThreshold = 1.5) # alpha controls FDR rate
    
    #result Log fold change shrinkage method (suitable for logfc based visualization)
    resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
    resLFC.Ordered<-resLFC[with(resLFC, order(abs(log2FoldChange), padj, decreasing = TRUE)), ]
    
    # converting Ensebl id to Gene symbole using biomart
    ens2symbol<-function(ids){
      mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
      genes <- getBM(filters= "ensembl_gene_id", 
                     attributes= c("ensembl_gene_id","hgnc_symbol"),
                     values=ids, mart= mart)
      return(genes)
    }
    
    print('something')
    df <- ens2symbol(row.names(res))
    
    res_df <- as.data.frame(res)                 
    res_df$ensembl_gene_id <- row.names(res_df)
    res_df <- merge(df,res_df, by = "ensembl_gene_id")
    resOrdered<-res_df[with(res_df, order(abs(log2FoldChange), padj, decreasing = TRUE)), ]
    
    #saving the results
    write.csv(res_df, 
              file= 'crazy.csv')
  }
}

data <- TCGADownload('TCGA-BLCA','Transcriptome Profiling', 'RNA-Seq', 'STAR - Counts', 'Gene Expression Quantification')

rna <- as.data.frame(SummarizedExperiment::assay(data))
clinical <- data.frame(data@colData)

DESEQDifferentialExpression(rna, clinical, 'definition', 'Solid_Tissue_Normal', 'OutputFIle.csv')


