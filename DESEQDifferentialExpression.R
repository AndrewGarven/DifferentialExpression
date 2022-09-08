library('apeglm')
library("DESeq2")
library("biomaRt")

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