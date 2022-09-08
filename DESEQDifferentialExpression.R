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
    
    # add symbol (gene Symbol) row to results table
    rownames(resLFC.Ordered) <- gsub("\\..*","",rownames(resLFC.Ordered))
    resLFC.Ordered$symbol <- (mapIds(org.Hs.eg.db, keys = rownames(resLFC.Ordered), keytype = 'ENSEMBL', column = 'SYMBOL', multiVals = "first"))
    
    #save results to 'OutputFileName'
    write.csv(resLFC.Ordered, 
              file= OutputFileName)
  }
}