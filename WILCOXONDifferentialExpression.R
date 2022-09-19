library('edgeR')
library('org.Hs.eg.db')

WILCOXONDifferentialExpression <- function(rna, clinical, clinicalFactor, referenceFactor, OutputFileName) {
  
  #INPUT
  # rna = data.frame object of transcriptional expression
  # clinical = data.frame object containing clinical information for subjects in 'rna'
  # clinicalFactor = string of 'clinical' column name used as factor for differential expression
  # referenceFactor = string of 'clinicalFactor' baseline comparator (ex. 'Healthy Tissue')
  # OutputFileName = string containing desired output file name ** must include '.csv'
  
  #Output
  # .csv file (located @ 'OutputFileName' contain per gene Wilcoxon differential expression analysis results
  
  if (clinicalFactor %in% colnames(clinical)){
    # obtain index of 'clinicalFactor'
    index = match(clinicalFactor,colnames(clinical))
    
    # assign 'clinicalFactor' to clinical$clinicalFactor and level into 2 groups
    clinical$clinicalFactor <- clinical[,index]
    clinical$clinicalFactor <-  gsub(" ", "_", clinical$clinicalFactor)
    clinical$clinicalFactor <- as.factor(clinical$clinicalFactor)
    clinical$clinicalFactor <- relevel(clinical$clinicalFactor, ref = referenceFactor)
    

    # Create a DGEList object with rna-counts and group (clinicalFactor)
    DGErna = DGEList(counts=rna, group=clinical$clinicalFactor)
    row.names(DGErna) = row.names(rna)
    
    # filter out lowly expressed genes
    keep <- filterByExpr(DGErna)
    DGErna <- DGErna[keep, , keep.lib.sizes=FALSE]
    
    # Library size normalization (scaling factor calculation)
    DGErna <- calcNormFactors(DGErna, method='TMM')
    rnaNORM<-as.data.frame(DGErna)
    
    row.names(rnaNORM) = row.names(DGErna)
    
    
    pvalues <- sapply(1:nrow(rnaNORM),function(i){
      data<-cbind.data.frame(gene=as.numeric(t(rnaNORM[i,])),clinical$clinicalFactor)
      p=wilcox.test(gene~clinical$clinicalFactor, data)$p.value
      return(p)
    })
    
    fdr=p.adjust(pvalues,method = "fdr")
    ref <- which(clinical$clinicalFactor == referenceFactor,arr.ind=T)
    nonref <- which(clinical$clinicalFactor != referenceFactor,arr.ind=T)
    
    dataCon1=rnaNORM[,row.names(clinical)[which(clinical$clinicalFactor==referenceFactor)]]
    dataCon2=rnaNORM[,row.names(clinical)[which(clinical$clinicalFactor!=referenceFactor)]]
    foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))
    
    rownames(rnaNORM) <- gsub("\\..*","",rownames(rnaNORM))
    
    outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr, symbol=(mapIds(org.Hs.eg.db, keys = rownames(rnaNORM), keytype = 'ENSEMBL', column = 'SYMBOL', multiVals = "first")))
    rownames(outRst)=rownames(rnaNORM)

    write.csv(outRst, file= OutputFileName)
    
  }
}