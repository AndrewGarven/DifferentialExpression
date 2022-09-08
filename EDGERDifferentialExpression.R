library('edgeR')
library('biomaRt')
library('org.Hs.eg.db')

EDGERDifferentialExpression <- function(rna, clinical, clinicalFactor, referenceFactor, OutputFileName) {
  #INPUT

  # checking that selected 'clinicalFactor' exists in clinical data.frame
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
    
    # filter out lowly expressed genes
    keep <- filterByExpr(DGErna)
    DGErna <- DGErna[keep, , keep.lib.sizes=FALSE]
    
    # Library size normalization (scaling factor calculation)
    DGErna <- calcNormFactors(DGErna)
    
    # model design 
    design <- model.matrix(~clinical$clinicalFactor)
    
    # estimating dispersion using cox-reid profile-adjusted likelihood
    DGErna <- estimateDisp(DGErna, design)
    
    # fit negative-binomial generalized linear models
    fit <- glmQLFit(DGErna, design)
    
    # perform DE and store results
    results <- glmQLFTest(fit)
    
    # add symbol (gene Symbol) row to results table 
    rownames(results$table) <- gsub("\\..*","",rownames(results$table))
    results$table$symbol <- (mapIds(org.Hs.eg.db, keys = rownames(results$table), keytype = 'ENSEMBL', column = 'SYMBOL', multiVals = "first"))
    
    # write resultant CSV file to directory 'OutputFileName'
    write.csv(results$table, 
             file= OutputFileName)
  }
}
