library('GGally')
library('svglite')


ParallelCoordinatePlot <- function(rna, clinical, clinicalFactor, referenceFactor, CasePerFactor, TopDiffExp, OutputFileName) {
  
  #INPUT
  # rna = data.frame object of transcriptional expression
  # clinical = data.frame object containing clinical information for subjects in 'rna'
  # clinicalFactor = string of 'clinical' column name used as factor for differential expression
  # referenceFactor = string of 'clinicalFactor' baseline comparator (ex. 'Solid Tissue Normal')
  # CasePerFactor = int object outlining the number of cases (patients, samples, etc.) you would like to visualize per group (per clinicalFactor)
  # TopDiffExp = list of strings containing the top differentially expressed Genes for visualization
  # OutputFileName = string containing desired output file name ** must include '.svg'
  
  #Output
  # gg parallel coordinate plot 
  # .pdf of parallel coordinate plot
  
  index = match(clinicalFactor,colnames(clinical))
  clinical$clinicalFactor <- clinical[,index]

  positiveID <- row.names(clinical[clinical$clinicalFactor != referenceFactor, ])
  negativeID <- row.names(clinical[clinical$clinicalFactor == referenceFactor, ])
  
  positiveRNA <- rna[ ,sample(ncol(rna[ ,positiveID]), CasePerFactor)]
  negativeRNA <- rna[ ,sample(ncol(rna[ ,negativeID]), CasePerFactor)]
  
  allRNAselectedCases <- data.frame(positiveRNA, negativeRNA)

  SelectedRNASelectedCases <- allRNAselectedCases[TopDiffExp, ]
  
  plot <- ggparcoord(
    data = SelectedRNASelectedCases,
  )
  
  plot
  
  ggsave(OutputFileName, width = 8, height = 8, units = "cm")
}
