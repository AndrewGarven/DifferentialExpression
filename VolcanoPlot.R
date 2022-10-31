library('EnhancedVolcano')

VolconoPlotAdvanced <- function(DEresult, GeneSymbol, adjustPCutoff, FCcutoff, OutputFileName) {
    
  #INPUT
  # DEresult = dataframe containing output from differential expresison analysis (Contains, gene name, fold-change, pval)
  # GeneSymbol = list object (subset from DEresult dataframe) // using current protocol DEresult$symbol contains gene symbols
  # adjustPCutoff = float value containing cutoff value for adjusted-pval (recommend 10e-6)
  # FCcutoff = float value containing a minimum log2 fold change value
  # OutputFileName = string containing desired output file name ** must include '.pdf'
  
  #OUTPUT 
  #Output
  # gg volcano plot 
  # .pdf of parallel coordinate plot
  
  x <- EnhancedVolcano(DEresult,
                  lab = GeneSymbol,
                  x = 'log2FoldChange',
                  y = 'padj',
                  xlab = bquote(~Log[2]~ 'fold change'),
                  pCutoff = adjustPCutoff,
                  FCcutoff = FCcutoff,
                  pointSize = 4.0,
                  labSize = 6.0,
                  colAlpha = 1,
                  legendPosition = 'right',
                  legendLabSize = 12,
                  legendIconSize = 4.0,
                  drawConnectors = TRUE,
                  widthConnectors = 0.75)
  x
  ggsave(OutputFileName)

}