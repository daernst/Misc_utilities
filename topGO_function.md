## Function to run topGO gene ontology enrichment analysis using DESeq2 results

### Function:

#### Required parameters:

*DESeq2_results*: a .csv file containing DESeq2 results (i.e., differentially expressed genes), with gene identifiers in the first column

*geneID2GO_file*: a .txt file containing a named list of character vectors, with gene identifiers as names, and character vectors containing the GO identifiers annotated to the gene (see topGO manual)

*category*: one of the three gene ontology categories; "BP", "MF", or "CC"

```
### Make function to run topGO GO enrichment analyses (Fisher's exact test, classic algorithm) ###
GO_enrich <- function(DESeq2_results, geneID2GO_file, category){
  require(topGO)
  require(dplyr)
  geneID2GO <- readMappings(file = geneID2GO_file) 
  geneNames <- names(geneID2GO)
  DEGs <- read.csv(DESeq2_results) %>% .[,1]
  geneList <- factor(as.integer(geneNames %in% DEGs))
  names(geneList) <- geneNames
  
  ### Build topGOdata object ###
  GOdata <<- new("topGOdata",
                 ontology = category,
                 allGenes = geneList,
                 annot = annFUN.gene2GO,
                 gene2GO = geneID2GO,
                 nodeSize = 5)
                 
  ### Run topGO ###
  resultFisher <- runTest(GOdata,
                          algorithm = "classic",
                          statistic = "fisher")
  allGO <- usedGO(object = GOdata)
  
  ### Get all results ###
  allRes <- GenTable(GOdata, 
                     p_value = resultFisher, 
                     orderBy = "resultFisher", 
                     ranksOf = "classic", 
                     topNodes = length(allGO))
  ### Correct p-values for multiple comparisons ###
  allRes$FDR <- p.adjust(allRes$p_value,
                         method = "fdr",
                         n = length(allRes$p_value))
  allRes <<- allRes
  sig_GOs <<- subset(allRes,
                     FDR <= 0.05)
}
```


### Run topGO for single GO category (i.e., BP, MF, or CC):

```
### Run topGO for biological process category ###
GO_enrich(DESeq2_results="DEGs_treatment_vs_control.csv",
          geneID2GO_file="geneID2GO.txt",
          category="BP")

```
