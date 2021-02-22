## Make boxplots of normalized gene counts from DESeq2

### Function:

#### Required parameters:
*dds*: a DESeqDataSet

*gene*: a character, specifying the name of the gene to plot

*intgroup*: character vector of names in colData(x) to use for grouping

#### Optional parameter:
*order*: character vector specifying plotting order of *intgroup*

```
geneBoxplot <- function(dds, gene, intgroup, order = NULL){
  require(DESeq2)
  require(ggplot2)
  
  ### Extract counts data for gene ###
  gcount <- plotCounts(dds = dds, 
                       gene = gene, 
                       intgroup = intgroup, 
                       normalized = TRUE,
                       transform = TRUE,
                       returnData = TRUE, 
                       replaced = FALSE)
                       
  ### Set plotting order ###
  if(!is.null(order)){
    gcount[,intgroup] <- factor(gcount[,intgroup], 
                            levels = order)
    }
  ### Plot boxplots ###
  p <- ggplot(data = gcount,
               aes(x = gcount[,intgroup], 
                  y = count)) + 
    stat_boxplot(aes(x = gcount[,intgroup], 
                     y = count), 
                 geom='errorbar', 
                 linetype = 1, 
                 width = 0.25) +
    geom_boxplot(aes(fill = gcount[,intgroup]), 
                 outlier.shape = NA) +
    labs(x = NULL, 
         y = "Normalized counts", 
         title = gene) +
    scale_fill_brewer(palette = "Paired") +
    expand_limits(y = 0) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(), 
          legend.position = "none", 
          plot.title = element_text(hjust = 0.5, 
                                    face = "bold")) +
    geom_jitter(width = 0.1)
  return(p)
}

```

### Plot single gene:

```
### Make boxplot for single gene ###
geneBoxplot(dds = dds,
            gene = "BANY.1.2.g00001",
            intgroup = "stage",
            order = c("Larva", "Adult"))
```

### Plot multiple genes in multipanel:

```
### Load required package ###
library(gridExtra)

### Specify gene list ###
genes <- c("BANY.1.2.g00001", "BANY.1.2.g00003", "BANY.1.2.g00004")

### Plot boxplots as multipanel ###
do.call(grid.arrange, 
        c(lapply(genes, 
                 geneBoxplot, 
                 dds = dds, 
                 intgroup = "stage", 
                 order = c("Larva", "Adult")), 
          ncol = 3))
```
