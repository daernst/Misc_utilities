## Function to make a volcano plot from DESeq2 results

### Function:

#### Required parameters:

*DESeq2_results*: a DESeqResults object

#### Optional parameters:

*highlight*: a vector containing the gene identifiers of genes to highlight in a different color (note: *highlight_color* must be specified)

*highlight_color*: the color of the points for the genes specified by *highlight*

- ***Note: genes with an FDR < 0.05 will be highlighted in red by default.***

```
### Create volcano plot function ###
volcano <- function(DESeq2_results, highlight = NULL, highlight_color = NULL){
  require(ggplot2)
  require(dplyr)
  require(tibble)
  
  ### Create dataframe with color column denoting differentially expressed genes (red =  FDR < 0.05; black = FDR > 0.05) ###
  dat <- as.data.frame(DESeq2_results) %>% 
    rownames_to_column() %>% 
    mutate(color = ifelse(DESeq2_results$padj < 0.05,
                          "red2",
                          "black"))
    
  ### Convert any padj=0 instances to smallest possible number for plotting (avoids Inf values; -log10(0)=Inf) ###
  dat$padj <- ifelse(dat$padj == 0,
                     .Machine$double.xmin,
                     dat$padj)
    
  ### Make volcano plot ###
  vol <- ggplot(dat,
                aes(x = log2FoldChange,
                    y = -log10(padj),
                    color = color)) +
    geom_point() + 
    scale_color_identity() +
    geom_hline(yintercept=0, 
               linetype="dashed") +
    geom_vline(xintercept=0, 
               linetype="dashed") +
    labs(x = expression("log"[2]*"FC"),
         y = expression("-log"[10]*"FDR")) +
    theme_bw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank()) +
    theme(legend.position="none") +
    theme(text = element_text(size=16))
  
  ### If specified, highlight specific genes ###
  if(!is.null(highlight) & !is.null(highlight_color)){
    gene_highlight <- dat %>%
      filter(rowname %in% highlight) %>%
      mutate(color = highlight_color)
    vol <- vol + 
      geom_point(gene_highlight, 
                 mapping = aes(x = log2FoldChange,
                               y = -log10(padj),
                               color=color))
  }
  vol <<- vol
  return(vol)
}
```

### Make volcano plot:

```
### Make vector of genes to highlight (optional) ###
intgenes <- c("BANY.1.2.g05189", "BANY.1.2.g10506", "BANY.1.2.g03920")


### Make volcano plot ###
volcano(DESeq2_results = res.lfcShrink, 
        highlight = intgenes, 
        highlight_color = "dodgerblue")
```
