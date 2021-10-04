## Make dot plot of log2FC values for genes of interest from DESeq2 results and cluster by gene expression

### Function:

#### Required parameters:
*results*: a dataframe containing DESeq2 results; must contain columns named "Gene" and "log2FoldChange"

*vsd*: a DESeqTransform object containing variance stabilized counts from DESeq2

*goi*: a vector containing the identifiers for genes of interest to be plotted

```
### Create log2FC dot plot function ###
log2FCdendro <- function(results, vsd, goi){
  library(ggtree)
  library(aplot)
  library(SummarizedExperiment)
  library(cowplot)
  library(tidyverse)
  ### Cluster genes by expression ###
  mat <- as.data.frame(SummarizedExperiment::assay(vsd)) %>%
    dplyr::filter(row.names(.) %in% goi)
  clust <- hclust(dist(mat))
  ### Reorder genes of interest levels so that the log2FC dot plot gene order matches the dendrogram branch order ###
  branch_order <- ggplot2::fortify(clust) %>%
    dplyr::filter(isTip == TRUE) %>%
    dplyr::arrange(desc(y)) %>%
    dplyr::pull(label)
  results$Gene <- factor(results$Gene,
                         levels = branch_order) %>%
    fct_rev()
  ### Make log2FC dot plot ###
  l2fc_plot <- results %>%
    dplyr::filter(Gene %in% goi) %>%
    ggplot2::ggplot(aes(x = Gene,
                        y = log2FoldChange,
                        color = log2FoldChange)) +
    scale_colour_gradient2(low = "blue",
                           mid = "lightgrey",
                           high = "red",
                           midpoint = 0) +
    geom_point(size = 5) +
    geom_segment(aes(y = 0,
                     yend = log2FoldChange,
                     xend = Gene),
                 size = 1) +
    coord_flip() +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = "none") +
    labs(x = "Gene",
         y = expression("log"[2]*"FC")) +
    geom_hline(yintercept = 0,
               linetype = "dashed") +
    xlab(NULL) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0), "cm")) +
    scale_x_discrete(position = "top")
  ### Make gene expression dendrogram ###
  ggtree_plot <- ggtree::ggtree(clust) +
    aplot::ylim2(l2fc_plot) +
    theme(plot.margin = unit(c(0.5, -0.15, 0.5, 0.5), "cm"))
  ### Combine dendrogram with log2FC dot plot ###
  cowplot::plot_grid(ggtree_plot,
                     l2fc_plot,
                     ncol = 2,
                     align = 'h',
                     rel_widths = c(0.4, 1))
}
```

### Make log2FC dot plot:
```
### Make log2FC dot plot ###
log2FCdendro(results = results,
             vsd = vsd,
             goi = goi)
```
