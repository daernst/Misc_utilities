## Make plot of log2FC values for genes of interest from DESeq2 results

### Function:

#### Required parameter:
*data*: a dataframe containing columns with gene identifiers and corresponding log2FC values; column names must be "Gene" and "log2FC"


#### Optional parameter:
*ordered*: whether to order the plot by effect size (i.e., the absolute value of log2FC)

```
### Create log2FC plot function ###
log2FCplot <- function(data, ordered = NULL){
  require(dplyr)
  require(ggplot2)
  data %>%
    ggplot(aes(x = if(is.null(ordered)){Gene} else{reorder(Gene, abs(log2FC))},
               y = log2FC,
               color = log2FC)) +
    scale_colour_gradient2(low = "blue",
                           mid = "lightgrey",
                           high = "red",
                           midpoint = 0) +
    geom_point(size = 5) +
    geom_hline(yintercept = 0,
               linetype = "dashed") +
    geom_segment(aes(y = 0, 
                     yend = log2FC, 
                     xend = Gene)) +
    coord_flip() +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = "none") +
    xlab("Gene")
}
```

### Make log2FC plot:
```
log2FCplot(data, ordered = TRUE)
```
