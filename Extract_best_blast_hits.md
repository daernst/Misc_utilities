## Function to extract best blast hit by evalue for each subject_id from full tabular report

### Function:

#### Required parameter:

*blast_file*: a .tsv file containing a full tabular blast report

#### Optional parameter:

*select_columns*: a vector containing the columns to extract from the *blast_file*

```
### Create best blast hit function ###
bb_hits <- function(blast_file, select_columns = NULL){
  require(readr)
  require(dplyr)
  
  ### Read in blast results and set column names ###
  blast_res <- read_tsv(blast_file, 
                        col_names = FALSE, 
                        comment = "#") %>%
    setNames(c("query_id", "query gi", "query acc.", "subject_id", "subject ids", 
               "subject gi", "subject gis", "subject_acc.", "subject accs.", "q. start", 
               "q. end", "s. start", "s. end", "query seq", "subject seq", 
               "evalue", "bit_score", "score", "alignment length", "perc_identity", 
               "identical", "mismatches", "positives", "gap opens", "gaps", 
               "% positives", "query/sbjct frames", "query frame", "BTOP", "subject tax ids", 
               "subject sci names", "subject com names", "subject blast names", "subject super kingdoms", "subject title", 
               "subject titles", "subject strand", "perc_query_coverage_per_subject", "% query coverage per hsp"))
  
  ### Extract best blast hit for each subject_id, breaking evalue ties by highest bit_score ###
  best_hits <- blast_res %>% 
    group_by(subject_id) %>%
    filter(evalue == min(evalue)) %>%
    slice(which.max(bit_score))
    
  ### Extract specified columns (optional) ###
  if(!is.null(select_columns)){
    best_hits <- best_hits %>% 
      select(select_columns)
  }
  best_hits <<- best_hits
}
```

### Extract best blast hits and desired columns:

```
### Extract best blast hits and columns with bb_hits function ###
bb_hits(blast_file = "Blastp_results.tsv",
        select_columns = c("query_id", "subject_id", "evalue", "bit_score"))
```
