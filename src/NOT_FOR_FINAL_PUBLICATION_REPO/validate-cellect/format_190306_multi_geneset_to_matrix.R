############### SYNOPSIS ###################
# Format MB 190306 multigeneset into ESmu matrix to be able to run CELLECT and replicate previous results.


### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# =============================== SETUP ================================= #
# ======================================================================= #

library(tidyverse)
library(here)

# source(here("src/lib/load_functions.R")) # load sc-genetics library

# ======================================================================= #
# =============================== REFORMAT ================================= #
# ======================================================================= #

### Load data
file.multigeneset <- here("src/validate-cellect/data-190306/multi_geneset.mousebrain_all_190306_es_fix.txt.gz")
df.raw <- read_tsv(file.multigeneset, col_names=c("annotation_full_string", "gene", "value"))

### Wrangle
df <- df.raw %>% separate(col=annotation_full_string, into=c("prefix", "annotation", "es_metric"), sep="\\.")
df %>% count(es_metric) # all sem_mean. Good...
df.wide <- df %>% select(gene, annotation, value) %>% spread(key="annotation", value="value")
df.wide

df.wide[is.na(df.wide)] <- 0 # ALT:  df %>% mutate_all(~replace(., is.na(.), 0)) 
df.wide

### Export
file.out <- here("src/validate-cellect/data-190306/specificity_matrix.mousebrain_all_190306_es_fix.csv.gz")
df.wide %>% write_csv(file.out)

# ======================================================================= #
# =============================== XXXXX ================================= #
# ======================================================================= #



# ======================================================================= #
# =============================== XXXXX ================================= #
# ======================================================================= #

