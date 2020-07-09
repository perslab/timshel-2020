############### SYNOPSIS ###################
### Produce table of difference of campbelll 2017 lvl2 ESmu scores when
# running CELLEX on the whole dataset versus just neurons
# for high-confidence obesity genes across neurons

### OUTPUT: 
# celltype * gene table of delta ESmu with only neurons as baseline

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #

library(here)
library(tidyverse)

source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))

# ======================================================================= #
# ================================ LOAD DATA ============================ #
# ======================================================================= #


### Get ES data
file_es_campbell_lvl2_all <- here("out","es", "campbell2017_lvl2.mu.csv.gz")
file_es_campbell_lvl2_neur <- here("out","es","campbell2017_lvl2_neur.mu.csv.gz")
file_obesitygenes <- here("data","genes_obesity","prot_mono_early_extr_obesity_23genes.csv")

df.es_campbell_lvl2_all <- read_csv(file_es_campbell_lvl2_all)
df.es_campbell_lvl2_neur <- read_csv(file_es_campbell_lvl2_neur)
df.obesitygenes <- read_csv(file_obesitygenes)

df.markers = get_metadata("hypothalamus") %>% select(annotation_fmt, annotation_marker)

# ======================================================================= #
# ================================= PROCESS ============================= #
# ======================================================================= #


df.es_campbell_lvl2_all <- df.es_campbell_lvl2_all %>% 
  filter(gene %in% df.obesitygenes$ensembl_gene_ID) %>% 
  select(eval(colnames(df.es_campbell_lvl2_neur))) 
df.es_campbell_lvl2_neur <- df.es_campbell_lvl2_neur %>% 
  filter(gene %in% df.obesitygenes$ensembl_gene_ID) 

all.equal(colnames(df.es_campbell_lvl2_all), colnames(df.es_campbell_lvl2_neur))
# [1] TRUE

all.equal(df.es_campbell_lvl2_all$gene, df.es_campbell_lvl2_neur$gene)
#[1] "Lengths (18, 11) differ (string compare on first 11)" "7 string mismatches"

# add NA rows to neuron-only ESmu df
for (missing_gene in setdiff(df.es_campbell_lvl2_all$gene, df.es_campbell_lvl2_neur$gene)){
  df.es_campbell_lvl2_neur <- df.es_campbell_lvl2_neur %>% add_row(gene=missing_gene)
}

# reorder neuron-only ESmu df
df.es_campbell_lvl2_neur <- df.es_campbell_lvl2_neur[match(df.es_campbell_lvl2_all$gene, df.es_campbell_lvl2_neur$gene),]

all.equal(df.es_campbell_lvl2_all$gene, df.es_campbell_lvl2_neur$gene)
# [1] TRUE

# replace NAs with 0s
df.es_campbell_lvl2_neur[is.na(df.es_campbell_lvl2_neur)] <- 0
 
colnames(df.es_campbell_lvl2_all)[-1] <- 
  colnames(df.es_campbell_lvl2_neur)[-1] <- 
  utils.rename_annotations.hypothalamus(
    colnames(df.es_campbell_lvl2_all)[-1], 
    specificity_id = "campbell2017_lvl2", 
    check_all_matches=T)

df.es_difference = tibble("gene"=df.es_campbell_lvl2_all$gene,
                          df.es_campbell_lvl2_all[,-1]-df.es_campbell_lvl2_neur[,-1])

# ======================================================================= #
# ================================= EXPORT ============================== #
# ======================================================================= #

file.out <- here("src","publication","tables","table-es_cf_campbell_only_neurons.csv")
write_csv(df.es_difference, path=file.out)

