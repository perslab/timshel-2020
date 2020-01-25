############### SYNOPSIS ###################
### AIM: Make table for ESmu genes for selected annotations

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #

library(tidyverse)
library(here)

source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))


# ======================================================================= #
# ================================ LOAD DATA ============================ #
# ======================================================================= #

# dataset_prefix <- "mousebrain"
dataset_prefix <- "hypothalamus"

filter.annotations <- get_prioritized_annotations_bmi(dataset_prefix)

### Get ES data
if (dataset_prefix != "hypothalamus") {
  df.es <- read_csv(here(sprintf("out/es/%s.mu.csv.gz", dataset_prefix)))
} else {
  df.es <- get_combined_hypo_es(merge="full")
}

### Add gene symbols
df.es <- hs_add_gene_symbol_from_ensembl_ids(df.es, colname_geneids_from="gene", colname_geneids_to="gene_name")
## Hypothalamus
# [1] "Number of genes mapped: 14594"
# [1] "Number of genes not mapped: 116"
# [1] "Number of genes with a NON-unique mapping (genes with duplicated ensembl gene IDs after mapping): 0"
# [1] "Total mapping stats: 116 genes have no mapping (not mapped + duplicates) out of 14710 input genes."
# [1] "Total genes mapped (non NA genes): 14594"
# [1] "Returning tibble with the column 'gene_name' added where all gene identifiers unique. Unmapped genes have NA values"

### ES table for all genes (SOM)
df.all_table <- df.es %>% select(gene, gene_name, !!filter.annotations)
file.out <- sprintf("tables/table-es_mu.%s.csv.gz", dataset_prefix)
df.all_table %>% write_csv(file.out)



# ======================================================================= #
# ============= Mousebrain - *WORKS*, uses s.es_obj.RData ================ #
# ======================================================================= #

# dataset_prefix <- "mousebrain"
# 
# ### Get ES data
# load(here(sprintf("out/es/%s.es_obj.RData", dataset_prefix)))
# 
# # filter.annotations <- es_obj[["annotations"]]
# filter.annotations <- get_prioritized_annotations_bmi("mousebrain")
# 
# ### ES table for all genes (SOM)
# df.all_table <- get_annotation_es.table(es_obj, annotations=filter.annotations, es_metric="es_mu") # ---> contains "gene" and "gene_name" column
# file.out <- sprintf("tables/table-es_mu.%s.csv.gz", dataset_prefix)
# df.all_table %>% write_csv(file.out)
# 
