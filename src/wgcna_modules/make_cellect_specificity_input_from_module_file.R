############### SYNOPSIS ###################
# Generate specificity input for CELLECT from WGCNA module output table


### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# =============================== SETUP ================================ #
# ======================================================================= #

library(tidyverse)
library(here)

source(here("src/lib/load_functions.R")) # load sc-genetics library

# ======================================================================= #
# =============================== CONSTANTS ================================ #
# ======================================================================= #

wd <- here("src/wgcna_modules")
setwd(wd)

# ======================================================================= #
# =============================== WRANGE DATA ========================== #
# ======================================================================= #

file.table <- here("out/wgcna/modules.mousebrain_bmicelltypes.rwgcna_table.csv.gz")
df.long <- read_csv(file.table)
df.long %>% dplyr::distinct(module) %>% nrow() # 571

### Export module 'geneset': add human genes + symbols [PNT]
df.geneset <- df.long 
df.geneset <- df.geneset %>% mutate(ensembl_gene_id_mouse=ensembl)
df.geneset <- df.geneset %>% mutate(ensembl_gene_id=mouse_to_human_ortholog_gene_mapping(gene_ids=ensembl, type_mouse_gene_ids="ensembl")) %>%
  filter(!is.na(ensembl_gene_id)) # keep only mapping genes
df.geneset <- hs_add_gene_symbol_from_ensembl_ids(df.geneset, colname_geneids_from="ensembl_gene_id", colname_geneids_to="gene_symbol")
# df.geneset %>% filter(!toupper(df.geneset$hgnc) == df.geneset$gene_symbol) # PNT mapping is slightly different from JT.
df.geneset <- df.geneset %>% select(-hgnc)
file.out <- here("out/wgcna/modules.mousebrain_bmicelltypes.rwgcna_table.human_only_geneset.csv.gz")
df.geneset %>% write_csv(file.out)

dfx <- df.geneset %>% dplyr::count(module)

## long --> matrix
df <- df.long %>% select(module, ensembl, pkME) %>% spread(key='module', value='pkME')
df[is.na(df)] <- 0 # ALT:  df %>% mutate_all(~replace(., is.na(.), 0)) 

### map mouse to human orthologs
df <- df %>% as.data.frame() %>% column_to_rownames("ensembl")
df.human <- mouse_to_human_ortholog_gene_expression_mapping(df, type_mouse_gene_ids="ensembl")
# [1] "Reading gene mapping files..."
# [1] "Subsetting on orthologs..."
# [1] "Number of genes not mapped to human ortholog: 729 out of 8912 genes (8.18 pct)"
df.human <- df.human %>% rownames_to_column(var="ensembl") %>% as_tibble()


### Export specificity matrix
file.out <- here("out/wgcna/modules.mousebrain_bmicelltypes.specificity_human.csv.gz")
df.human %>% write_csv(file.out)



### Make module 'metadata' file
# file.out.meta <- here("out/wgcna/modules.mousebrain_bmicelltypes.metadata.csv.gz") 
# use ---> results/modules--metadata.txt instead


