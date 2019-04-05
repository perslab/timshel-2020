############### SYNOPSIS ###################
# FIT SEM models
# Export Mousebrain enrichment scores



### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

library(tidyverse)
library(here)

source(here("src/lib/load_functions.R")) # load sc-genetics library

wd <- "/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain"
setwd(wd)


# ======================================================================= #
# ==============================  LOAD DATA  =============================== #
# ======================================================================= #

### We use CPM and log-transformed averaged data (as in the 'standard' workflow)
### File has mouse ensembl geneIDs
file.avg_expr <- "/projects/timshel/sc-genetics/sc-genetics/data/expression/mousebrain/mousebrain-L5.cpm_average_log.ensembl_musculus.csv.gz"
df.avg_expr <- read_csv(file.avg_expr) %>% column_to_rownames(var="X1")

### Ortholog mapping
df.human <- mouse_to_human_ortholog_gene_expression_mapping(df.avg_expr, type_mouse_gene_ids="ensembl")
# "Output dimension of final data table (ortholog mapped): n_genes=16125 x n_features=265"

# ======================================================================= #
# ==============================  SEMs  =============================== #
# ======================================================================= #

wrapper_calculate_and_write_sems(df.human, name.dataset="mousebrain")
# "OBS: N=417 genes removed because they contained NA values after zscore calculation. "

# ======================================================================= #
# ==============================  ***Enrichment scores**  =============================== #
# ======================================================================= #

### Read data
file.data <- "/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain-agg_L5.enrichment_values.csv.gz" # mouse gene symbols
df <- read_csv(file.data) %>% rename(gene=X1)
# ^^ n duplicate genes ~= 65
df <- df %>% distinct(gene, .keep_all = TRUE) # remove duplicate gene names
df <- df %>% column_to_rownames(var="gene") # set rownames


### Write file: enrichment scores
write_csv(df %>% rownames_to_column(var="gene"), path="mousebrain.celltype_expr.enrichment.csv.gz") # mouse genes
df.human <- mouse_to_human_ortholog_gene_expression_mapping(df)
write_csv(df.human %>% rownames_to_column(var="gene"), path="mousebrain.celltype_expr.enrichment.hsapiens_orthologs.csv.gz")  # orthologs

### Write log10(x+1) transformed ortholog
df.human.log10 <- log10(df.human+1)
class(df.human.log10) # --> data frame
any(is.nan(as.matrix(df.human.log10))) # ---> FALSE, no nan values
min(df.human.log10) > 0 # ---> TRUE, all values are positive (we want to know this for POS_TRANSFORMATION)
write_csv(df.human.log10 %>% rownames_to_column(var="gene"), path="mousebrain.celltype_expr.enrichment_log.hsapiens_orthologs.csv.gz")


### Write QUANTILES transformed ortholog
df.human.rank <- df.human %>% transmute_all(funs(transform_to_quantiles), threshold_bin_zero=1) # *OBS* output tibble looses rownames
# ^ OBS: threshold_bin_zero=1 because of distribution of Linnarson enrichments
head(df.human.rank)
rownames(df.human.rank) <- rownames(df.human) # set rownmaes *IMPORTANT*
write_csv(df.human.rank %>% rownames_to_column(var="gene"), path="mousebrain.celltype_expr.enrichment_quantiles.hsapiens_orthologs.csv.gz")

