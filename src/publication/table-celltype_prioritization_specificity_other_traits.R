############### SYNOPSIS ###################
### AIM: Check if MB prioritized cell-types are FDR significant in other traits

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
# ================================= RUN ================================= #
# ======================================================================= #


### PARAMS
filter.annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")
filter.gwas <- utils.get_gwas_ids_for_meta_analysis()
dataset_prefix <- "mousebrain"

### Load LDSC results
file.data <- here("results", sprintf("prioritization_celltypes--%s.multi_gwas.csv.gz", dataset_prefix))
df <- suppressMessages(read_csv(file.data))
df

# ================================= FILTER ================================= #
### ALL GWAS - P-values filtered
df.fdr <- df %>% 
  filter(annotation %in% filter.annotations) %>% 
  filter(gwas %in% filter.gwas) %>%
  filter(p.value<0.05/265) # *OBS*: hardcoded number of cell-types
df.fdr  


# ================================= EXPORT ================================= #
### prepare Rename GWAS
tmp_gwas_vector <- filter.gwas # convenience selector. If you want a specific order of the result, this vector should contain (unique) ordered values
newnames <- utils.rename_gwas(tmp_gwas_vector, style="fullname_author_year") # "fullname_author_year","fullname","abrv_author_year","abrv_year","abrv" 
rename_vector <- newnames; names(rename_vector) <- tmp_gwas_vector

### For each annotation, list FDR significant GWAS
df.fdr.by_annotation <- df.fdr %>% 
  mutate(gwas = recode_factor(gwas, !!!rename_vector)) %>%
  group_by(annotation) %>%
  summarise(gwas_fdr = paste0(gwas, collapse="\n"))

### Add metadata
df.metadata <- get_metadata(dataset_prefix)
df.fdr.by_annotation <- df.fdr.by_annotation %>% left_join(df.metadata, by="annotation")

### Write output file
file.out <- here("src/publication/tables/table-celltype_mb_prioritized_fdr_traits.csv")
df.fdr.by_annotation %>% write_csv(file.out)



# ======================================================================= #
# ================================= LEFTOVERS ================================= #
# ======================================================================= #



### ALL GWAS - P-values filtered | SPREAD
# df.fdr <- df %>% 
#   filter(annotation %in% filter.annotations) %>% 
#   select(gwas, p.value, annotation) %>% 
#   mutate(p.value = if_else(p.value<0.05/265, p.value, NA_real_)) %>%
#   spread(key=gwas, value=p.value) # format file
