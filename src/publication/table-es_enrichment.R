############### SYNOPSIS ###################
### Combine and export ES enrichment (Mousebrain + Hypothalamus)

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

data_prefixes <- get_scrna_seq_dataset_prefixes("brain")

### Get ES data
dir.data <- here("out/es_enrichment_test")
files_all <- list.files(dir.data, pattern="*.csv.gz")
files_all
files_keep <- sapply(data_prefixes, function(prefix){files_all[base::startsWith(files_all, prefix)]}) # returns named char vec. names==data_prefixes
files_keep
# grepl(paste0("^", data_prefixes, collapse="|"), files_all) # ALT for extracting correct files
list.dfs <- lapply(file.path(dir.data, files_keep), read_csv)
names(list.dfs) <- names(files_keep)
df <- bind_rows(list.dfs, .id="specificity_id")

### Rename 'cell-type' to annotation
df <- df %>% rename(annotation = cell_type)

# ======================================================================= #
# ================================= PROCESS ============================= #
# ======================================================================= #

### add dataset
df <- df %>% mutate(dataset = case_when(
  specificity_id=="mousebrain"~"MSN",
  specificity_id %in% get_scrna_seq_dataset_prefixes("hypo") ~"Hypothalamus"))

### Load hypothalamus metadata with 'clean names'
df.metadata.hyp <- get_metadata("hypothalamus")

### Map hypothalamus annotations to 'clean names'
df <- df %>% mutate(annotation_fmt = df.metadata.hyp$annotation_fmt[match(annotation, df.metadata.hyp$annotation)]) # mousebrain annotation will get NA values
df <- df %>% mutate(annotation_fmt = if_else(is.na(annotation_fmt), annotation, annotation_fmt)) # replace NA with original value (cell_type)

# ======================================================================= #
# ================================= EXPORT ============================== #
# ======================================================================= #
  

# df.export <- df %>% select(Dataset=dataset, Annotation=annotation_fmt, `Enrichment P-value`=p.value)
df.export <- df %>% select(dataset, annotation, annotation_fmt, p.value)

file.out <- here("src/publication/tables/table-es_enrichment.combined.csv")
df.export %>% write_csv(file.out)

