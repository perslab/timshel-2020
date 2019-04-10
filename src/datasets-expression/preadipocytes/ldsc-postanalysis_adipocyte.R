############### SYNOPSIS ###################
# Export LDSC results for adipocyte


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

source(here("src/lib/load_functions.R")) # load sc-genetics library

setwd(here("src/datasets-expression/preadipocytes"))

# ======================================================================= #
# ============================ PARAMETERS =============================== #
# ======================================================================= #


### PARAMS
genomic_annotation_prefix <- "celltypes.preadipocyte_developing_1808_branch_pc2_quantile"
dataset_prefix <- "preadipocyte_developing_1808_branch_pc2_quantile"

# ======================================================================= #
# =============================== LOAD LDSC CTS RESULTS ================================= #
# ======================================================================= #

### Load - MULTI GWAS
dir.data <- "/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/"
filenames <- list.files(path=dir.data,  pattern=sprintf("%s.(.*).cell_type_results.txt", genomic_annotation_prefix))
filenames <- filenames[!grepl(pattern="__CONDITIONAL__", filenames, perl=T)] # exclude any conditional results
list.dfs <- lapply(file.path(dir.data, filenames), load_ldsc_cts_results, dataset_prefix)
names(list.dfs) <- stringr::str_match(filenames, pattern=sprintf("%s__(.*).cell_type_results.txt", genomic_annotation_prefix))[,2] # ALT: filenames
names(list.dfs)
df.ldsc_cts <- list.dfs %>% bind_rows(.id="gwas")
df.ldsc_cts <- df.ldsc_cts %>% filter(sem=="sem_mean")

### format file [*not* compatible with multi-GWAS plotting]
df.ldsc_cts.export <- df.ldsc_cts %>% select(gwas, p.value, annotation) %>% spread(key=gwas, value=p.value) # format file

# ======================================================================= #
# =========================== EXPORT to results ========================= #
# ======================================================================= #

# ### Export (selected columns)
# file.out <- here("results", sprintf("primary-%s.multi_gwas.csv.gz", dataset_prefix))
# file.out
# df.ldsc_cts %>% select(-sem, -dataset, -n_obs_sem, -fdr_significant, -p.value.adj) %>% write_csv(file.out)
# 

# ======================================================================= #
# =========================== EXPORT Pval summary ======================== #
# ======================================================================= #

### ALL GWAS
file.out <- sprintf("out/genetic_prioritization.1808_branch_deciles.multi_gwas.csv")
df.ldsc_cts.export %>% write_csv(file.out)


### ALL GWAS - P-values filtered
df.ldsc_cts.export.mod <- df.ldsc_cts %>% 
  select(gwas, p.value, annotation) %>% 
  mutate(p.value = if_else(p.value<0.05/30, p.value, NA_real_)) %>%
  spread(key=gwas, value=p.value) # format file

df.ldsc_cts.export.mod <- df.ldsc_cts.export.mod %>% select(annotation, 
                                                            contains("WHR"), 
                                                            contains("LIPIDS"), 
                                                            contains("BMI"), 
                                                            contains("CARDIOVASCULAR"), 
                                                            contains("T2D"), 
                                                            everything())

# write
file.out <- sprintf("out/genetic_prioritization.1808_branch_deciles.multi_gwas.pval_adj.csv", dataset_prefix)
df.ldsc_cts.export.mod %>% write_csv(file.out)

