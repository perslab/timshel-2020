############### SYNOPSIS ###################
# Analyze and export LDSC results for multiple datasets (tabula muris, mousebrain, ...)


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

setwd(here("src/GE-adipocyte"))

# ======================================================================= #
# ============================ PARAMETERS =============================== #
# ======================================================================= #


### PARAMS
genomic_annotation_prefix <- "celltypes.preadipocyte_developing_1808_branch_pc2_quantile"
dataset_prefix <- "preadipocyte_developing_1808_branch_pc2_quantile"
# ======================================================================= #
# =============================== LOAD LDSC CTS RESULTS ================================= #
# ======================================================================= #

### Single GWAS
# file.ldsc_cts <- sprintf("/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/%s__BMI_UPDATE_Yengo2018.cell_type_results.txt", genomic_annotation_prefix)
# df.ldsc_cts <- load_ldsc_cts_results(file.ldsc_cts, dataset_prefix)
# df.ldsc_cts <- df.ldsc_cts %>% filter(sem=="sem_mean")


### Load - MULTI GWAS
dir.data <- "/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/"
filenames <- list.files(path=dir.data,  pattern=sprintf("%s.(.*).cell_type_results.txt", genomic_annotation_prefix))
filenames <- filenames[!grepl(pattern="__CONDITIONAL__", filenames, perl=T)] # exclude any conditional results
list.dfs <- lapply(file.path(dir.data, filenames), load_ldsc_cts_results, dataset_prefix)
names(list.dfs) <- stringr::str_match(filenames, pattern=sprintf("%s__(.*).cell_type_results.txt", genomic_annotation_prefix))[,2] # ALT: filenames
names(list.dfs)
df.ldsc_cts <- list.dfs %>% bind_rows(.id="gwas")
df.ldsc_cts <- df.ldsc_cts %>% filter(sem=="sem_mean")

### format file [*not* compatible with multi-GWAS plotting]


# ======================================================================= #
# =========================== EXPORT to results ========================= #
# ======================================================================= #

# ### Export (selected columns)
# file.out <- here("results", sprintf("primary-%s.multi_gwas.csv.gz", dataset_prefix))
# file.out
# df.ldsc_cts %>% select(-sem, -dataset, -n_obs_sem, -fdr_significant, -p.value.adj) %>% write_csv(file.out)
# 

# ======================================================================= #
# =============================== EXPORT Pval summary ================================= #
# ======================================================================= #

### Multi GWAS
file.out <- sprintf("out.190306.%s.sem_mean.multi_gwas.csv", dataset_prefix)
df.ldsc_cts.export %>% write_csv(file.out)


### Filter
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
file.out <- sprintf("out.190306.%s.sem_mean.multi_gwas.pval_adj.csv", dataset_prefix)
df.ldsc_cts.export.mod %>% write_csv(file.out)

# ======================================================================= #
# =============================== PLOT ================================= #
# ======================================================================= #

trait <- "WHRadjBMI_UKBB_Loh2018"
trait <- "WHR_Pulit2019"
trait <- "BMI_UKBB_Loh2018"
trait <- "LIPIDS_HDL_Teslovich2010"
df.ldsc_cts.export.mod$LIPIDS_HDL_Teslovich2010

df.plot <- df.ldsc_cts %>% filter(gwas==trait)
ggplot(df.plot, aes(x=fct_reorder(annotation,-log10(p.value)), y=-log10(p.value))) + 
  geom_col() + 
  geom_hline(yintercept = -log10(0.05/nrow(df.plot))) + 
  labs(x="annotation", title=trait) +
  coord_flip()
  # theme(axis.title.x = element_text(angle=70))


