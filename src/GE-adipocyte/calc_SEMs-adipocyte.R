############### SYNOPSIS ###################
# Calculate adipocyte ES

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


dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib/"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library


setwd(here("src/GE-adipocyte"))


# ======================================================================= #
# ================================== HUMAN ============================== #
# ======================================================================= #

# dataset_prefix <- "preadipocyte_developing_1808_branch"
dataset_prefix <- "preadipocyte_developing_1808_branch_pc2_quantile"

sem_obj <- create_sem_object(sprintf("/projects/timshel/sc-scheele_lab_adipose_fluidigm_c1/data-preadipocytes_developing/%s", dataset_prefix))
sem_obj <- exclude_sporadic_expressed_genes(sem_obj)

### preadipocyte_developing_1808_branch
# [1] "Number of genes with NA in anova pvalue: 63. (This number is not used for anything. We just report it for your information)"
# [1] "Number of sporadic_expressed_genes: 11652"
# [1] "Number of zero-mean genes: 63"
# [1] "Total number of genes marked for exclusion: 11715"
# [1] "Number of genes in object: 11264"
### preadipocyte_developing_1808_branch_pc2_quantile
# [1] "Number of genes with NA in anova pvalue: 63. (This number is not used for anything. We just report it for your information)"
# [1] "Number of sporadic_expressed_genes: 11252"
# [1] "Number of zero-mean genes: 63"
# [1] "Total number of genes marked for exclusion: 11315"
# [1] "Number of genes in object: 11664"

sem_obj <- calc_sem_wrapper(sem_obj) # NEW
sem_obj <- calc_empirical_pvalues_wrapper(sem_obj) # NEW
sem_obj <- transform_sems(sem_obj, method="rank_normalize")
sem_obj <- set_group_by_annotation_slots(sem_obj)
sem_obj <- calc_sem_meta(sem_obj)

print("Saving object...")
# save(sem_obj, file=sprintf("%s.sem_obj.RData", dataset_prefix))

# ======================================================================= #
# ================== Write multi_geneset_file for LDSC ================= #
# ======================================================================= #

# load(file=sprintf("%s.sem_obj.RData", dataset_prefix))

### Write file
df_multi_geneset <- write_multi_geneset_file(sem_obj, 
                                             dataset_prefix=dataset_prefix,
                                             write_file=F)
df_multi_geneset

### ***** IMPORTANT *****
### OBS: the dataset contains a few duplicate gene IDs (because of Pytrik's 10x mapping). 
### make_annot_from_geneset_all_chr.py will fail if there are duplicates.
### Let's get rid of them quickly.
n_distinct(sem_obj$genes) # 11660
length(sem_obj$genes) # 11664
genes.duplicated <- sem_obj$genes[duplicated(sem_obj$genes)]
genes.duplicated

### Filter on sem_mean and write file
file.out <- sprintf("multi_geneset.%s.sem_mean.txt.gz", dataset_prefix)
df_multi_geneset.es_mean <- df_multi_geneset %>% 
  filter(.sem_name=="sem_mean") %>%
  filter(!gene %in% genes.duplicated) %>% # remove duplicated
  select(geneset_name, gene, value) %>% # *IMPORTANT*: re-arrange column order to agree with make_annot.py file format
  write_tsv(file.out, col_names=F) # no header

# ======================================================================= #
# ============================== LOAD =============================== #
# ======================================================================= #


# ======================================================================= #
# ======================== Hierarchical analysis ======================== #
# ======================================================================= #
