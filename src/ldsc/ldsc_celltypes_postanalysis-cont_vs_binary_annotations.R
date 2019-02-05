############### SYNOPSIS ###################
# Assess robustness of LDSC CTS results for tabula muris and mousebrain
# - compare results binary and continuous annotations
# - compare individuals SEM results.

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# =============================== SETUP ================================= #
# ======================================================================= #

library(tidyverse)


dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib/"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/ldsc/"
setwd(wd)


# ======================================================================= #
# ================================ FUNCTIONS ================================= #
# ======================================================================= #

load_ldsc_cts_results <- function(file.ldsc_cts, dataset_prefix) {
  # dataset_prefix: character vector given the data set prefix. This is used to seperate the sem_name from the annotation_name.
  
  df.ldsc_cts <- read_tsv(file.ldsc_cts) # The last column gives a P-value from a one-sided test that the coefficient is greater than zero. 
  ### Name column examples in CTS file:
  # celltypes.mousebrain.all.binary.mousebrain_all.DEGLU5.tstat
  # celltypes.campbell_lvl1.all.campbell_lvl1.a18.Neurons6.tstat
  # celltypes.tabula_muris.all.binary.tabula_muris.Brain_Non-Myeloid.Bergmann_glial_cell.specificity
  pattern <- sprintf(".*\\.%s\\.(.*)\\.(.*?)$",dataset_prefix)
  mat.match <- str_match(df.ldsc_cts$Name, pattern) # returns matrix where the first column is the complete match, followed by one column for each capture group
  mat.match
  df.ldsc_cts <- df.ldsc_cts %>% mutate(sem=mat.match[,3],
                                         annotation=mat.match[,2],
                                         dataset=dataset_prefix) %>%
     select(-Name) # drop Name col now

  ### split string - ONLY works if there are no dots (.) in the annotation_name. (Works for mousebrain)
  # tmp_split <- stringr::str_split_fixed(df.ldsc_cts$Name,pattern="\\.", n=Inf)
  # df.ldsc_cts <- df.ldsc_cts %>% mutate(sem=tmp_split[,length(tmp_split[1,])],
  #                                       annotation=tmp_split[,length(tmp_split[1,])-1],
  #                                       dataset=tmp_split[,length(tmp_split[1,])-2]) %>% 
  #   select(-Name) # drop Name col now
  
  df.ldsc_cts <- df.ldsc_cts %>% rename(estimate=Coefficient, std.error=Coefficient_std_error, p.value=Coefficient_P_value)
  df.ldsc_cts <- df.ldsc_cts %>% mutate(fdr_significant = if_else(p.value <= 0.05/n(), true=T, false=F))
  # df.ldsc_cts <- df.ldsc_cts %>% left_join(df.metadata, by="annotation") # add meta data
  return(df.ldsc_cts)
}


# ======================================================================= #
# =============================== LOAD LDSC CTS RESULTS ================================= #
# ======================================================================= #

### *UPDATE*: binary result files live in out/out.ldsc/collection.celltypes.before_new_naming_scheme/

### Mousebrain
# dataset_prefix <- "mousebrain_all"
# files.ldsc_cts <- c("binary"="/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/celltypes.mousebrain.all.binary.BMI_Yengo2018.baseline_v1.1_all_genes.cell_type_results.txt",
#                      "continuous"="/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/celltypes.mousebrain.all.BMI_Yengo2018.baseline_v1.1_all_genes.cell_type_results.txt")

### Tabula Muris
dataset_prefix <- "tabula_muris"
files.ldsc_cts <- c("binary"="/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/celltypes.tabula_muris.all.binary.BMI_Yengo2018.baseline_v1.1_all_genes.cell_type_results.txt",
                    "continuous"="/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/celltypes.tabula_muris.all.BMI_Yengo2018.baseline_v1.1_all_genes.cell_type_results.txt")

df.ldsc_cts <- bind_rows(lapply(files.ldsc_cts, load_ldsc_cts_results, dataset_prefix), .id="annotation_type")
df.ldsc_cts


# ======================================================================= #
# ====================== COMPARE BINARY TO CONTINUOUS =================== #
# ======================================================================= #

df.ldsc_cts.bin_vs_cont <- df.ldsc_cts %>% 
  mutate(p.value.mlog10 = -log10(p.value)) %>%
  select(sem, p.value.mlog10, annotation, annotation_type) %>% # need to remove all columns with 'SEM unique' values before calling spread()
  spread(key=annotation_type, value=p.value.mlog10)

df.ldsc_cts.bin_vs_cont

ggplot(df.ldsc_cts.bin_vs_cont, aes(binary, continuous)) + geom_point() + geom_abline() + ggpubr::stat_cor(method = "pearson")

# ======================================================================= #
# ====================== COMPARE SEMS RESULTS =================== #
# ======================================================================= #

library(GGally)

df.ldsc_cts.spread <- df.ldsc_cts %>% 
  filter(annotation_type=="continuous") %>% 
  # filter(annotation_type=="binary") %>% 
  mutate(p.value.mlog10 = -log10(p.value)) %>%
  select(sem, p.value.mlog10, annotation) %>% # need to remove all columns with 'SEM unique' values before calling spread()
  spread(key=sem, value=p.value.mlog10)
  
df.ldsc_cts.spread

p <- GGally::ggpairs(df.ldsc_cts.spread, columns=c("ges", "sem_mean", "si", "specificity", "tstat"))
p

# ======================================================================= #
# =============================== XXXXX ================================= #
# ======================================================================= #

# ======================================================================= #
# =============================== XXXXX ================================= #
# ======================================================================= #

