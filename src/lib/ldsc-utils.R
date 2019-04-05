### Helper functions for LDSC post analysis

library(tidyverse)

# ======================================================================= #
# ============================== MISC =============================== #
# ======================================================================= #

get_genomic_annotation_prefix <- function(dataset_prefix) {
  if (dataset_prefix == "campbell_lvl1") {
    genomic_annotation_prefix <- "celltypes.campbell_lvl1.all"
  } else if (dataset_prefix == "campbell_lvl2") {
    genomic_annotation_prefix <- "celltypes.campbell_lvl2.all"
  } else if (dataset_prefix == "mousebrain_all") {
    genomic_annotation_prefix <- "celltypes.mousebrain.all"
  } else if (dataset_prefix == "tabula_muris") {
    genomic_annotation_prefix <- "celltypes.tabula_muris.all"
  } else {
    stop("wrong dataset_prefix")
  }
  return(genomic_annotation_prefix)
}


# ======================================================================= #
# ============================== LOAD CTS FILE =============================== #
# ======================================================================= #

load_ldsc_cts_results <- function(file.ldsc_cts, dataset_prefix) {
  # dataset_prefix: character vector given the data set prefix. This is used to seperate the sem_name from the annotation_name.
  
  df.ldsc_cts <- read_tsv(file.ldsc_cts) # The last column gives a P-value from a one-sided test that the coefficient is greater than zero. 
  ### Name column examples in CTS file:
  # celltypes.mousebrain.all.binary__mousebrain_all.DEGLU5.tstat
  # celltypes.campbell_lvl1.all__campbell_lvl1.a18.Neurons6.tstat
  # celltypes.tabula_muris.all.binary__tabula_muris.Brain_Non-Myeloid.Bergmann_glial_cell.specificity
  pattern <- sprintf(".*__%s\\.(.*)\\.(.*?)$",dataset_prefix)
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
  df.ldsc_cts.tmp.summary <- df.ldsc_cts %>% group_by(sem) %>% summarise(n_obs_sem=n())
  df.ldsc_cts <- left_join(df.ldsc_cts, df.ldsc_cts.tmp.summary, by="sem")
  df.ldsc_cts <- df.ldsc_cts %>% mutate(fdr_significant = if_else(p.value <= 0.05/n_obs_sem, true=T, false=F),
                                        p.value.adj = p.value*n_obs_sem)
  # df.ldsc_cts <- df.ldsc_cts %>% left_join(df.metadata, by="annotation") # add meta data
  return(df.ldsc_cts)
}

