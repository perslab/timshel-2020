############### SYNOPSIS ###################
### Generate dataset summary tables: N cells per cell-type + N ES genes

### OUTPUT: 
# ....

### REMARKS:


### REFERENCE:


# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #

library(tidyverse)
library(here)

library(patchwork)

source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))

# ======================================================================= #
# ================================ PARAMS ============================== #
# ======================================================================= #

# dataset_prefix <- "tabula_muris"
# dataset_prefix <- "mousebrain"
# dataset_prefix <- "campbell2017_lvl2"

data_prefixes <- get_scrna_seq_dataset_prefixes("all")
## *SOURCE BELOW FUNCTION BEFORE RUNNING THIS:*
map(data_prefixes, wrapper_n_es_genes_summary_stats)


wrapper_n_es_genes_summary_stats <- function(dataset_prefix) {
  # ======================================================================= #
  # ================================ LOAD DATA ============================ #
  # ======================================================================= #
  message(dataset_prefix)
  load(here(sprintf("out/es/%s.es_obj.RData", dataset_prefix)))
  
  
  # ======================================================================= #
  # =========================== Number of ES genes ========================= #
  # ======================================================================= #
  
  ### Summary of number of ES genes for each metric
  df.n_es <- get_empirical_pvalues_summary(es_obj, threshold=0.05, slot="es_pvalues")
  # annotation tstat   ges    si specificity
  # 1 ABC         2701  2477  1265        2351
  # 2 ACBG        1502  1056  2126        2170
  # 3 ACMB        2939  2126  2169        2588
  # 4 ACNT1       1721  1343  2414        3364
  
  ### Calculate number of ES mu genes
  df.es_mu_genes <- es_obj[["es_meta"]][["mu"]] %>% 
    gather(key="annotation", value="es_weight") %>% 
    group_by(annotation) %>% 
    summarise(es_mu = sum(es_weight>0)) %>%
    ungroup()
  
  ### Join
  df.n_es <- df.n_es %>% left_join(df.es_mu_genes, by="annotation")
  df.n_es
  
  ### Rename
  # vec.rename <- c(DET="tstat", GES="ges", EP="specificity", NSI="si", ESmu="es_mu")
  vec.rename <- c(ESmu="es_mu")
  df.n_es <- df.n_es %>% rename(!!vec.rename)
  
  ### Get number of cells
  df.n_cells <- es_obj[["ncells"]] # 1 x n_annotations
  df.n_cells <- df.n_cells %>% gather(key="annotation", value="N_cells_es")
  df.n_es <- df.n_es %>% left_join(df.n_cells, by="annotation")
  
  # ======================================================================= #
  # ============================ SUMMARIZE ============================= #
  # ======================================================================= #
  
  df.n_es.sum_stats <- sapply(df.n_es %>% select(-annotation), summary) %>% # matrix with rownames from summary (Min., 1st Qu., Median, ..)
    as.data.frame() %>%
    rownames_to_column(var="statistic") %>%
    as.tibble()
  df.n_es.sum_stats
  # df.n_es.sum_stats %>% select(-statistic) %>% summarise_all(.funs=list(ratio ~max(.)/min(.))) # min/max ratio
  
  ### Example output
  # statistic tstat   ges    si specificity es_mu
  # 1 Min.       567   320   127         412  1329 
  # 2 1st Qu.   2276  1497   440         803  2955 
  # 3 Median    3393  2702   672        1078  4020 
  # 4 Mean      3452. 3055.  963.       1348. 4119.
  # 5 3rd Qu.   4455  4436  1295        1780  5174 
  # 6 Max.      7160  7625  2796        3761  8147 
  
  # ======================================================================= #
  # ============================ EXPORT TABLE(S) ============================= #
  # ======================================================================= #
  
  ### Write file
  file.out <- sprintf("tables/table-annotation_summary.%s.csv", dataset_prefix)
  df.n_es %>% write_csv(file.out)
  
  # file.out <- sprintf("tables/table-n_es_genes.%s.sum_stats.csv", dataset_prefix)
  # df.n_es.sum_stats %>% write_csv(file.out)
  message(sprintf("Done with %s", dataset_prefix))
}

# ======================================================================= #
# =============== Empirical pval distribution [NOT IMPORTANT] =========== #
# ======================================================================= #

# ### Empirical pval distribution
# df.plot <- get_empirical_distribution(es_obj, es_name="ges", annotation=NULL) # ACBG (median=1)
# df.plot.gather <- df.plot %>% gather(key="distribution", value="es_value")
# df.plot.gather %>% ggplot(aes(x=es_value, fill=distribution)) + geom_density(alpha=0.2) + xlim(c(-10,10))
# df.plot.gather %>% ggplot(aes(x=es_value, fill=distribution)) + geom_density(alpha=0.2) + xlim(c(10,1000))
