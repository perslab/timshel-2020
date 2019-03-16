############### SYNOPSIS ###################
### AIM: investigate the difference between a 'global' and 'hierarcial' ES

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

dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib/"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library

setwd(here("src/expression_specificity"))


# ======================================================================= #
# ========================== LOAD ES HIERARCHICAL ======================== #
# ======================================================================= #

### load data
load(file="mousebrain.sem_obj.hier.Class.RData") # sem_obj.hier.Class

# ======================================================================= #
# ============================== INVESTIGATE HIERARCHICAL ================================ #
# ======================================================================= #

# selected_annotation <- "SCINH3"
selected_annotation <- "DEINH3"
df.hier <- sem_obj.hier.Class[["Neurons"]][["sem_meta"]][["mean"]] %>% mutate(gene=sem_obj.hier.Class[["Neurons"]][["genes"]]) %>% select(gene, hier=!!sym(selected_annotation))
df.global <- sem_obj[["sem_meta"]][["mean"]] %>% mutate(gene=sem_obj[["genes"]]) %>% select(gene, global=!!sym(selected_annotation))
print(nrow(df.hier)); print(nrow(df.global)) # check that they have similar number of genes.
df.cmp <- inner_join(df.hier, df.global, by="gene") # *OBS* inner join.
df.cmp
df.cmp <- df.cmp %>% mutate(
  cases = case_when( # Like an if statement, the arguments are evaluated in order, so you must proceed from the most specific to the most general.
    (hier == 0) & (global == 0) ~ "both.zero",
    (hier > 0) & (global) ~ "both.non_zero",
    (hier > 0) ~ "hier.non_zero",
    (global > 0) ~ "global.non_zero"
  ))
df.cmp %>% count(cases)
# DEINH3
# 1 both.non_zero    1884
# 2 both.zero       10753
# 3 global.non_zero  1886
# 4 hier.non_zero     522

library("ggpubr")

ggplot(df.cmp, aes(x=hier, y=global)) + geom_point() + 
  geom_rug(col=rgb(.5,0,0,alpha=.2)) + 
  geom_smooth(method="lm") +
  geom_abline() +
  ggpubr::stat_cor(method = "pearson") +
  # geom_density_2d() + # does not work: Computation failed in `stat_density2d()`: bandwidths must be strictly positive
  labs(title=selected_annotation)

# ======================================================================= #
# ================================ GSEA ================================= #
# ======================================================================= #


do_gprofiler <- function(input, background){
  df.res = tryCatch({ # to avoid gprofiler error "transfer closed with outstanding read data remaining"
    df.res <- gProfileR::gprofiler(input, organism="hsapiens", ordered_query=F, significant=T, custom_bg=background)
    df.res <- df.res %>% mutate(call_status="Successful")
    return(df.res)
  }, warning = function(w) {
    message(sprintf("warning: %s", w))
  }, error = function(e) {
    message(sprintf("error: %s", e))
    return(data.frame(call_status="Failed")) # return single-row data frame
  }, finally = {
    message("done with function do_gprofiler")
  })
  return(df.res)
}

### GO test
# TODO: add function
df.cmp %>% filter(cases == "hier.non_zero") %>% pull(gene) %>% cat() # non neuronal specific
df.cmp %>% filter(cases == "global.non_zero") %>% pull(gene) %>% cat() # synaptic signaling, synapse
x <- df.cmp %>% filter(cases == "both.non_zero") %>% pull(gene) # synaptic signaling, synapse

df.go <- do_gprofiler(input=x, background=sem_obj[["genes"]])



# ======================================================================= #
# ======================== Hierarchical analysis ======================== #
# ======================================================================= #

# TaxonomyRank2              n
# <chr>                  <int>
#   1 CNS neurons              181
# 2 PNS neurons               33
# 3 CNS glia                  23
# 4 Neural crest-like glia    13
# 5 Vascular cells            10
# 6 Immune cells               5

# TaxonomyRank1      n
# <chr>          <int>
#   1 Neurons          214
# 2 Glia              36
# 3 Vascular cells    10
# 4 Immune cells       5

# Class              n
# <chr>          <int>
#   1 Neurons          214
# 2 Vascular          11
# 3 Astrocytes        10
# 4 Oligos            10
# 5 PeripheralGlia    10
# 6 Ependymal          5
# 7 Immune             5


