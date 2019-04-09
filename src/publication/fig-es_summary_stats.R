############### SYNOPSIS ###################
### AIM: generate ES summary stats (figures and tables)

### OUTPUT: 
# ....

### REMARKS:
# Most content is copied and modified from "src/expression_specificity/es_obj_summary_stats.R"

### REFERENCE:

# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #

library(tidyverse)
library(here)

source(here("src/lib/load_functions.R")) # load sc-genetics library

setwd(here("src/publication"))

# ======================================================================= #
# ============================ PARAMETERS =============================== #
# ======================================================================= #

### Tabula muris
# dataset_prefix <- "tabula_muris"
### Mousebrain
dataset_prefix <- "mousebrain_all"

# ======================================================================= #
# =========================== LOAD DATA ========================= #
# ======================================================================= #

# load("/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain.sem_obj.RData")
# sem_obj.mb <- sem_obj
# load("/projects/timshel/sc-genetics/sc-genetics/src/GE-maca/tabula_muris.sem_obj.RData")
# sem_obj.tm <- sem_obj


# ======================================================================= #
# =========================== Number of ES genes ========================= #
# ======================================================================= #

### Summary of number of ES genes
df.summary <- get_empirical_pvalues_summary(sem_obj, threshold=0.05, slot="sem_pvalues")
df.summary.stats <- sapply(df.summary %>% select(-annotation), summary) %>% # matrix with rownames from summary (Min., 1st Qu., Median, ..)
  as.data.frame() %>%
  rownames_to_column(var="statistic") %>%
  as.tibble()
df.summary.stats
df.summary.stats %>% select(-statistic) %>% summarise_all(.funs=funs(ratio = max(.)/min(.))) # min/max ratio

### Example output
# statistic tstat   ges    si specificity
# 1 Min.       643   471   137         390 
# 2 1st Qu.   1138   922.  526.        962.
# 3 Median    1692  1476   817        1517 
# 4 Mean      1940. 1687.  877.       1519.
# 5 3rd Qu.   2526. 2265  1228.       1892.
# 6 Max.      5178  4764  2010        3822 

### TODO: write out table

# =============================================================================================== #
# =============================== meta-SEM mean: counts per bin ================================= #
# =============================================================================================== #
# ***OBS*** this code is only valid for "binned" es objects

### meta-SEM mean: counts per bin
# df.sem_meta.bin_counts <- sem_obj.hier.Class[["Neurons"]][["sem_meta"]][["mean"]] %>% gather(key="annotation", value="sem_meta_bin") %>% group_by(annotation) %>% count(sem_meta_bin)
df.sem_meta.bin_counts <- sem_obj[["sem_meta"]][["mean"]] %>% gather(key="annotation", value="sem_meta_bin") %>% group_by(annotation) %>% count(sem_meta_bin)
df.sem_meta.bin_counts <- df.sem_meta.bin_counts %>% filter(!sem_meta_bin == 0) # remove zero-bin

# plot: number of non-zero genes per cell-type (with meta-data) [*this is a nice plot*]
df.sem_meta.bin_counts.summary <- df.sem_meta.bin_counts %>% 
  group_by(annotation) %>% 
  summarise(sum = sum(n)) %>% 
  left_join(df.metadata, by=c("annotation"))
df.sem_meta.bin_counts.summary %>% group_by(Class) %>% summarise(mean = mean(sum))

df.sem_meta.bin_counts.summary %>% 
  ggplot(aes(x=annotation, y=sum, fill=Class)) + geom_col() + geom_hline(yintercept = mean(df.sem_meta.bin_counts.summary$sum), color="red") + labs(y="number of non-zero bin genes per cell-type")


# ======================================================================= #
# =============== Empirical pval distribution [NOT IMPORTANT] =========== #
# ======================================================================= #

# ### Empirical pval distribution
# df.plot <- get_empirical_distribution(sem_obj, sem_name="ges", annotation=NULL) # ACBG (median=1)
# df.plot.gather <- df.plot %>% gather(key="distribution", value="sem_value")
# df.plot.gather %>% ggplot(aes(x=sem_value, fill=distribution)) + geom_density(alpha=0.2) + xlim(c(-10,10))
# df.plot.gather %>% ggplot(aes(x=sem_value, fill=distribution)) + geom_density(alpha=0.2) + xlim(c(10,1000))
