############### SYNOPSIS ###################
### AIM:
# Compare ES modifications (made March 2019: fixed tstat, GES, nSI) to their implementation before modifications (first implementation)
# March 2019 ES modifications: 
# - tstat: added var.pooled_scaling_factor to tstat_sem() function.
# - GES: fixed bug: FROM  (frac.other+epsilon_2)/(frac.other+epsilon_2) TO (frac.self+epsilon_2)/(frac.other+epsilon_2)
# - nSI: subtracting 1 in the both denominator and numerator of the final calculation to remove arbitrary 'bimodality'.

### OK TO DELETE THIS SCRIPT LATER. IT DOES NOT CONTAIN IMPORTANT CODE.

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

setwd(here("src/expression_specificity"))


# ======================================================================= #
# ================================ LOAD ES ================================ #
# ======================================================================= #

file.rdata.es <- here("src/GE-maca/tabula_muris.sem_obj-190306.RData") # newest calculations
load(file.rdata.es)
sem_obj.new <- sem_obj
rm(sem_obj)

file.rdata.es <- here("src/GE-maca/tabula_muris.sem_obj.RData") # first calculations
load(file.rdata.es)
sem_obj.old <- sem_obj
rm(sem_obj)

# ======================================================================= #
# ================================ COMPARE ================================ #
# ======================================================================= #

### SWTICH
sem_obj <- sem_obj.new
sem_obj <- sem_obj.old

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

### NEW
# statistic tstat   ges    si specificity
# 1 Min.       643   402   137         390 
# 2 1st Qu.   1138   865   526.        962.
# 3 Median    1692  1521   817        1517 
# 4 Mean      1940. *2089.*  877.       1519.
# 5 3rd Qu.   2526. 2968  1228.       1892.
# 6 Max.      5178  6685  2010        3822

### OLD
# statistic tstat   ges    si specificity
# <chr>     <dbl> <dbl> <dbl>       <dbl>
# 1 Min.       643   471   137         390 
# 2 1st Qu.   1138   922.  526.        962.
# 3 Median    1692  1476   817        1517 
# 4 Mean      1940. *1687.*  877.       1519.
# 5 3rd Qu.   2526. 2265  1228.       1892.
# 6 Max.      5178  4764  2010        3822 

# ======================================================================= #
# ============================ PLOTS SCATTER ========================= #
# ======================================================================= #

df.es.old <- get_es_annotation_centric(sem_obj.old, annotation=annotation)
df.es.new <- get_es_annotation_centric(sem_obj.new, annotation=annotation)
df.es.new

## join
df.es.join <- inner_join(df.es.old,df.es.new, by=c("gene", "es_metric"), suffix=c(".old", ".new"))
df.es.join

### plot
plot_df <-df.es.join %>% group_by(es_metric) %>%
  do(
    plots = ggplot(data = .) + aes(x = es_weight.old, y = es_weight.new) +
      geom_point() + ggtitle(.$es_metric)
  )

# plot_df$plots # show plots

library(gridExtra)

grid.arrange(grobs = plot_df$plots, ncol = 2) ## display plot

# ======================================================================= #
# ============================ PLOTS ========================= #
# ======================================================================= #

# ============================ ANNO+GENE PLOT ========================= #
annotation <- "Pancreas.pancreatic A cell" # GCG
gene_highlight <- "GCG"
df.es <- get_es_annotation_centric(sem_obj, annotation)

p <- es_plot_annotation_centric(df.es, genes_highlight=gene_highlight)
p + labs(title="Pancreatic Alpha-cells (glucagon+)")

df.es.gene <- get_mean_nes_gene_centric(sem_obj, genes_select=gene_highlight)
es_plot_gene_centric(df.es.gene, annotations_highlight=annotation)


# ============================ ANNO DENSITY PLOT ========================= #
df.es <- get_es_annotation_centric(sem_obj, annotation=annotation)
df.es <- df.es %>% mutate(
  es_weight_pseudo_log=sign(es_weight)*log10(abs(es_weight+1)),
  es_weight_log=log10(es_weight)
)
# df.es <- df.es %>% filter(es_metric=="tstat")
# df.es <- df.es %>% filter(es_metric == "mean nES")
df.es <- df.es %>% filter(es_weight == 0)
p <- ggplot(df.es, aes(x=es_weight)) + 
  geom_density() + 
  scale_x_log10() +
  facet_wrap(~es_metric, scales = "free")
p

df.es %>% count(es_metric)

## NEW
# es_metric       n
# <chr>       <int>
#   1 ges         12756
# 2 mean nES     4896
# 3 si          12756
# 4 specificity 12756
# 5 tstat        5639

### OLD
# es_metric       n
# <chr>       <int>
#   1 ges         12756
# 2 mean nES     3400
# 3 si          14303
# 4 specificity 12756
# 5 tstat        5639

# ---> NEW more ES genes > 0
# why did si have 14303 genes in the old? it makes sense because >0. OK DELETE


