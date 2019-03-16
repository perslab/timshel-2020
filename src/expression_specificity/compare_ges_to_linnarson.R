############### SYNOPSIS ###################
### AIM: Compare our GES calculations to Linnarson published ES values

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
# ================================ LOAD ES ================================= #
# ======================================================================= #

### load data
load(file="mousebrain.sem_obj.RData") # human

# ======================================================================= #
# ================= Compare to Linnarson original GES values ============ #
# ======================================================================= #
### Conclusion
# 1: My calculation works. There is a high correlation between my GES and Linnarson GES. The difference is due to normalization
# 2: My GES values are distributed 'close' to normal but with a long tail. The distribution is 'nicer' (less 'spiky', less variance) compared to Linnarson GES.
# 3: The median of my GES for each cell-type seems to be VERY close to 1. This is a nice attribute.

df.ges.linnarson <- read_file_fast("/raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/mousebrain/mousebrain-agg_L5.enrichment_values.csv.gz")
df.ges.my <- df.ges %>% mutate(gene=sem.ex[["genes"]])
df.ges.my.join <- df.ges.my %>% left_join(df.ges.linnarson, by="gene")
ggplot(df.ges.my.join, aes(x=ENT8.x, y=ENT8.y)) + geom_point()
with(df.ges.my.join, cor.test(x=ENT9.x, y=ENT9.y, method="spearman"))

ggplot(df.ges.my.join, aes(x=ENT8.x)) + geom_density() # my GES are more 'gently' (less spike) distributed
ggplot(df.ges.my.join, aes(x=ENT8.y)) + geom_density() + xlim(0,100)
summary(df.ges.my.join %>% pull(ENT8.x)) # ---> all cell-types have median GES VERY close to 1. (my GES)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.09231  0.95519  0.99979  1.16738  1.23329 17.68845 
summary(df.ges.my.join %>% pull(ENT8.y))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000    0.951    2.185    7.261    6.102 6017.554

df <- df.ges.my.join %>% summarise_if(is.numeric, .funs=funs(mean, median, min, max))
df.x <- as.data.frame(df %>% t()) %>% rownames_to_column(var="stat")


