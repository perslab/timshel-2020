############### SYNOPSIS ###################
### AIM: Summary stats for es objects
# -Empirical pval distribution
# meta-SEM mean: counts per bin
# SCRIPT INCLUDE code for using the 'skimr' library

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


library(skimr) # REF: https://cran.r-project.org/web/packages/skimr/vignettes/Using_skimr.html


source(here("src/lib/load_functions.R")) # load sc-genetics library

setwd(here("src/expression_specificity"))


# ======================================================================= #
# ============================ PARAMETERS =============================== #
# ======================================================================= #

### Tabula muris
# dataset_prefix <- "tabula_muris"

### Mousebrain
# dataset_prefix <- "mousebrain_all"


# ======================================================================= #
# =========================== LOAD DATA ========================= #
# ======================================================================= #

# load("/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain.sem_obj.RData")
# sem_obj.mb <- sem_obj
# load("/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-maca/tabula_muris.sem_obj.RData")
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

# ======================================================================= #
# ====================== Empirical pval distribution ===================== #
# ======================================================================= #

### Empirical pval distribution
df.plot <- get_empirical_distribution(sem_obj, sem_name="ges", annotation=NULL) # ACBG (median=1)
df.plot.gather <- df.plot %>% gather(key="distribution", value="sem_value")
df.plot.gather %>% ggplot(aes(x=sem_value, fill=distribution)) + geom_density(alpha=0.2) + xlim(c(-10,10))
df.plot.gather %>% ggplot(aes(x=sem_value, fill=distribution)) + geom_density(alpha=0.2) + xlim(c(10,1000))

### qvals summary
df.summary <- get_empirical_pvalues_summary(sem_obj, threshold=0.05, slot="sem_qvalues")




# =============================================================================================== #
# =============================== OLD CODE NOT UPDATED YEY ================================= #
# =============================================================================================== #


# =============================================================================================== #
# =============================== meta-SEM mean: counts per bin ================================= #
# =============================================================================================== #
# ***OBS*** this code is only valid for "binned" es objects

### meta-SEM mean: counts per bin
# df.sem_meta.bin_counts <- sem_obj.hier.Class[["Neurons"]][["sem_meta"]][["mean"]] %>% gather(key="annotation", value="sem_meta_bin") %>% group_by(annotation) %>% count(sem_meta_bin)
df.sem_meta.bin_counts <- sem_obj[["sem_meta"]][["mean"]] %>% gather(key="annotation", value="sem_meta_bin") %>% group_by(annotation) %>% count(sem_meta_bin)
df.sem_meta.bin_counts <- df.sem_meta.bin_counts %>% filter(!sem_meta_bin == 0) # remove zero-bin
# plot: single-cell type
ggplot(df.sem_meta.bin_counts %>% filter(annotation == "DEINH3"), aes(x=sem_meta_bin, y=n)) + geom_col() # single cell-type
ggplot(df.sem_meta.bin_counts %>% filter(annotation == "MEINH2"), aes(x=sem_meta_bin, y=n)) + geom_col() # single cell-type
# plot: number of non-zero genes per cell-type (with meta-data)
df.sem_meta.bin_counts.summary <- df.sem_meta.bin_counts %>% 
  group_by(annotation) %>% 
  summarise(sum = sum(n)) %>% 
  left_join(df.metadata, by=c("annotation"))
df.sem_meta.bin_counts.summary %>% group_by(Class) %>% summarise(mean = mean(sum))
df.sem_meta.bin_counts.summary %>% 
  ggplot(aes(x=annotation, y=sum, fill=Class)) + geom_col() + geom_hline(yintercept = mean(df.sem_meta.bin_counts.summary$sum), color="red") + labs(y="number of non-zero bin genes per cell-type")

# plot: mean number of genes per bin
df.sem_meta.bin_counts %>% group_by(sem_meta_bin) %>% summarise(mean = mean(n)) %>% ggplot(aes(x=sem_meta_bin, y=mean)) + geom_col() + labs(y="mean number of genes in bin across all cell-types")


# =============================================================================================== #
# =============================== skimr: CV, SD, p100 ================================= #
# =============================================================================================== #
# ***TODO***: port this code to new workflow

DIR.X <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/mousebrain"
filenames <- list.files(path=DIR.X,  pattern="*.hsapiens_orthologs.csv.gz")
list.x <- lapply(file.path(DIR.X, filenames), read_csv)
list.x <- lapply(list.x, function(df) df %>% select(-gene))
names(list.x) <- filenames
list.sum <- lapply(list.x, skim) # summarise USING ***skim()***
list.sum
list.sum.sum <- lapply(list.sum, function(df) {
  df %>% select(-formatted) %>% 
    tidyr::spread(key="stat", value="value") %>%
    select(variable,mean,sd,p100) %>%
    mutate(cv=sd/mean)
})
list.sum.sum[[1]]

df.sum.sum <- bind_rows(list.sum.sum, .id = "file") # combine
df.sum.sum

### Inspect
df.x <- df.sum.sum %>% group_by(file) %>%
  top_n(5, p100)
df.x

### Plot
ggplot(df.sum.sum, aes(x=cv, fill=file)) + geom_density(alpha=0.4) + xlim(0,15) # CV
ggplot(df.sum.sum, aes(x=file, y=log10(cv), fill=file)) + geom_boxplot() # CV
ggplot(df.sum.sum, aes(x=file, y=log(p100), fill=file)) + geom_boxplot() # p100


############## TMP ##################
df <- list.sum[[1]]
head(df)
df %>% 
  select(-formatted) %>% 
  tidyr::spread(key="stat", value="value") %>%
  select(variable,mean,sd,p100) %>%
  mutate(cv=sd/mean)
# group_by(variable) %>%
# summarise(mean_mean=mean())
?spread

skim(list.x[[1]])
skim(list.x[[2]])
skim(list.x[[3]])

