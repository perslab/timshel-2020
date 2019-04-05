############### SYNOPSIS ###################
### AIM: filter modules based on pair-wise kme correlation => 0.6 and retain only module with the most LDSC significant p-value.

### Inspired by Kim 2018: "Among significantly enriched pathway-trait pairs we calculated a pairwise correlation
# for every pair of annotations and retained the more significant pathway for correlated pathways with r ≥ 0.5 (as in a previous study (Li et al., 2018)). 
# If correlated pathways were enriched for different traits, we retained both pathways."

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

library(corrr)
library(GGally)


source(here("src/lib/load_functions.R")) # load sc-genetics library

setwd(here("src/wgcna_modules"))

# ======================================================================= #
# ================================ PARAMETERS ================================ #
# ======================================================================= #

### LOAD kME
file.kme <- "/projects/jonatan/mousebrain_7/tables/Neurons_sub_ClusterName_7.2_run1_kMs_full_join.csv.gz" # this file correspons to RUN_ID="wgcna.mousebrain-190213.fdr_sign_celltypes.continuous"
df.kme <- read_csv(file.kme)

### LOAD LDSC CTS RESULTS
gwas <- "BMI_UKBB_Loh2018"
genomic_annotation_prefix <- "wgcna.mousebrain-190213.fdr_sign_celltypes.continuous"
file.ldsc_cts <- sprintf("/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/%s__%s.cell_type_results.txt", genomic_annotation_prefix, gwas)
df.ldsc_cts <- read_tsv(file.ldsc_cts)
df.ldsc_cts <- df.ldsc_cts %>% 
  mutate(module_id = str_split_fixed(Name, "__", n=Inf)[,2]) %>%
  select(module_id, p.value=Coefficient_P_value)
df.ldsc_cts
# ======================================================================= #
# ================================ Correlation ================================ #
# ======================================================================= #

### Compute correlation
df.corrr <- df.kme %>% select(-genes) %>% correlate(use = "pairwise.complete.obs", method="pearson")

### Convert to long format
df.corrr.long <- df.corrr %>%
  shave() %>% # # Convert the upper or lower triangle of a correlation data frame (cor_df) to missing values.
  stretch(na.rm = TRUE) # convert to long format. na.rm = TRUE: matrix diagonal (self-correlation) have NA values and will be dropped.
# x             y                      r
# <chr>         <chr>              <dbl>
#   1 antiquewhite3 antiquewhite3    NA     
# 2 antiquewhite3 aquamarine2      -0.375 
# 3 antiquewhite3 blue1            -0.134 


### Add ldsc data
df.corrr.long <- df.corrr.long %>% left_join(df.ldsc_cts, by=c("x"="module_id")) %>% rename(p.value.x = p.value)
df.corrr.long <- df.corrr.long %>% left_join(df.ldsc_cts, by=c("y"="module_id")) %>% rename(p.value.y = p.value)
df.corrr.long
df.corrr.long


### How many modules have a correlation r > threshold
df.corrr.long %>% filter(abs(r) > 0.7) %>% nrow() # 180
df.corrr.long %>% filter(abs(r) > 0.6) %>% nrow() # 453
df.corrr.long %>% filter(abs(r) > 0.5) %>% nrow() # 878


# ======================================================================= #
# ================================ Pruning/Filtering ================================ #
# ======================================================================= #

r_threshold <- 0.7
df.corrr.long.r_threshold <- df.corrr.long %>% filter(abs(r) >= r_threshold)

### Run pruning
## TODO: make into a function: prune_modules(df.corr.long, r_threshold)
## return list.prunned
## ALTERNATIVE ALGORITHM using square correlation matrix as input: https://stackoverflow.com/questions/29294983/how-to-calculate-correlation-between-all-columns-and-remove-highly-correlated-on

### DOCS
## This is a greedy algorithm, and may give different results for different input orders.
## E.g. it is possible that moduleA is clumped under moduleB, but could also be clumped under moduleC (and moduleB is not correlated enough to also be clumped under moduleC).
## The function garantuees that none of the 'lead' modules will have a pairwise correlation >r_threshold.
list.prunned <- list() # key=lead module; value=prunned modules
continue_iteration <- T
iter_counter <- 1
df.leftover <- as.data.frame(df.corrr.long.r_threshold) # convert to data frame to allow for df[i,j] operations
while (nrow(df.leftover)>0) {
  print(sprintf("iteration = %s | nrow = %s", iter_counter, nrow(df.leftover)))
  idx_win <- which.min(c(df.leftover[1, "p.value.x"], df.leftover[1, "p.value.y"])) # returns 1 if x is most significant and 2 if y is most significant
  if (idx_win == 1) { # x most sign.
    colname_win <- "x"
    colname_loose <- "y"
  } else { # y most sign.
    colname_win <- "y"
    colname_loose <- "x"
  }
  name_win <- df.leftover[1, colname_win]
  name_loose <- df.leftover[1, colname_loose]
  print(sprintf("%s won over %s and got [%s]", name_win, name_loose, paste(list.prunned[[name_loose]], collapse=",")))
  list.prunned[[name_win]] <- unique(c(name_loose, list.prunned[[name_loose]], list.prunned[[name_win]])) # update list. list.prunned[[<uninitialized_value>]] returns NULL
  # ^ we add the loose (name_loose) and any of the looser's 'wins' (list.prunned[[name_loose]]) and the winner's previous 'wins' (list.prunned[[name_win]])
  list.prunned[[name_loose]] <- NULL # delete the looser
  df.leftover <- subset(df.leftover, x!=name_loose & y!=name_loose) # remove looser from data frame
  iter_counter <- iter_counter + 1 
}
df.leftover
list.prunned
length(list.prunned)
list.prunned$lightpink3
# lightpink3
# lavenderblush

# ======================================================================= #
# ============= Check that none of lead modules are correlated ========= #
# ======================================================================= #

### Plot
# df.corrr %>% focus(names(list.prunned), mirror=TRUE) %>% rplot()
### Inspect
df.safety_check <- df.corrr %>% focus(names(list.prunned), mirror=TRUE) %>% shave() %>% stretch(na.rm = TRUE) %>% arrange(desc(abs(r)))
head(df.safety_check)
### Check
stopifnot(max(abs(df.safety_check$r)) <= r_threshold) # ok

# ======================================================================= #
# ================================ Finalize ================================ #
# ======================================================================= #

### Make df of lead modules
df.modules_pruned <- tibble(module_id = names(list.prunned), modules_pruned = sapply(list.prunned, paste, collapse=";"), n_modules_pruned = sapply(list.prunned, length))
### Get vector of modules that were prunned (i.e. they were not the lead module)
modules_pruned <- unlist(list.prunned, use.names=FALSE) # use.names=FALSE not needed | REF: http://r.789695.n4.nabble.com/How-to-union-the-elements-in-a-list-tp818444p818449.html
### Make a complete list of modules
df.modules_pruned <- df.ldsc_cts %>% filter(!module_id %in% modules_pruned) %>% left_join(df.modules_pruned, by="module_id") # add all modules

### Export
file.out <- sprintf("%s-%s-wgcna_modules_prunned-r_threshold_%s.csv", r_threshold, genomic_annotation_prefix, gwas)
file.out
df.modules_pruned %>% write_csv(file.out)


# ======================================================================= #
# ================================ caret findCorrelation() ================================ #
# ======================================================================= #

# library(caret)
# mat.cor = cor(df.kme %>% select(-genes), use="pairwise.complete.obs")
# hc = findCorrelation(mat.cor, cutoff=0.8) # putt any value as a "cutoff" 


# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #

# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #


