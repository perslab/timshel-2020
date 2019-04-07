

library(tidyverse)


wd <- "/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/"
setwd(wd)

# =========================================================================================================== #
# ======================= BRAIN AREA ROLYPOLY *different 'SEM'* MOUSEBRAIN  RESULTS ================ #
# =========================================================================================================== #


### Load 'table.pvals.ALL_DATA.csv'
file.in <- "/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/export-combined.rp_hm3.v3.nboot100/inference_rp_hm3.body_BMI_Yengo2018.tss.10kb.hm3.none.protein_coding_only.nboot100/table.pvals.ALL_DATA.csv"
df.res_all <- read_csv(file.in)
df.res_all <- df.res_all %>% filter(grepl("mousebrain",dataset))
df.res_all <- df.res_all %>% mutate(dataset=stringr::str_replace(dataset, ".*\\.", "")) # shorten dataset names (rename)
# ALTERNATIVE SHORTEN? del: lapply(stringr::str_split(colnames(df.clean), "\\."), '[[', -1)
# ALTERNATIVE SHORTEN? del: stringr::str_match(colnames(df.clean), "\\.?(.*?)$") # match
df.res_all


# for (category.now in df %>% distinct(dataset) %>% pull()) {
enrichment_category <- function(category.now, df) {
  # print(category.now)
  # df.now <- df %>% mutate(category.tmp = if_else(category==category.now, category.now, "Other")) # SWITCH = region/category
  df.now <- df %>% mutate(category.tmp = if_else(Class==category.now, category.now, "Other")) # SWITCH = Class
  wtest <- wilcox.test(bt_value ~ category.tmp, data = df.now, alternative="greater")
  return(wtest$p.value)
}

list.res <- list()
# for (dataset.now in df.res_all %>% distinct(dataset) %>% pull()) { # all 'data sets'/SEM
for (dataset.now in c("enrichment_quantiles", "enrichment_log", "z_score_pos", "z_score_pos_log")) { 
  print(dataset.now)
  # df.dataset <- df.res_all %>% filter(dataset=="z_score_pos_log")
  # df.dataset <- df.res_all %>% filter(dataset=="specificity_quantiles")
  
  df.dataset <- df.res_all %>% filter(dataset==dataset.now) # filter(dataset==dataset) <-- DOES NOT WORK
  
  # df.wtest_res <- data.frame(enrichment=sapply(df.dataset %>% distinct(category) %>% pull(), enrichment_category, df=df.dataset)) # SWITCH = region/category
  df.wtest_res <- data.frame(enrichment=sapply(df.dataset %>% distinct(Class) %>% pull(), enrichment_category, df=df.dataset)) # SWITCH = Class
  df.wtest_res <- df.wtest_res %>% rownames_to_column(var="category")
  list.res[[dataset.now]] <- df.wtest_res
  
  # ggplot(df.wtest_res, aes(x=category, y=-log10(enrichment))) + geom_col() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title=dataset)
} 

df.list.res <- bind_rows(list.res, .id="dataset")

### Wrap data-sets
ggplot(df.list.res, aes(x=category, y=-log10(enrichment))) + geom_col() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~dataset)
### Select specific dataset
ggplot(df.list.res %>% filter(dataset=="z_score_pos_log"), aes(x=category, y=-log10(enrichment))) + geom_col() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





# area postrema of the medulla oblongata

### tmp 
# ggplot(df, aes(x=category.tmp, y=bt_value)) + geom_boxplot()


