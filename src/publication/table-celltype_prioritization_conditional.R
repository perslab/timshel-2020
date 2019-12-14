############### SYNOPSIS ###################
### AIM: Write conditional MB cell-type prioritization result tables

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
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))


# ======================================================================= #
# ================================= RUN ================================== #
# ======================================================================= #


### conditional data
# file.data <- here("results/prioritization_celltypes_conditional--mousebrain.BMI_UKBB_Loh2018.csv.gz")
file.data <- here("results/cellect_ldsc/conditional.csv")
df.cond <- read_csv(file.data)
df.cond <- format_cellect_ldsc_results(df.cond)

### baseline (prioritization)
file.data <- here("results/cellect_ldsc/prioritization.csv")
df.ldsc_cts <- read_csv(file.data)
df.ldsc_cts <- format_cellect_ldsc_results(df.ldsc_cts)
df.ldsc_cts <- df.ldsc_cts %>% filter(gwas == "BMI_UKBB_Loh2018")  # filter BMI results
df.ldsc_cts <- df.ldsc_cts %>% mutate(conditional_annotation="baseline")
df.ldsc_cts <- df.ldsc_cts %>% filter(specificity_id %in% unique(df.cond$specificity_id)) # keep only specific datasets


### combine + split
df <- bind_rows(df.cond, df.ldsc_cts)
list.df <- split(df, df$specificity_id)

### export
export_results <- function(df, name) {
  filter.annotations <- get_prioritized_annotations_bmi(dataset=name)
  df <- df %>% filter(conditional_annotation %in% c(filter.annotations, "baseline")) # keep only conditional_annotation for prioritized annotations
  # ^ this is necessary if file contains 'extra results'
  
  df.spread <- df %>% 
    select(annotation, conditional_annotation, p.value) %>%
    spread(key=conditional_annotation, value=p.value) %>%
    arrange(baseline)
  # print(df.spread, n=15)
  
  file.out <- here(sprintf("src/publication/tables/table-celltype_ldsc_conditional_results.%s.csv", name))
  message(sprintf("exporting file: %s", file.out))
  df.spread %>% write_csv(file.out)
}

### Run
imap(list.df, export_results)



