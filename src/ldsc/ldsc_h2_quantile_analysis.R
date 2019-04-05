############### SYNOPSIS ###################
### AIM: Perform h2 quantile analysis (h2 enrichment for each quantile of a continuous annotation)

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


setwd(here("src/ldsc"))

# ======================================================================= #
# ================================ LOAD DATA ================================ #
# ======================================================================= #

# This output file has one row for each quantile (starting with lowest values) 
# column summarizing the heritability explained by each quantile, its enrichment and corresponding standard error and P value.

### Load - MULTI GWAS AND ANNOTATIONS
dir.data <- "/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc_h2/"
filenames <- list.files(path=dir.data,  pattern="(.*).results_quantile_qfixed_h2")
# filenames <- filenames[!grepl(pattern="__CONDITIONAL__", filenames, perl=T)] # exclude any conditional results
list.dfs <- lapply(file.path(dir.data, filenames), read_tsv)
list.dfs <- lapply(list.dfs, function(df) {df %>% mutate(q=paste0("q", seq(0,n()-1)))}) # add quantile
names(list.dfs) <- stringr::str_match(filenames, pattern="(.*).results_quantile_qfixed_h2")[,2] 
names(list.dfs)
df.ldsc <- list.dfs %>% bind_rows(.id="run_str_full")

### Split sting
df.ldsc <- df.ldsc %>% separate(col=run_str_full, into=c("run_name", "annotation", "gwas"), sep="__") # e.g. run_str = "wgcna.mousebrain-190111__wheat3__T2D_UKBB_DIAMANTE_Mahajan2018"
### Filter
# df.ldsc <- df.ldsc %>% filter(gwas %in% c("BMI_UKBB_Loh2018", "T2D_DIAMANTE_Mahajan2018"))


# ======================================================================= #
# =========================== EXPORT to results ========================= #
# ======================================================================= #

### Export (selected columns)
file.out <- here("results/h2_annotation_intervals-multi_gwas.csv.gz")
# file.out
# df.ldsc %>% arrange(gwas) %>% write_csv(file.out)


# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #



# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #
