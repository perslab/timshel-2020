############### SYNOPSIS ###################
### AIM: Compare h2 (enrichment) for top cell-types/annotation for different traits.

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

library(GGally)
library(corrr) # devtools::install_github("drsimonj/corrr")
library(gghighlight)


dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib/"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library


setwd(here("src/ldsc"))

# ======================================================================= #
# ================================ LOAD DATA ================================ #
# ======================================================================= #

### Load - MULTI GWAS AND ANNOTATIONS
dir.data <- "/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc_h2/"
filenames <- list.files(path=dir.data,  pattern="(.*).results$")
list.dfs <- lapply(file.path(dir.data, filenames), read_tsv)
list.dfs <- lapply(list.dfs, function(df) {df %>% filter(row_number() == n())}) # extract last row
  # df %>% last() # --> ok, but returns vector
  # df %>% filter(row_number() == n())
names(list.dfs) <- stringr::str_match(filenames, pattern="(.*).results")[,2] 
names(list.dfs)
df.ldsc <- list.dfs %>% bind_rows(.id="run_str_full")

### Split sting
df.ldsc <- df.ldsc %>% separate(col=run_str_full, into=c("run_name", "annotation", "gwas"), sep="__") # e.g. run_str = "wgcna.mousebrain-190111__wheat3__T2D_UKBB_DIAMANTE_Mahajan2018"
### Filter
# df.ldsc <- df.ldsc %>% filter(gwas %in% c("BMI_UKBB_Loh2018", "T2D_DIAMANTE_Mahajan2018"))
df.ldsc
# ======================================================================= #
#============================= COMBINE TABLE ========================== #
# ======================================================================= #

df.export <- df.ldsc %>% group_by(gwas) %>% 
  mutate(rank = rank(-`Coefficient_z-score`)) %>% 
  top_n(n=5, wt=desc(rank)) %>%
  arrange(gwas, rank) %>%
  select(-Category)
df.export
df.export %>% write_csv("out.h2_trait_annotation_table.csv")

### export WGCNA
df.export <- df.ldsc %>% 
  filter(run_name=="wgcna.mousebrain-190111") %>%
  select(-Category)
df.export %>% write_csv("out.h2_trait_annotation_table.wgcna.csv")


# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #



# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #
