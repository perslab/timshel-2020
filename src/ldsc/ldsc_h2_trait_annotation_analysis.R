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


source(here("src/lib/load_functions.R")) # load sc-genetics library


setwd(here("src/ldsc"))

# ======================================================================= #
# ================================ LOAD DATA ================================ #
# ======================================================================= #

### Load - MULTI GWAS AND ANNOTATIONS
dir.data <- "/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc_h2/"
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
# =========================== EXPORT to results ========================= #
# ======================================================================= #

### Export (selected columns)
file.out <- here("results/h2_annotations.multi_gwas.csv.gz")
file.out
# df.ldsc %>% select(-Category) %>% arrange(gwas) %>% write_csv(file.out)


# ======================================================================= #
#============================= COMBINE TABLE ========================== #
# ======================================================================= #

df.export <- df.ldsc %>% group_by(gwas) %>% 
  mutate(rank = rank(-`Coefficient_z-score`)) %>% 
  arrange(gwas, rank) %>%
  select(-Category)
df.export
df.export %>% write_csv("out.h2_annotations.all.csv")

df.export <- df.ldsc %>% group_by(gwas) %>% 
  mutate(rank = rank(-`Coefficient_z-score`)) %>% 
  top_n(n=5, wt=desc(rank)) %>%
  arrange(gwas, rank) %>%
  select(-Category)
df.export
df.export %>% write_csv("out.h2_annotations.top5.csv")

### export WGCNA
df.export <- df.ldsc %>% 
  filter(run_name=="wgcna.mousebrain-190111") %>%
  select(-Category)
df.export %>% write_csv("out.h2_annotations.wgcna.csv")


# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #



# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #
