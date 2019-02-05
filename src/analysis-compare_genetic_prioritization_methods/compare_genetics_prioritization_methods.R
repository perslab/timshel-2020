############### SYNOPSIS ###################

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

dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib/"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library


setwd(here("src/analysis-compare_genetic_prioritization_methods"))

# ======================================================================= #
# ================================ CONSTANTS ================================ #
# ======================================================================= #

# ONLY MAKE FIGURE FOR BMI_UKBB_Loh2018_no_mhc

gwas_name <- "BMI_UKBB_Loh2018_no_mhc"
# gwas_name <- "BMI_UPDATE_Yengo2018_no_mhc"
dataset_prefix <- "mousebrain"
# dataset_prefix <- "tabula_muris"

# ======================================================================= #
# ================================ LDSC ================================ #
# ======================================================================= #
file.ldsc <- here("out/out.ldsc/celltypes.mousebrain.all__BMI_UKBB_Loh2018.cell_type_results.txt")
df.ldsc <- load_ldsc_cts_results(file.ldsc, dataset_prefix="mousebrain.all")
df.ldsc <- df.ldsc %>% filter(sem=="sem_mean") # just make sure we only have sem_mean results
df.ldsc <- df.ldsc %>% select(annotation, pval=p.value)
df.ldsc
# ======================================================================= #
# ================================ MAGMA ================================ #
# ======================================================================= #
file.magma <- here("src/magma/out.cell_prioritization.mousebrain.BMI_UKBB_Loh2018_no_mhc.sem_meta_mean.csv")
df.magma <- read_csv(file.magma)
df.magma <- df.magma %>% select(annotation, pval=p.value)
df.magma

# ======================================================================= #
# ================================ DEPICT ================================ #
# ======================================================================= #
file.depict <- here("src/depict/results/BMI_UKBB_Loh2018_no_mhc.5e-8.mousebrain_sem_mean_tissueenrichment.txt")
df.depict <- read_tsv(file.depict)
# problems(df.depict) # ---> we accept the parsing errors and correct for them below
df.depict <- df.depict %>% select(annotation=`MeSH term`,
                     pval=Name) # fdr=`MeSH first level term`

file.depict_zscore <- here("src/depict/results/BMI_UKBB_Loh2018_no_mhc.5e-8.mousebrain_z_score_two_step_tissueenrichment.txt")
df.depict_zscore <- read_tsv(file.depict_zscore)
df.depict_zscore <- df.depict_zscore %>% select(annotation=`MeSH term`,
                                  pval=Name)



# ======================================================================= #
# ================================ ROLYPPOLY ================================ #
# ======================================================================= #
# read CSV in this folder


# ======================================================================= #
# ================================ Combine data ================================ #
# ======================================================================= #

list.comb <- list("magma"=df.magma,
                  "ldsc"=df.ldsc,
                  "depict_zscore"=df.depict_zscore,
                  "depict"=df.depict)
df <- bind_rows(list.comb, .id="method")


# ======================================================================= #
# ================================ PLOT CORRELATE / MATRIXPLOT ================================ #
# ======================================================================= #

### Matrix plot
df.matplot <- df %>% mutate(pval_mlog10 = -log10(pval)) 
df.matplot <- df.matplot %>% select(-pval) %>% spread(key="method", value="pval_mlog10")
ggpairs(df.matplot %>% select(-annotation))

# ggcorr()

ggplot(df, aes(x=annotation, y=-log10(pval), fill=method)) + geom_col(position=position_dodge())



# ======================================================================= #
# ================================ PLOT BARPLOT ================================ #
# ======================================================================= #

# TODO: filter only FDR significant cell-types

ggplot(df, aes(x=annotation, y=-log10(pval), fill=method)) + geom_col(position=position_dodge())


# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #




# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #




# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #




