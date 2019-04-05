############### SYNOPSIS ###################
# Analyze and export LDSC results for multiple datasets (tabula muris, mousebrain, ...)
# JOINT!!!


### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# =============================== SETUP ================================= #
# ======================================================================= #

library(tidyverse)
library(here)

source(here("src/lib/load_functions.R")) # load sc-genetics library

setwd(here("src/ldsc"))

# ======================================================================= #
# ============================ PARAMETERS =============================== #
# ======================================================================= #


### PARAMS
dataset_prefix <- "mousebrain_all"
# dataset_prefix <- "tabula_muris"

genomic_annotation_prefix <- get_genomic_annotation_prefix(dataset_prefix)

# ======================================================================= #
# =============================== LOAD LDSC CTS RESULTS ================================= #
# ======================================================================= #

file.ldsc_cts <- sprintf("/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/%s__%s.cell_type_results.txt", genomic_annotation_prefix, gwas)
df.ldsc_cts <- load_ldsc_cts_results(file.ldsc_cts, dataset_prefix)
df.ldsc_cts <- df.ldsc_cts %>% filter(sem=="sem_mean")
df.ldsc_cts

# ======================================================================= #
# =============================== LOAD LDSC H2 JOINT ================================= #
# ======================================================================= #

file.h2_joint <- "/nfsdata/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc_h2_joint/celltypes.mousebrain.all__BMI_UKBB_Loh2018__JOINT__TEGLU23--DEINH3--MEGLU1--MEINH2--DEGLU5--MEGLU10--TEGLU17--MEGLU11--TEGLU4--DEGLU4--TEINH12.results"
df.h2_joint <- read_tsv(file.h2_joint) %>%
  rename(tau_zscore=`Coefficient_z-score`) %>%
  mutate(tau_pval=pnorm(abs(tau_zscore),lower.tail=F))

df.h2_joint %>% arrange(desc(tau_zscore))

### PLOT

ggplot(df.h2_joint, aes(x=Category, y=tau_zscore)) + 
  geom_col(width=.5, fill="tomato3") + 
  geom_hline(yintercept = 1.96, color="black") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6))


ggplot(df.h2_joint, aes(x=Category, y=tau_zscore, label=round(tau_zscore,digits=2))) + 
  geom_point(stat='identity', fill="black", size=6)  +
  geom_segment(aes(y = 0,
                   yend = tau_zscore, 
                   x = Category, 
                   xend = Category), 
               color = "black") +
  geom_hline(yintercept = 1.96, color="gray", linetype="dashed") +
  geom_text(color="white", size=2) +
  labs(title="LDSC joint model", subtitle="with N=11 FDR cell-types - BMI_UKBB_Loh2018") + 
  coord_flip()
pnorm(1.98, lower.tail = F)


# ======================================================================= #
# =============================== COMPARE JOINT TO 'NON-JOINT' (CTS) ================================= #
# ======================================================================= #

### df.h2_joint: extract Mousebrain cell-types and clean names 
df.h2_joint.mb <- df.h2_joint %>% 
  filter(grepl("mousebrain_all",Category)) %>%
  mutate(Category=str_split_fixed(Category, "\\.", n=Inf)[,2])

### join dfs
df.cmp <- df.h2_joint.mb %>% left_join(df.ldsc_cts, by=c("Category"="annotation"))
df.cmp

### plot
ggplot(df.cmp, aes(x=-log10(tau_pval), y=-log10(p.value), label=Category)) + 
  geom_point() + 
  geom_text() + 
  geom_abline() + 
  labs(x="-log10(pval joint)", y="-log10(pval non-joint (cts))", title="'NON-JOINT' (CTS) vs JOINT")

# ======================================================================= #
# =============================== XXXXX ================================= #
# ======================================================================= #

