############### SYNOPSIS ###################
### AIM: Make hypothalamus plot: Campbell+Mousebrain + BMI geneset cell-type enrichment

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

library(patchwork)

source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))

# ======================================================================= #
# ================== LOAD LDSC CTS RESULTS (multi GWAS) ================ #
# ======================================================================= #

### Read LDSC results
file.results <- here("results/cellect_ldsc/prioritization.csv")
df.ldsc_cts <- read_csv(file.results)

df.ldsc_cts %>% count(specificity_id)
df.ldsc_cts <- format_cellect_ldsc_results(df.ldsc_cts)
*INSERT CORRECT DATASET SUBSETTING HERE*

# df.ldsc_cts <- df.ldsc_cts %>% filter(specificity_id %in% c("campbell2017_lvl2", "moffitt2018", "romanov2017", "chen2017"))
  # df.ldsc_cts <- df.ldsc_cts %>% filter(specificity_id %in% c("mousebrain", "campbell2017_lvl2", "moffitt2018", "romanov2017", "chen2017"))
  # df.ldsc_cts <- df.ldsc_cts %>% filter(!specificity_id %in% c("modules.mousebrain_bmicelltypes")) # all single-cell datasets

# ======================================================================= #
# ============================ BMI point plot =============================== #
# ======================================================================= #

### Create 'plotting' data frame
df.ldsc_cts <- df.ldsc_cts %>% filter(gwas == "BMI_UKBB_Loh2018")  # filter BMI results


### Add fdr_significant flag (within GWAS)
df.ldsc_cts.tmp.summary <- df.ldsc_cts %>% group_by(specificity_id) %>% summarise(n_obs=n())
df.ldsc_cts <- left_join(df.ldsc_cts, df.ldsc_cts.tmp.summary, by="specificity_id")
df.ldsc_cts <- df.ldsc_cts %>% mutate(fdr_significant = if_else(p.value <= 0.05/n_obs, true=T, false=F))
df.ldsc_cts

df.ldsc_cts %>% filter(fdr_significant) %>% write_csv(here("src/publication/figs/tmp-priori.fdr_annotations.csv"))

df.ldsc_cts %>% filter(fdr_significant) %>% count(specificity_id)

### Add pvalue
df.plot <- df.ldsc_cts
df.plot <- df.plot %>% mutate(p.value.mlog10 = -log10(p.value))
df.plot <- df.plot %>% arrange(specificity_id) %>% mutate(annotation=factor(annotation, levels=unique(annotation)))
p <- ggplot(df.plot, aes(x=annotation, y=p.value.mlog10, color=specificity_id)) + 
  geom_point(aes(size=fdr_significant))
p
pp <- ggplotly(p)
pp
htmlwidgets::saveWidget(pp, here("src/publication/figs/tmp-prioritization.html"))
library(plotly)

