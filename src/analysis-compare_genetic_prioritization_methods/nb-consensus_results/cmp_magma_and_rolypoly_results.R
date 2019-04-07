
# ======================= Mousebrain: Correlation between MAGMA and RolyPoly 'mousebrain_skene.quantiles' ======================= #

library(tidyverse)

rm(list=ls()) # clear workspace

### MAGMA
file.in <- "/projects/timshel/sc-genetics/sc-genetics/src/nb-magma_celltyping/out.magma_celltyping.BMI_Yengo2018.results.mousebrain.csv"
df.magma <- read_csv(file.in)
df.magma <- df.magma %>% mutate(COVAR=stringr::str_match(df.magma$COVAR ,"_(.*)")[,2]) # rename Neurons_HBGLU3 to HBGLU3.

### tmp quick MAGMA plots
# ggplot(df.magma, aes(x=BETA, y=BETA_STD)) + geom_point() # --> MAGMA BETA and BETA_STD correlate nicely
# ggplot(df.magma, aes(x=-log10(P), y=BETA)) + geom_point() # --> MAGMA P-val and BETA correlate pretty well
# ggplot(df.magma, aes(x=SE)) + geom_density() # --> MAGMA have roughly a bell-shaped SE distribution with low SD.

### RolyPoly - 'table.pvals.ALL_DATA.csv'
file.in <- "/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/export-combined.rp_hm3.v3.nboot100/inference_rp_hm3.body_BMI_Yengo2018.tss.10kb.hm3.none.protein_coding_only.nboot100/table.pvals.ALL_DATA.csv"
df.rolypoly <- read_csv(file.in)
df.rolypoly <- df.rolypoly %>% filter(dataset=="mousebrain_skene.quantiles") # filter 

### Join data frames
df.combined <- left_join(df.rolypoly, df.magma, by=c("annotation"="COVAR"))
df.combined

### Correlation
with(df.combined, cor.test(x=-log10(bp_value), y=-log10(P))) # rho=0.56, p-value < 2.2e-16

### Plot
ggplot(df.combined, aes(x=-log10(bp_value), y=-log10(P))) + geom_point()

