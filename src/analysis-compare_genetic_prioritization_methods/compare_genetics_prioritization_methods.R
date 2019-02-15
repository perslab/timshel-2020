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
library(gghighlight)


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
file.rolypoly <- "/nfsdata/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/export-combined.v3.nboot100/inference_rp.body_BMI_Locke2015.tss.10kb.pos_only.all_genes.nboot100/table.pvals.gtex.sub_tissue_gene_tpm.avg_expr.csv"
df.rolypoly <- read_csv(file.rolypoly)
df.rolypoly <- df.rolypoly %>% select(annotation, pval=bp_value)

# ======================================================================= #
# ================================ Combine data ================================ #
# ======================================================================= #

list.comb <- list("magma"=df.magma,
                  "ldsc"=df.ldsc,
                  "depict_zscore"=df.depict_zscore,
                  "depict"=df.depict)
df <- bind_rows(list.comb, .id="method")

# ======================================================================= #
# ================================ Export data ================================ #
# ======================================================================= #

df.export <- df %>% spread(key="method", value="pval")
df.export %>% write_csv("method_comparison.mousebrain.pvals.csv")
df.export.rank <- df %>% group_by(method) %>% mutate(pval=rank(pval)) %>% spread(key="method", value="pval")
df.export.rank %>% write_csv("method_comparison.mousebrain.rank.csv")


# ======================================================================= #
# ================================ PLOT CORRELATE / MATRIXPLOT ================================ #
# ======================================================================= #

### Matrix plot
df.matplot <- df %>% mutate(pval_mlog10 = -log10(pval)) 
df.matplot <- df.matplot %>% select(-pval) %>% spread(key="method", value="pval_mlog10")
ggpairs(df.matplot %>% select(-annotation))


# ======================================================================= #
# ================================ PLOT BARPLOT ================================ #
# ======================================================================= #

### Barplot - all cell-types
ggplot(df, aes(x=annotation, y=-log10(pval), fill=method)) + 
  geom_col(position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()

### Barplot - only FDR significant cell-types
filter.celltypes <- c("TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5","MEGLU10","TEGLU17","MEGLU11","TEGLU4","DEGLU4","TEINH12") # BMI_UKBB_Loh2018 FDR sign.
ggplot(df %>% filter(annotation %in%  filter.celltypes), aes(x=annotation, y=-log10(pval), fill=method)) + 
  geom_col(position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()



# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #

# p <- ggplot(df.ldsc_cts.spread, aes(x=continuous__sem_mean, y=raw__specificity_top10pct_binary)) + 
#   geom_point() +
#   gghighlight((continuous__sem_mean > 2), max_highlight = 4L, label_key = annotation) +
#   geom_abline()
# p



# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #




# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #




