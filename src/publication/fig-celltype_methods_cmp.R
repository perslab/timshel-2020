############### SYNOPSIS ###################
### AIM: compare LDSC and MAGMA results

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


source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))

# ======================================================================= #
# ================================ CONSTANTS ================================ #
# ======================================================================= #

# ONLY MAKE FIGURE FOR BMI_UKBB_Loh2018_no_mhc

gwas_name <- "BMI_UKBB_Loh2018_no_mhc"
dataset_prefix <- "mousebrain"


# ======================================================================= #
# ================================ LDSC ================================ #
# ======================================================================= #
file.ldsc <- here("out/ldsc/cts/celltypes.mousebrain.all__BMI_UKBB_Loh2018.cell_type_results.txt")
df.ldsc <- load_ldsc_cts_results(file.ldsc, dataset_prefix="mousebrain.all")
df.ldsc <- df.ldsc %>% filter(sem=="sem_mean") # just make sure we only have sem_mean results
df.ldsc <- df.ldsc %>% select(annotation, pval=p.value)
df.ldsc
# ======================================================================= #
# ================================ MAGMA ================================ #
# ======================================================================= #
# file.magma <- here("src/model-fit_gene_based_scores/out.cell_prioritization.mousebrain.BMI_UKBB_Loh2018_no_mhc.sem_meta_mean.csv")
file.magma <- here("out/magma/celltype_models/out.cell_prioritization.mousebrain.BMI_UKBB_Loh2018_no_mhc.sem_meta_mean.csv")
df.magma <- read_csv(file.magma)
df.magma <- df.magma %>% select(annotation, pval=p.value)
df.magma


# ======================================================================= #
# ============================= COMBINE DATA ============================ #
# ======================================================================= #
list.comb <- list("magma"=df.magma,
                  "ldsc"=df.ldsc)
df <- bind_rows(list.comb, .id="method")

# ======================================================================= #
# ================================ Export data ================================ #
# ======================================================================= #

df.export <- df %>% spread(key="method", value="pval")
file.out <- sprintf("tables/table_celltype_priori_method_comparison.mousebrain.pvals.csv")
# df.export %>% write_csv(file.out)

# ======================================================================= #
# ================================ PLOT BARPLOT ================================ #
# ======================================================================= #

### Order annotation by LDSC significance
order.annotation <- df %>% filter(method=="ldsc") %>% arrange(pval) %>% pull(annotation) 
df <- df %>% mutate(annotation=fct_relevel(annotation, order.annotation))
df$annotation %>% levels()

### Barplot - all cell-types
ggplot(df, aes(x=annotation, y=-log10(pval), fill=method)) + 
  geom_col(position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()

### Barplot - only FDR significant cell-types
filter.annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")
p <- ggplot(df %>% filter(annotation %in%  filter.annotations), aes(x=annotation, y=-log10(pval), fill=method)) + 
  geom_col(width=0.4, position=position_dodge(width=0.5)) +
  geom_hline(yintercept=-log10(0.05/265), linetype="dashed", color="gray") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip() +
  labs(y=expression(-log[10](P)))
p
file.out <- "figs/fig_celltypepriori_method_comparison.mb_fdr.pdf"
ggsave(plot=p, filename=file.out, width=10, height=6)



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




