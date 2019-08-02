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
# ================================ PLOT SCATTERPLOT ================================ #
# ======================================================================= #

df.plot <- df %>% 
  mutate(pval = -log10(pval)) %>%
  spread(key="method", value="pval")
df.plot

### Cell-types to highlight
filter.annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")
colormap.annotations <- get_color_mapping.prioritized_annotations_bmi(dataset="mousebrain")

### PLOT
fdr_threshold.mlog10 <- -log10(0.05/nrow(df.plot))
p <- ggplot(df.plot, aes(x=ldsc, y=magma)) + 
  geom_point(color="gray") + 
  ### Lines
  geom_abline() + 
  geom_hline(yintercept=fdr_threshold.mlog10, linetype="dashed", color="gray") +
  geom_vline(xintercept=fdr_threshold.mlog10, linetype="dashed", color="gray") +
  ### MAGMA FDR
  geom_point(data=df.plot %>% filter((magma > fdr_threshold.mlog10) & (ldsc < fdr_threshold.mlog10)), color="black", size=3) +
  geom_text_repel(data=df.plot %>% filter((magma > fdr_threshold.mlog10) & (ldsc < fdr_threshold.mlog10)), aes(label=annotation), color="black") + 
  # LDSC FDR
  geom_point(data=df.plot %>% filter(annotation %in%  filter.annotations), aes(color=annotation), size=3) +
  geom_text_repel(data=df.plot %>% filter(annotation %in%  filter.annotations), aes(label=annotation, color=annotation)) + 
  labs(x=expression(-log[10](P[S-LDSC])), y=expression(-log[10](P[MAGMA]))) + 
  scale_color_manual(values=colormap.annotations) + 
  guides(color=F) + 
  ggpubr::stat_cor(method="pearson")
p

file.out <- "figs/fig_celltypepriori_method_comparison.mb_fdr_scatter.pdf"
ggsave(plot=p, filename=file.out, width=8, height=5)


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




