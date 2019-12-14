############### SYNOPSIS ###################
# Plot conditional analysis for MB
# 


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

library(colorspace)

source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))

# ======================================================================= #
# ============================ LOAD DATA =============================== #
# ======================================================================= #

# file.data <- here("results/prioritization_celltypes_conditional--mousebrain.BMI_UKBB_Loh2018.csv.gz")
file.data <- here("results/cellect_ldsc/conditional.csv")
df.cond <- read_csv(file.data)
df.cond <- df.cond %>% rename(p.value=pvalue)

# filter.exclude <- "wgcna.mousebrain-190111-dodgerblue"
# df <- df %>% filter(!annotation %in% filter.exclude,
#                     !condition %in% filter.exclude)

### baseline (prioritization)
file.data <- here("results/cellect_ldsc/prioritization.csv")
df.ldsc_cts <- read_csv(file.data)
df.ldsc_cts <- format_cellect_ldsc_results(df.ldsc_cts)
df.ldsc_cts <- df.ldsc_cts %>% filter(gwas == "BMI_UKBB_Loh2018")  # filter BMI results
df.ldsc_cts <- df.ldsc_cts %>% mutate(conditional_annotation="baseline")
df.ldsc_cts <- df.ldsc_cts %>% filter(specificity_id %in% unique(df.cond$specificity_id)) # keep only specific datasets


### combine
df <- bind_rows(df.cond, df.ldsc_cts)

# ======================================================================= #
# ============================= PLOT - HEATMAP ========================== #
# ======================================================================= #

### SELECTED ANNOTATIONS
filter.annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")
df.heatmap <- df %>% filter(annotation %in% filter.annotations)

### ADD "block"
df.heatmap <- df.heatmap %>% mutate(block = if_else(conditional_annotation=="baseline", "A", "B"))

### Set 'diagonal' to NA ('conditional on self')
df.heatmap <- df.heatmap %>% mutate(p.value = if_else(as.character(conditional_annotation) == as.character(annotation), NA_real_, p.value))

### ANNOTATIONS ORDERED BY HIERARCHICAL CLUSTERING OF ESmu
# annotations.ordered <- c(
#                  "TEGLU23", # HC
#                  "TEGLU4", # CORTEX
#                  "TEGLU17", # CORTEX
#                  "TEINH12", # HC/CORTEX
#                  "DEGLU5", # MIDBRAIN
#                  "DEINH3", # HYPOTHALAMUS
#                  "MEGLU1", # MIDBRAIN
#                  "MEINH2", # MIDBRAIN
#                  "MEGLU11", # MIDBRAIN
#                  "MEGLU10", # MIDBRAIN
#                  "DEGLU4" # THALAMUS
#                  )

annotations.ordered <- c(
  "TEGLU21",
  "TEGLU23",
  "TEGLU4",
  "TEGLU17",
  "TEGLU19",
  "TEINH12",
  "DEGLU4",
  "DEINH3",
  "DEGLU5",
  "MEGLU1",
  "MEINH2",
  "MEGLU11",
  "MEGLU10",
  "MEGLU2",
  "MEGLU3",
  "MEGLU7",
  "MEINH3",
  "TEINH1",
  "HBINH8",
  "HBGLU2",
  "HBGLU5",
  "HBSER5")

### Order "annotation" by their hierarchical clustering 
df.heatmap <- df.heatmap %>% mutate(annotation = factor(annotation, levels=annotations.ordered))
### Order "conditional_annotation" as above +move "baseline" first
df.heatmap <- df.heatmap %>% mutate(conditional_annotation = factor(conditional_annotation, levels=c("baseline", annotations.ordered)))

### Get annotation colors
colormap.annotations <- get_color_mapping.prioritized_annotations_bmi(dataset="mousebrain")

PVAL_THRESHOLD <- 0.05/265 # pvalue cut-off in main MB fig
# PVAL_THRESHOLD <- 0.05/11 
p <- ggplot(df.heatmap, aes(x=conditional_annotation, y=annotation, fill=-log10(p.value))) +
  geom_tile() + 
  geom_text(data=df.heatmap %>% filter(p.value <= PVAL_THRESHOLD), label="*", color="black", hjust=0.5, vjust=0.75) + # add asterisk if fdr significant
  colorspace::scale_fill_continuous_sequential(palette="Blues 2", rev=TRUE, na.value = "white") + # "Blues 2","Blues 3","Purples 3"
  labs(x="Conditional cell-type", y="", fill=expression(-log[10](P[S-LDSC]))) +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) + # remove grid
  theme(axis.text.x=element_text(angle=45, hjust=1))
### Add text color to annotations
colxaxis <- colormap.annotations[match(levels(df.heatmap$annotation), names(colormap.annotations))]
p <- p + theme(axis.text.y = element_text(color=colxaxis), 
               axis.text.x = element_text(color=colxaxis))
### Facet spacing of baseline
p <- p + facet_grid(cols=vars(block), space="free_x", scales="free_x", drop=T) + 
  theme(strip.background=element_blank(), # remove facet strip bg REF: https://stackoverflow.com/questions/14185754/remove-strip-background-keep-panel-border
        strip.text.x=element_blank())
p

file.out <- sprintf("figs/fig_celltypepriori_mb_conditional-4x6.pdf")
ggsave(p, filename=file.out, height=4, width=6)
file.out <- sprintf("figs/fig_celltypepriori_mb_conditional.pdf")
ggsave(p + theme(legend.position="top"), filename=file.out, height=6.5, width=5)

# ======================================================================= #
# ========= PLOT: cell priori [FACET WRAP MULTI-CONDITIONAL] ============== #
# ======================================================================= #


# df.plot <- df.ldsc_cts
# fdr_threshold <- 0.05/nrow(df.plot %>% filter(condition==(df.plot %>% pull(condition))[1])) # *OBS*: TMP DIRTY CODE FOR CONDITIONAL
# p <- ggplot(df.plot, aes(x=annotation, y=-log10(p.value))) +
#   geom_point(aes(color=color_by_variable, size=estimate)) + 
#   ggrepel::geom_text_repel(data=df.plot %>% filter(fdr_significant), aes(x=annotation, y=-log10(p.value), label=annotation, color=color_by_variable), hjust = 0, nudge_x = 1.5) + # OBS geom_text_repel
#   geom_hline(yintercept=-log10(fdr_threshold), color="red") + 
#   facet_wrap(~condition, ncol = 2) + 
#   labs(x="Cell Type", y=expression(-log[10](P))) +
#   theme_classic() + 
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())
# p
# file.out <- sprintf("out.tmp.190212.plot.cell_prioritization.%s.%s.es_mean.color_by_class.pdf", dataset_prefix, "CONDITIONAL_BMI_UKBB_Loh2018")
# ggsave(p, filename=file.out, width=15, height=20)
