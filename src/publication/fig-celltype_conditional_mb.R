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

file.data <- here("results/prioritization_celltypes_conditional--mousebrain_all.BMI_UKBB_Loh2018.csv.gz")
df <- read_csv(file.data)

filter.exclude <- "wgcna.mousebrain-190111-dodgerblue"
df <- df %>% filter(!annotation %in% filter.exclude,
                    !condition %in% filter.exclude)

# ======================================================================= #
# ============================= PLOT - HEATMAP ========================== #
# ======================================================================= #

### SELECTED ANNOTATIONS
filter.annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")
df.heatmap <- df %>% filter(annotation %in% filter.annotations)

### ADD "block"
df.heatmap <- df.heatmap %>% mutate(block = if_else(condition=="baseline", "A", "B"))

### Set 'diagonal' to NA ('conditional on self')
df.heatmap <- df.heatmap %>% mutate(p.value = if_else(as.character(condition) == as.character(annotation), NA_real_, p.value))

### ANNOTATIONS ORDERED BY HIERARCHICAL CLUSTERING OF ESmu
annotations.ordered <- c(
                 "TEGLU23", # HC
                 "TEGLU4", # CORTEX
                 "TEGLU17", # CORTEX
                 "TEINH12", # HC/CORTEX
                 "DEGLU5", # MIDBRAIN
                 "DEINH3", # HYPOTHALAMUS
                 "MEGLU1", # MIDBRAIN
                 "MEINH2", # MIDBRAIN
                 "MEGLU11", # MIDBRAIN
                 "MEGLU10", # MIDBRAIN
                 "DEGLU4" # THALAMUS
                 )

### Order "annotation" by their hierarchical clustering 
df.heatmap <- df.heatmap %>% mutate(annotation = factor(annotation, levels=annotations.ordered))
### Order "condition" as above +move "baseline" first
df.heatmap <- df.heatmap %>% mutate(condition = factor(condition, levels=c("baseline", annotations.ordered)))

### Get annotation colors
colormap.annotations <- get_color_mapping.prioritized_annotations_bmi(dataset="mousebrain")

# PVAL_THRESHOLD <- 0.05/265 # pvalue cut-off in main MB fig
# PVAL_THRESHOLD <- 0.05/11 
p <- ggplot(df.heatmap, aes(x=condition, y=annotation, fill=-log10(p.value))) +
  geom_tile() + 
  # geom_text(data=df.heatmap %>% filter(p.value >= PVAL_THRESHOLD), label="*", color="black", hjust=0.5, vjust=0.75) + # add asterisk if fdr significant
  colorspace::scale_fill_continuous_sequential(palette="Blues 2", rev=TRUE, na.value = "white") + # "Blues 2","Blues 3","Purples 3"
  labs(x="Conditional cell-type", y="", fill=expression(-log[10](P))) +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) + # remove grid
  theme(axis.text.x=element_text(angle=45, hjust=1))
### Add text color to annotations
colxaxis <- colormap.annotations[match(levels(df.heatmap$annotation), names(colormap.annotations))]
p <- p + theme(axis.text.y = element_text(color=colxaxis))
### Facet spacing of baseline
p <- p + facet_grid(cols=vars(block), space="free_x", scales="free_x", drop=T) + 
  theme(strip.background=element_blank(), # remove facet strip bg REF: https://stackoverflow.com/questions/14185754/remove-strip-background-keep-panel-border
        strip.text.x=element_blank())
p

file.out <- sprintf("figs/fig_celltypepriori_mb_conditional.pdf")
ggsave(p, filename=file.out, width=6, height=4)

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
# file.out <- sprintf("out.tmp.190212.plot.cell_prioritization.%s.%s.sem_mean.color_by_class.pdf", dataset_prefix, "CONDITIONAL_BMI_UKBB_Loh2018")
# ggsave(p, filename=file.out, width=15, height=20)
