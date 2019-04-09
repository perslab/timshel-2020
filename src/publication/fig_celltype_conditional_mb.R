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

source(here("src/lib/load_functions.R")) # load sc-genetics library

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
filter.annotations <- c("TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5","MEGLU10","TEGLU17","MEGLU11","TEGLU4","DEGLU4","TEINH12")
df.heatmap <- df %>% filter(annotation %in% filter.annotations)

# x = condition
# y = annotation
### TODO
# TODO: color: same as main figure
# TODO: mark asterisk for FDR cutoff from main fig.
# TOOD: order conditions by their CLUSTERING (move "baseline" first)
# TODO: add NA values for "self comparisons"
# NOT: DO NOT ADD FACET SPACER. DO IT IN ILLUSTRATOR.

p <- ggplot(df.heatmap, aes(x=condition, y=annotation, fill=-log10(p.value))) +
  geom_tile() + 
  theme(axis.text.x = element_text(angle=45, hjust=1))
p

# ======================================================================= #
# ========= PLOT: cell priori [FACET WRAP MULTI-CONDITIONAL] ============== #
# ======================================================================= #


df.plot <- df.ldsc_cts
fdr_threshold <- 0.05/nrow(df.plot %>% filter(condition==(df.plot %>% pull(condition))[1])) # *OBS*: TMP DIRTY CODE FOR CONDITIONAL
p <- ggplot(df.plot, aes(x=annotation, y=-log10(p.value))) +
  geom_point(aes(color=color_by_variable, size=estimate)) + 
  ggrepel::geom_text_repel(data=df.plot %>% filter(fdr_significant), aes(x=annotation, y=-log10(p.value), label=annotation, color=color_by_variable), hjust = 0, nudge_x = 1.5) + # OBS geom_text_repel
  geom_hline(yintercept=-log10(fdr_threshold), color="red") + 
  facet_wrap(~condition, ncol = 2) + 
  labs(x="Cell Type", y="-log10(P-value)") + # expression(-log[10]("P-value"))
  # ggtitle(sprintf("%s - %s", "BMI_UKBB_Loh2018", name_elem)) +
  # ggtitle(sprintf("%s - %s", dataset_prefix, name_elem)) +
  theme_classic() + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p
# file.out <- sprintf("out.tmp.190212.plot.cell_prioritization.%s.%s.sem_mean.color_by_class.pdf", dataset_prefix, "CONDITIONAL_BMI_UKBB_Loh2018")
# ggsave(p, filename=file.out, width=15, height=20)
