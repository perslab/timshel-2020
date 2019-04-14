############### SYNOPSIS ###################
### Module cell-type enrichment plot

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

library(RColorBrewer)

source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))


# ======================================================================= #
# =============================== LOAD DATA ============================== #
# ======================================================================= #

### JT April 8th 2019
# Module-cell-type enrichment - wilcox analytical p-values:
# matrix of p-values (celltypes vs all WGCNA modules) 
file.data <- "/projects/jonatan/applied/18-mousebrain_7/tables/mb_GSA_SEMall_WGCNAtop_2_analyt_mat_wilcoxonPval_genesetTests.txt"
df <- read.table(file.data, sep="\t", header=T) %>% rownames_to_column("module_id") %>% as_tibble
df


# ======================================================================= #
# ============================= SELECT DATA ============================= #
# ======================================================================= #

### SELECTED ANNOTATIONS
filter.annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")

### SELECTED modules
filter.modules <- "lavenderblush" # "lightpink3"

### Select module + cell-types
df.gather <- df %>% 
  select(module_id, filter.annotations) %>%
  filter(module_id == filter.modules) %>%
  gather(key="annotation", value="p.value", -module_id)

### Add columns
df.gather <- df.gather %>% mutate(p.value.mlog10=-log10(p.value))

# ======================================================================= #
# =============================== MAKE PLOT ============================== #
# ======================================================================= #

p <- ggplot(df.gather, aes(x=annotation, y=p.value.mlog10)) + 
  # geom_col() + 
  geom_hline(yintercept=-log10(0.05/n_distinct(df.gather$annotation)), linetype="dashed", color="gray") + 
  geom_segment(aes(x=annotation, xend=annotation, y=0, yend=p.value.mlog10), color="grey") +
  geom_point(aes(color=p.value.mlog10), size=5) + 
  scale_color_viridis_c(direction=-1) + 
  labs(x="", y=expression(-log[10](P[enrichment]))) + # title="M1 genes enrichment in cell-type ES genes"
  guides(color=FALSE) +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  theme(plot.title = element_text(size = rel(0.6))) +
  coord_flip()
p


file.out <- "figs/fig_module_celltype_enrichment.pdf"
ggsave(file.out, plot=p, width=4, height=5)


