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
filter.annotations <- c("TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5","MEGLU10","TEGLU17","MEGLU11","TEGLU4","DEGLU4","TEINH12")

### SELECTED modules
filter.modules <- "lavenderblush" # "lightpink3"

### Select module + cell-types
df.gather <- df %>% 
  select(module_id, filter.annotations) %>%
  filter(module_id == filter.modules) %>%
  gather(key="annotation", value="p.value", -module_id)

# ======================================================================= #
# =============================== MAKE PLOT ============================== #
# ======================================================================= #

p <- ggplot(df.gather, aes(x=annotation, y=-log10(p.value))) + 
  geom_col() + 
  geom_hline(yintercept=-log10(0.05/n_distinct(df.gather$annotation))) + 
  theme(axis.text.x=element_text(angle=45, hjust=1))
p



