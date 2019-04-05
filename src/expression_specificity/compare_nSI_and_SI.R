############### SYNOPSIS ###################
### AIM:
# Compare nSI and SI values

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

source(here("src/lib/load_functions.R")) # load sc-genetics library

setwd(here("src/expression_specificity"))


# ======================================================================= #
# ================================ LOAD ES ================================ #
# ======================================================================= #

file.rdata.es <- here("src/datasets-expression/tabula_muris/tabula_muris.sem_obj-190306.RData") # sem_obj
load(file.rdata.es)

df.nsi <- sem_obj$sem$si %>% mutate(gene=sem_obj$genes)

# ======================================================================= #
# ================================ LOAD pSI ================================ #
# ======================================================================= #

file.rdata.psi <- here("/src/expression_specificity/pSI/psi.tabula_muris_tissuecell_type.pmax_1.e_min_1e-5.RData") # list.si + df.avg_expr
load(file.rdata.psi)

df.si <- list.si$SI %>% rownames_to_column(var="gene")


# ======================= COPY OF CALL TO specificity.index.timshel ======================= #
### *OBS*:  we here use 'tabula_muris.tissue_celltype.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz' which is preproccesed in a slightly different why than our 'new' python preprocessing

### DOCS: specificity.index.timshel() is JUST a COPY of CRAN specificity.index() that return both pSI and SI values as a list [no SI = FALSE argument]
# source("/raid5/projects/timshel/sc-genetics/sc-genetics/src/nb-sem/specificity.index.timshel.R") # loads 'specificity.index.timshel' function
### CALCULATION
# file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/tabula_muris/tabula_muris.tissue_celltype.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz" # tabula_muris
# df.avg_expr <- read_csv(file.in) %>% column_to_rownames("gene") %>% as.data.frame()
# system.time(list.si <- specificity.index.timshel(pSI.in=df.avg_expr, p_max = 1, e_min=1e-5))
# list(SI=df.SI, pSI=df.pSI) # named list
# save.image("psi.tabula_muris_tissuecell_type.pmax_1.e_min_1e-5.RData")
### RUNTIME
# MACA  p_max = 1, e_min=1e-5: 13.7h


# ======================================================================= #
# ================================ JOIN ================================ #
# ======================================================================= #

names(df.si) # column names have been 'make.names': e.g. "Thymus.DN1.thymic.pro.T.cell"
names(df.nsi) <- make.names(names(df.nsi)) # do the same for df.nsi

# ======================================================================= #
# ================================ COMPARE ================================ #
# ======================================================================= #

### Get random annotation
annotation.selected <- sample(names(df.nsi), size=1)
annotation.selected

### Join
df.selected <- inner_join(df.nsi %>% select(gene, !!sym(annotation.selected)),
           df.si %>% select(gene, !!sym(annotation.selected)),
           by="gene", suffix=c(".nsi", ".si")) %>% 
  rename(nsi=paste0(annotation.selected, ".nsi"),
         si=paste0(annotation.selected, ".si"))
df.selected

### Plot
p <- ggplot(df.selected, aes(x=nsi, y=si)) + 
  geom_point() + 
  geom_density_2d() +
  labs(title=annotation.selected) + 
  ggpubr::stat_cor()
p 
p + xlim(c(0.8, 1))

# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #


# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #
