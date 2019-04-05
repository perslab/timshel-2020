


wd <- "/projects/timshel/sc-genetics/sc-genetics/src/nb-sem/"
setwd(wd)

library(tidyverse)

#load("psi.maca_tissue_cell_type.pmax_1.e_min_1e-5.RData")
load("psi.tabula_muris_tissuecell_type.pmax_1.e_min_1e-5.RData")
names(list.si)
df.SI <- list.si[["SI"]]
df.pSI <- list.si[["pSI"]]

i <- 4
### cell-types
qplot(unlist(df.SI[i,]), unlist(df.pSI[i,]))

### genes
qplot(unlist(df.SI[,i]), unlist(df.pSI[,i]))
qplot(unlist(df.SI[,i]), -log10(unlist(df.pSI[,i])+1))
qplot(rank(unlist(df.SI[,i])), rank(unlist(df.pSI[,i])))

cor.test(unlist(df.SI[,i]), unlist(df.pSI[,i]), method="spearman") # ---> perfect correlation



