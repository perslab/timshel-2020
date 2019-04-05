############### SYNOPSIS ###################
# EWCE generate cell-type data

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


### Runtime full MACA (n=10 cores) ==> ~15 min

# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

wd <- "/projects/timshel/sc-genetics/sc-genetics/src/nb-magma_celltyping"
setwd(wd)


library(tidyverse)
library(MAGMA.Celltyping) # also loads EWCE package

source("magma_celltyping_modified_functions.R")

# ======================================================================= #
# ============================  LOAD DATA  =============================== #
# ======================================================================= #

file.data <- "/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain-agg_L5.csv.gz"
df.data_raw <- read_csv(file.data) %>% rename(gene=X1)
# ^^ duplicates: ‘3110039M20Rik’, ‘4930556M19Rik’, ‘Apoc2’, ‘Atn1’, ‘C1s2’, ‘Ccdc142’, ‘Ccl19’, ‘Ccl21a’, ‘Ccl21b’, ‘Ccl21c’, ‘Ccl27a’, ‘Cd37’, ‘D130017N08Rik’, ‘Dancr’, ‘Fam205a2’, ‘Fbxw14’, ‘Flg’, ‘Gbp6’, ‘Gm15853’, ‘Gm16701’, ‘Gm2464’, ‘Gm3286’, ‘Hist2h2bb’, ‘Il11ra2’, ‘Itgam’, ‘l7Rn6’, ‘Ltbp4’, ‘Map2k7’, ‘Nova2’, ‘Ntn5’, ‘Olfr108’, ‘Olfr126’, ‘Olfr1284’, ‘Olfr1309’, ‘Olfr1316’, ‘Olfr1366’, ‘Olfr1396’, ‘Olfr1496’, ‘Olfr170’, ‘Olfr730’, ‘Olfr790’, ‘Olfr809’, ‘Pcdha11’, ‘Pcdhga8’, ‘Pik3c2g’, ‘Rp1’, ‘Schip1’, ‘Sgsm3’, ‘Smim20’, ‘Syngr4’, ‘Tead2’, ‘Tgfb1i1’, ‘Tulp2’, ‘U2af1l4’, ‘Umad1’, ‘Zfand4’ 
df.data_raw <- df.data_raw %>% distinct(gene, .keep_all = TRUE) # remove duplicate gene names
df.data_raw <- df.data_raw %>% column_to_rownames(var="gene")
head(df.data_raw)

# ======================================================================= #
# ============================  RUN  =============================== #
# ======================================================================= #

calculate.specificity.for.level <- function(df.mean_exp){
  normalised_meanExp = t(t(df.mean_exp)*(1/colSums(df.mean_exp)))
  specificity = normalised_meanExp/(apply(normalised_meanExp,1,sum)+0.000000000001)
  return(specificity)
}

ctd = list(list(
  "mean_exp"=as.matrix(df.data_raw), # matrix, genes x cells. with rownames and columnnames
  "annot"=colnames(df.data_raw),
  "specificity"=calculate.specificity.for.level(as.matrix(df.data_raw))
  ))

ctd
### STRUCTURE OF CTD FILES
# list[<data_set>] --> list of 3
# list[<data_set>]["mean_exp"] # data.frame. genes x cell_types
# list[<data_set>]["specificity"] # matrix genes x cell_types
# list[<data_set>]["annot"] # character. cell labels
# x <- ctd[[1]][["mean_exp"]]
# x <- ctd[[1]][["specificity"]]
# x <- ctd[[1]][["annot"]]

save(ctd,file="CellTypeData_mousebrain.rda")


### Load
# load("CellTypeData_mousebrain.rda") # loads ctd file
# str(ctd)

