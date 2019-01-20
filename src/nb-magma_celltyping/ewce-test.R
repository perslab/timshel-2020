############### SYNOPSIS ###################
# EWCE generate cell-type data

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/nb-magma_celltyping"
setwd(wd)

library(Seurat) # only needed if fiddling with Seurat objects
library(MAGMA.Celltyping) # also loads EWCE package

source("magma_celltyping_modified_functions.R")

# ======================================================================= #
# ==========================  EXAMPLE DATA  ============================= #
# ======================================================================= #

# Load the single cell expression data
data("cortex_mrna")

str(cortex_mrna) # --> list of 2. 
# "exp": matrix, genes x cells
# "annot": data.frame, cells x meta_data

expData = cortex_mrna$exp # genes x cells
l1=cortex_mrna$annot$level1class # n=cells, cell_type labels
l2=cortex_mrna$annot$level2class # n=cells, cell_type labels
annotLevels = list(l1=l1,l2=l2)


### Runs in parallel
### https://github.com/NathanSkene/EWCE/blob/master/R/generate.celltype.data.r
# fNames_ALLCELLS = generate.celltype.data(exp=expData,annotLevels=annotLevels,groupName="allKImouse") # returns filename for saved celltype_data file
# ^^: saves file "CellTypeData_<groupName>.rda" using 'save(ctd, ..)' command

### Load
# load("CellTypeData_allKImouse.rda") # loads ctd file
# str(ctd)

### STRUCTURE OF CTD FILES
# list[<data_set>] --> list of 3
# list[<data_set>]["mean_exp"] # data.frame. genes x cell_types
# list[<data_set>]["specificity"] # matrix genes x cell_types
# list[<data_set>]["annot"] # character. cell labels
# x <- ctd[[1]][["mean_exp"]]
# x <- ctd[[1]][["specificity"]]
# x <- ctd[[1]][["annot"]]


### Experiment
# str(cortex_mrna$annot$level1class)
# str(cortex_mrna$annot$level2class)
# annotLevels = list(level1class=cortex_mrna$annot$level1class,level2class=cortex_mrna$annot$level2class)
# fNames_CortexOnly = generate.celltype.data(exp=exp_CortexOnly_DROPPED,annotLevels=annotLevels,groupName="kiCortexOnly")
# print(fNames_CortexOnly)
# fNames_CortexOnly = filter.genes.without.1to1.homolog(fNames_CortexOnly)
# print(fNames_CortexOnly)


