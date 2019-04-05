
wd <- "/projects/timshel/sc-genetics/sc-genetics/src/nb-sem/"
setwd(wd)

### Install
# install.packages("pSI")

# library(pSI)
library(tidyverse)

# ======================= READ DATA ======================= #
#file.in <- "/projects/timshel/sc-genetics/sc-genetics/data/expression/mousebrain/mousebrain_skene.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz" # mousebrain
#file.in <- "/projects/timshel/sc-genetics/sc-genetics/data/expression/maca/maca.per_tissue_celltype.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz" # maca
file.in <- "/projects/timshel/sc-genetics/sc-genetics/data/expression/tabula_muris/tabula_muris.tissue_celltype.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz" # tabula_muris
df.avg_expr <- read_csv(file.in) %>% column_to_rownames("gene") %>% as.data.frame()

# ======================= CALC pSI/SI ======================= #
source("/projects/timshel/sc-genetics/sc-genetics/src/nb-sem/specificity.index.timshel.R") # loads 'specificity.index.timshel' function
system.time(list.si <- specificity.index.timshel(pSI.in=df.avg_expr, p_max = 1, e_min=1e-5))
# list(SI=df.SI, pSI=df.pSI) # named list
# save.image("psi.mousebrain.pmax_1.e_min_1e-5.RData")
#save.image("psi.maca_tissue_cell_type.pmax_1.e_min_1e-5.RData")
save.image("psi.tabula_muris_tissuecell_type.pmax_1.e_min_1e-5.RData")

print("DONE")


# ======================= RUNTIMES ======================= #
### RUNTIMES MOUSEBRAIN (14945 genes x 265 cell-types) (mousebrain_skene.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz)
# default	~0.7 h
# p_max_0.5	~ 1 h
# p_max_1	~ 2 h
# pmax_0.1_e_min_0	~5 h
# Mousebrain p_max = 1, e_min=1e-5: 25 h
# MACA  p_max = 1, e_min=1e-5: 13.7h

# ======================= CALC pSI/SI - original function ======================= #
# pSI.in data frame with expresion values for genes in rows, and samples or cell types in columns (at this point replicate arrays have been averaged, so one column per cell type)
#SI logical. option to output SI value instead of default pSI value
# system.time(df.pSI <- specificity.index(pSI.in=df.avg_expr, SI = FALSE, p_max = 1, e_min=0))
# save.image("psi.mousebrain.pmax_0.1_e_min_0.RData")
# system.time(df.SI <- specificity.index(pSI.in=df.avg_expr, SI = T, p_max = 1, e_min=0))
# save.image("psi.mousebrain.pmax_0.1_e_min_0.RData")



### Arguments
# specificity.index(pSI.in, pSI.in.filter, bts = 50, p_max = 0.1, e_min = 0.3, hist = FALSE, SI = FALSE)
#bts	: numeric. number of distributions to average for permutation testing
#p_max	:numeric. maximum pvalue to be calculated ----> set pvalue to NA if pvalue is larger than p_max.

# ======================= Write files ======================= #

# df.pSI %>% rownames_to_column(var="gene") %>% write_csv("arc_lira.pSI.csv.gz")
# df.SI %>% rownames_to_column(var="gene") %>% write_csv("arc_lira.SI.csv.gz")


# ======================= pSI DOCS ======================= #
# pSI-package	Specificity Index Statistic
# candidate.genes	Candidate Gene Lists
# candidate.overlap	Candidate Gene Overlap
# dataset.s1	Supplementary Data Set 1
# fisher.iteration	Fisher's Exact Test Across All Cell Types & pSI Thresholds
# pSI	Specificity Index Statistic
# pSI.count	Convert pSI output to gene count list
# pSI.list	Convert pSI output to gene list
# sample.data	Sample Input & Output pSI Data Sets
# specificity.index	Specificity Index Statistic

### pSI.count This functions counts number of genes specific to each sample type
# pSI.count(pSIs, write.csv = FALSE) | pSIs:	data frame output from specificity.index function with the number of columns equal to the number of samples and genes as rows.

### pSI.list returns list consisting of 6 data frames, one for each pSI threshold.
# pSI.list(pSIs, write.csv = TRUE) # | pSIs:	data frame output from specificity.index function with the number of columns equal to the number of samples and genes as rows.

### fisher.iteration will test a candidate gene list for overrepresenation in the various cell type/pSI threshold combinations produced by the specificty.index function. 
# fisher.iteration(pSIs, candidate.genes, background = "data.set", p.adjust = TRUE)

# ======================= Install (OLD) ======================= #
### REF:http://genetics.wustl.edu/jdlab/psi_package/
# wget http://genetics.wustl.edu/jdlab/files/2014/01/pSI_1.1.tar_.gz # pSI_1.1.tar
# install.packages('pSI_1.1.tar_.gz', repos = NULL, type="source")