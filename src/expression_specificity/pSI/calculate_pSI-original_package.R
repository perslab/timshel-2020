
wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/nb-sem/"
setwd(wd)

### Install
# install.packages("pSI")

library(pSI)
library(tidyverse)

# ======================= READ DATA ======================= #
file.in <- "/projects/timshel/sc-genetics/sc-genetics/data/expression/arc_lira/arc_lira.celltype_expr.avg_expr.csv.gz"
df.avg_expr <- read_csv(file.in) %>% column_to_rownames("gene") %>% as.data.frame()
# df.avg_expr <- df.avg_expr %>% column_to_rownames("gene") %>% as.data.frame()

# ======================= CALC pSI/SI ======================= #
# pSI.in data frame with expresion values for genes in rows, and samples or cell types in columns (at this point replicate arrays have been averaged, so one column per cell type)
#SI logical. option to output SI value instead of default pSI value
df.pSI <- specificity.index(pSI.in=df.avg_expr, SI = FALSE)
df.SI <- specificity.index(pSI.in=df.avg_expr, SI = T)
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