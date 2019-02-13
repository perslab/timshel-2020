############### SYNOPSIS ###################
# Explore and plot espression specificity correlations
# Dendrogram
# Correlation heatmap

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


dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib/"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library

setwd(here("src/expression_specificity"))


# ======================================================================= #
# ============================ PARAMETERS =============================== #
# ======================================================================= #

dataset_prefix <- "mousebrain_all"

filter.celltypes <- c("TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5","MEGLU10","TEGLU17","MEGLU11","TEGLU4","DEGLU4","TEINH12") # BMI_UKBB_Loh2018 FDR sign.

# ======================================================================= #
# ============================ LOAD DATA =============================== #
# ======================================================================= #

file.es <- here(sprintf("data/genes_cell_type_specific/%s.mean.csv.gz", dataset_prefix))
df.es <- read_csv(file.es) # genes x cell-types

# ======================================================================= #
# ============================ Calculate correlation =============================== #
# ======================================================================= #

cormat.es <- cor(df.es %>% select(-gene), method="pearson") 

# ======================================================================= #
# ============================ Calculate dendrogram =============================== #
# ======================================================================= #

### Compute distances and hierarchical clustering
# dd.eucledian <- dist(t(df.es %>% select(-gene)), method = "euclidean") # computes distances between the rows of a data matrix
dd.corr <- as.dist(1-cormat.es)
hc <- hclust(dd.corr, method = "ward.D2") # Hierarchical clustering 
dend <- as.dendrogram(hc) # Turn the object into a dendrogram.
plot(dend)

# ======================================================================= #
# ============================ [corrplot] Correlation plots =============================== #
# ======================================================================= #

### https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
library(corrplot)

pdf("plot.corrplot.mixed.pdf", width=20, height=20)
corrplot.mixed(cormat.es)
dev.off()

pdf("plot.corrplot.hclust.pdf", width=20, height=20)
corrplot(cormat.es, order = "hclust", addrect = 10)
dev.off()

### selected FDR using corrplot
tmp.cor <- cor(df.es %>% select(one_of(filter.celltypes), -gene), method="pearson")  # one_of(): variables in character vector. You need this when mixing unquoted names and character vector
pdf("plot.corrplot.mixed_hclust_fdr_only.pdf", width=10, height=10)
corrplot.mixed(tmp.cor, order = "hclust")
dev.off()


### selected FDR using corrplot --> GIVES WRONG RESULT BECAUSE focus() does not return diagonal matrix
# cormat.es.seleced <- as_cordf(cormat.es) %>% focus(filter.celltypes, mirror=T) %>% as_matrix()
# corrplot.mixed(cormat.es.seleced, is.corr=FALSE) # is.corr=FALSE convert to correlation matrix


# ======================================================================= #
# ============================ [corrr] Correlation plots =============================== #
# ======================================================================= #

### corrr: exploring correlations
### Code: https://github.com/drsimonj/corrr/
### Blog: https://drsimonj.svbtle.com/exploring-correlations-in-r-with-corrr
library(corrr)

df.es_cordf <- as_cordf(cormat.es) # instead of calling correlate(df.es %>% select(-gene))

### selecting cell-types
# df.es_cordf %>% focus(matches("^MSN"), mirror=T)
# df.es_cordf %>% focus(filter.celltypes, mirror=T) %>% rearrange(method="HC") %>% shave() %>% rplot(print_cor=T) # shave()
df.es_cordf %>% focus(filter.celltypes, mirror=T) %>% rearrange(method="HC") %>% rplot(print_cor=T)
# *OBS*: focus does NOT return a diagonal matrix. That is, it is not 'ordered'.

# df.es_cordf %>% focus(filter.celltypes, mirror=T) %>% select(sort(current_vars())) # ... trying to order the colnames but then rowname is included

### INDEX
# as_cordf	Coerce lists and matrices to correlation data frames
# as_matrix	Convert a correlation data frame to matrix format
# correlate	Correlation Data Frame
# fashion	Fashion a correlation data frame for printing.
# first_col	Add a first column to a data.frame
# focus	Focus on section of a correlation data frame.
# focus_	Focus on section of a correlation data frame.
# focus_if	Conditionally focus correlation data frame
# network_plot	Network plot of a correlation data frame
# pair_n	Number of pairwise complete cases.
# rearrange	Re-arrange a correlation data frame
# rplot	Plot a correlation data frame.
# shave	Shave off upper/lower triangle.
# stretch	Stretch correlation data frame into long format.

# ======================================================================= #
# ============================ [ggcor] Correlation plots =============================== #
# ======================================================================= #
library(GGally)
# ggcorr (used for lira lab stats):
# part of GGally
# intro: https://briatte.github.io/ggcorr/ 
# Repo: https://github.com/briatte/ggcorr

# ggcorr(data = NULL, cor_matrix = cormat.es)
# ggcorr(df.lab, geom = "circle", hjust = 1, size = 2, layout.exp = 5) # label = TRUE
# ggcorr(df.lab, label = TRUE, label_size = 1.5, label_round = 2, label_alpha = TRUE, hjust = 1, size = 2, layout.exp = 5)



# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #


