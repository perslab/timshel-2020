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
# ============================ Correlation heatmap =============================== #
# ======================================================================= #

### corrr: exploring correlations
### Code: https://github.com/drsimonj/corrr/
### Blog: https://drsimonj.svbtle.com/exploring-correlations-in-r-with-corrr
library(corrr)

x <- as_cordf(cormat.es) %>% # instead of calling correlate(df.es %>% select(-gene))
  focus(matches("^MSN"), mirror=T)
x
x %>% rplot()

dim(cormat.es)

# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #




# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #




# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #


