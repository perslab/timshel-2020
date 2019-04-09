############### SYNOPSIS ###################
# AIM: Explore and plot espression specificity correlations
# *OBS*: THERES PLOTS WERE LEFT OUT OF THE PUBLICATION FIGURE SCRIPT BECAUSE THEY WERE LESS GOOD.
# *OBS*: SEE "fig-annotation_clusting" for the publication version


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


source(here("src/lib/load_functions.R")) # load sc-genetics library

setwd(here("src/expression_specificity"))

# ======================================================================= #
# ============================ COMPUTE CORMAT =============================== #
# ======================================================================= #

# .... SEE "fig-annotation_clusting" for how to do this


# ======================================================================= #
# ============================ [corrr] Correlation plots =============================== #
# ======================================================================= #

### corrr: exploring correlations
### Code: https://github.com/drsimonj/corrr/
### Blog: https://drsimonj.svbtle.com/exploring-correlations-in-r-with-corrr
library(corrr)
# *OBS*: focus() does NOT return a diagonal matrix. That is, it is not 'ordered'.

df.es_cordf <- as_cordf(cormat.es) # instead of calling correlate(df.es %>% select(-gene))

### selecting cell-types
# df.es_cordf %>% focus(matches("^MSN"), mirror=T)
# df.es_cordf %>% focus(filter.annotations, mirror=T) %>% rearrange(method="HC") %>% shave() %>% rplot(print_cor=T) # shave() looks weird when the matrix is not diagonal
df.es_cordf %>% focus(filter.annotations, mirror=T) %>% rearrange(method="HC") %>% rplot(print_cor=T)
ggsave(sprintf("plot.%s.corrr_hc.pdf",dataset_prefix), w=10, h=10)

# df.es_cordf %>% focus(filter.annotations, mirror=T) %>% select(sort(current_vars())) # ... trying to order the colnames but then rowname is included, so it becomes to much of a mess to solve

### corrr INDEX
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
# library(GGally)
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


