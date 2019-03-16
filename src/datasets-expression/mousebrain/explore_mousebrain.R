############### SYNOPSIS ###################
# Explore Linnarson mouse brain data

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

library(tidyverse)
# library(GGally)

devtools::install_github(repo = 'hhoeflin/hdf5r', args = "--configure-vars='LIBS=-L/tools/anaconda/3-4.4.0/lib'") # --> configure: error: The version of hdf5 installed on your system is not sufficient. Please ensure that at least version 1.8.13 is installed

devtools::install_github(repo = 'hhoeflin/hdf5r', args = "--configure-vars='HDF5R_LIBS=-L/tools/anaconda/3-4.4.0/lib -L. -lhdf5_cpp -lhdf5 -lz -lm'") # ---> sh: -L.: command not found

devtools::install_github(repo = 'hhoeflin/hdf5r', args = "--configure-vars='HDF5R_LIBS=-L/tools/anaconda/3-4.4.0/lib -lhdf5 -lhdf5_hl -lz -lm'") # ---> sh: -L.: command not found
devtools::install_github(repo = 'hhoeflin/hdf5r', args = "--configure-vars='HDF5R_LIBS=-L/tools/anaconda/3-4.4.0/lib'") # ---> sh: -L.: command not found


### Fails
# devtools::install_github(repo = 'hhoeflin/hdf5r') # Update hdf5r for stability improvements
# devtools::install_github(repo = 'mojaveazure/loomR', ref = 'develop') # Update loomR for up-to-date functions

### devtools::install_github(repo = "hhoeflin/hdf5r") # FAILS
  ### ^^ Installing hdf5r from github repo FAILS! (Ygg has HDF5 version 1.8.12): "configure: error: The version of hdf5 installed on your system is not sufficient. Please ensure that at least version 1.8.13 is installed"
# devtools::install_github(repo = "mojaveazure/loomR") # WORKS.
  # ^ NOW hdf5r can be installed: trying URL 'https://cran.rstudio.com/src/contrib/hdf5r_1.0.0.tar.gz'
library(loomR)
library(Seurat)

# ======================================================================= #
# ============================ LOOMR tutorial =========================== #
# ======================================================================= #

# #1: https://satijalab.org/seurat/mca_loom.html
# #2: http://satijalab.org/loomR/loomR_tutorial.html

# ======================================================================= #
# =============================== PARAMS ================================ #
# ======================================================================= #

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/"
setwd(wd)

# ======================================================================= #
# =============================== XXXXXX ================================ #
# ======================================================================= #

# Connect to the loom file in read/write mode
lfile <- loomR::connect(filename = "pbmc.loom", mode = "r+")
lfile


# ======================================================================= #
# =============================== XXXXXX ================================ #
# ======================================================================= #


# ======================================================================= #
# =============================== XXXXXX ================================ #
# ======================================================================= #


