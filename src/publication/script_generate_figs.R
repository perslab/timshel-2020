############### SYNOPSIS ###################
### Generate publication figs

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

setwd(here("src/publication"))


# ======================================================================= #
# ================================= RUN ================================= #
# ======================================================================= #

files_fig <- list.files(pattern="fig_(.*).R")
for (file_fig in files.fig) {
  print(file_fig)
  source(file_fig)
}



