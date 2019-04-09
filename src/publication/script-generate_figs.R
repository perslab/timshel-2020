############### SYNOPSIS ###################
### Generate publication figs by sourcing all fig-(.*).R files in src/publication

# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #
library(here)
source(here("src/lib/load_functions.R")) # load sc-genetics library

setwd(here("src/publication"))
files_source <- list.files(pattern="fig-(.*).R") # returns filenames (not full filepaths)
# ======================================================================= #
# ================================= RUN ================================= #
# ======================================================================= #
for (file_source in files_source) {
  print(sprintf("Sourcing file: %s...", file_source))
  source(file_source)
}

