############### SYNOPSIS ###################
### Load all publication library functions

# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #
library(here)
setwd(here("src/publication"))
files_source <- list.files(pattern="lib-(.*).R") # returns filenames (not full filepaths)
SCRIPT_NAME <- "lib-load_functions.R" # *MUST* be kept up-to-date with the name of this script
# ======================================================================= #
# ================================= RUN ================================= #
# ======================================================================= #
for (file_source in files_source) {
  if (file_source == SCRIPT_NAME) { # don't source this file. We need this to avoid ending up in a infinity 'file sourcing' loop!
    next
  }
  source(file_source)
  print(sprintf("Sourced file: %s", file_source))
}

### Alternative
# R.utils::sourceDirectory(path=".", pattern="lib-(.*).R", modifiedOnly=TRUE) # REF: https://stat.ethz.ch/pipermail/r-help/2008-September/173606.html