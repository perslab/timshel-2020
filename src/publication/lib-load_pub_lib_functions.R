### DESCRIPTION
# Loads all functions in the lib_path into the R workspace. This makes it convenient to load all functions in one go.
# By sourcing this script, we ensure that all functions in this folder are accessible in all files. That is, functions are SHARED across files.

# Author: Pascal N. Timshel

# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #
library(here)

SCRIPT_NAME <- "lib-load_pub_lib_functions.R" # *MUST* be kept up-to-date with the name of this script
lib_path <- here("src/publication")
lib_pattern <- "lib-(.*).R" # *OBS*: specific pattern

# ======================================================================= #
# ================================ RUN ================================ #
# ======================================================================= #

files_source <- list.files(path=lib_path, pattern=lib_pattern) # returns filenames (not full filepaths)

for (file_source in files_source) {
  if (file_source == SCRIPT_NAME) { # don't source this file. We need this to avoid ending up in a infinity 'file sourcing' loop!
    next
  }
  file_source.filepath <- file.path(lib_path, file_source) # creating absolute path
  source(file_source.filepath)
  print(sprintf("Sourced file: %s", file_source))
}


# we don't want these variables to clutter our workspace.
rm(lib_path, lib_pattern,
   SCRIPT_NAME,
   files_source, file_source) 

# ======================================================================= #
# ======================= NEW METHOD [no print statement] ================ #
# ======================================================================= #

# R.utils::sourceDirectory(path=lib_path, pattern=lib_pattern, modifiedOnly=TRUE) # REF: https://stat.ethz.ch/pipermail/r-help/2008-September/173606.html

