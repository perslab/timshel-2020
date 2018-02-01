# Title: Loads all functions in the directory into the R workspace. This makes it convenient to load all functions in one go.
# By sourcing this script, we ensure that all functions in this folder are accessible in all files. That is, functions are SHARED across files.

# Author: Pascal N. Timshel, Pers Lab
# Date: December 2017



####################################################################################################
########################################### FUNCTIONS ##############################################
####################################################################################################


# COPY FROM: https://github.com/molgenis/molgenis-pipelines/wiki/How-to-source-another_file.R-from-within-your-R-script
# How to source another_file.R from within your R script, mdijkstra edited this page on Sep 30, 2014

LocationOfThisScript = function() # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
  this.file = NULL
  # This file may be 'sourced'
  for (i in -(1:sys.nframe())) {
    if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
  }
  
  if (!is.null(this.file)) return(dirname(this.file))
  
  # But it may also be called from the command line
  cmd.args = commandArgs(trailingOnly = FALSE)
  cmd.args.trailing = commandArgs(trailingOnly = TRUE)
  cmd.args = cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
  res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)
  
  # If multiple --file arguments are given, R uses the last one
  res = tail(res[res != ""], 1)
  if (0 < length(res)) return(dirname(res))
  
  # Both are not the case. Maybe we are in an R GUI?
  return(NULL)
}

####################################################################################################
####################################### CONSTANTS ##################################################
####################################################################################################

SCRIPT_NAME <- "load_functions.R" # this script name

####################################################################################################
########################################## MAIN ####################################################
####################################################################################################

print(sprintf("Started loading functions via %s...", SCRIPT_NAME))

dir.lib_functions = LocationOfThisScript() # get the directory where this file lives. 
print(sprintf("Lib function dir: %s", dir.lib_functions))

files.lib_functions <- list.files(dir.lib_functions, pattern="*.R") # get all .R files in the current directory as a vector. Returns basename (e.g. "plot_seurat.R") and not full path

for (file2source in files.lib_functions) {
  if (file2source == SCRIPT_NAME) { # don't source this file. We need this to avoid ending up in a infinity 'file sourcing' loop!
    next
  }
  tmp.filepath <- file.path(dir.lib_functions, file2source) # creating absolute path
  source(tmp.filepath)
  print(sprintf("Sourced file: %s", tmp.filepath))
}

print("Done loading functions")


####################################################################################################
######################################## CLEAN UP ##################################################
####################################################################################################

rm(SCRIPT_NAME, tmp.filepath, LocationOfThisScript) # we don't want these variables to clutter our workspace.


