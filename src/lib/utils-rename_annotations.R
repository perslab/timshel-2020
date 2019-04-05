################ SYNOPSIS ###################
# Functions to rename dataset annotation names


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



# ======================================================================= #
# =============================== DATA ================================= #
# ======================================================================= #


# ======================================================================= #
# =============================== TM RENAME FUNCTION ================================= #
# ======================================================================= #

utils.rename_annotations.tabula_muris <- function(annotation_ids, style, check_all_matches=F, return_as_df=F) {
  ### DESCRIPTION
  # This function is an alternative version of utils.rename_gwas.drop_no_matches()
  # gwas_id's not found in the data base, will not be renamed.
  # *Advantage of this function*: this function guarantees the same order and length of input-to-output.
  ### INPUT
  # gwas_ids:      a vector of gwas_ids's
  # style:        style.
  # check_all_matches:  if true, the function will raise an exception if the input gwas_ids are not all found in the database.
  ### OUTPUT
  # a vector of formatted gwas names. 
  # the length of the vector is equal to the input gwas_ids.
  stop("NOT IMPLEMENTED YET")

}

### EXAMPLE CALL
# style <- "fullname_author_year"
# gwas_ids <- c("CROHNS_Jostins2012", "LIPIDS_HDL_Teslovich2010") # must pass
# gwas_ids.fmt <- utils.rename_gwas(gwas_ids, style, return_as_df=F)
# gwas_ids.fmt

# ======================================================================= #
# ========================= CAMPBELL RENAME FUNCTION ===================== #
# ======================================================================= #


utils.rename_annotations.campbell2015 <- function(annotation_ids, style, check_all_matches=F, return_as_df=F) {
  ### DESCRIPTION
  # This function is an alternative version of utils.rename_gwas.drop_no_matches()
  # gwas_id's not found in the data base, will not be renamed.
  # *Advantage of this function*: this function guarantees the same order and length of input-to-output.
  ### INPUT
  # gwas_ids:      a vector of gwas_ids's
  # style:        style.
  # check_all_matches:  if true, the function will raise an exception if the input gwas_ids are not all found in the database.
  ### OUTPUT
  # a vector of formatted gwas names. 
  # the length of the vector is equal to the input gwas_ids.
  stop("NOT IMPLEMENTED YET")

}

