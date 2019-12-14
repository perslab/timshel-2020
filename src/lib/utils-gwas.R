################ SYNOPSIS ###################
# Functions to rename GWAS trait names


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
# ===================== DATA GWAS_ID META-ANALYSIS ====================== #
# ======================================================================= #



# ======================================================================= #
# =========================== DATA RENAMING ============================= #
# ======================================================================= #

file.gwas_rename <- here("data/gwas/gwas_trait_naming.csv")
### HEADER
# trait_name_full,trait_name_abrv,author,year,gwas_id
# Multiple Sclerosis,MS,Patsopoulos,2011,MS_Patsopoulos2011
# Alzheimer's Disease,AD,Lambert,2013,AD_Lambert2013

# ======================================================================= #
# ========================= FUNCTION: RENAME ============================= #
# ======================================================================= #


utils.get_gwas_ids_for_meta_analysis <- function() {
  # N_GWAS = 39
  # GWAS represent represent the lastest GWAS for each trait.
  # Traits are selected to be as independent/non-redundant as possible.
  # LAST MODIFIED: 04.05.2019
  file.gwas_ids_meta_analysis <- here("data/gwas/gwas_traits.csv")
  gwas_ids.meta_analysis <- suppressMessages(read_csv(file.gwas_ids_meta_analysis, col_names=F)) %>% pull(X1)
  print(sprintf("Returning n=%s gwas_ids", length(gwas_ids.meta_analysis)))
  return(gwas_ids.meta_analysis)
}
  

# ======================================================================= #
# ======================== FUNCTION: RENAME GWAS ========================= #
# ======================================================================= #

utils.rename_gwas <- function(gwas_ids, style, check_all_matches=F, return_as_df=F) {
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
  
  ### DEV
  # style <- "fullname_author_year"
  # gwas_ids <- c("CROHNS_Jostins2012", "LIPIDS_HDL_Teslovich2010", "asdasd") # must fail
  # gwas_ids <- c("CROHNS_Jostins2012", "LIPIDS_HDL_Teslovich2010") # must pass
  # return_as_df <- FALSE
  # check_all_matches <- FALSE

  ### Check input
  allowed_styles <- c("fullname_author_year",
                      "fullname",
                      "abrv_author_year",
                      "abrv_year",
                      "abrv" # OBS: abrv may not be unique!
                      )
  if (!style %in% allowed_styles) {
    stop("Got wrong style. Allowed styles are: [%s]", paste(allowed_styles, collapse=", "))
  }
  if (style == "abrv") {
    print("Warning: trait name abreviations are not unique.")
  }
  ### Read data
  df.gwas_database <- suppressMessages(read_csv(file.gwas_rename))
  if (any(duplicated(df.gwas_database$gwas_id))) {
    stop("internal renaming data frame contains duplicates. Fix the .csv file")
  }
  ### Check that all gwas_ids are present in "data base"
  bool.found <- gwas_ids %in% df.gwas_database$gwas_id
  if (check_all_matches && !all(bool.found)) {
    gwas_ids_not_found <- gwas_ids[!bool.found]
    stop(sprintf("Some gwas_ids not found: [%s]", paste(gwas_ids_not_found, collapse=",")))
  }
  ### Format output conditional on style
  df.gwas_database.select <- df.gwas_database %>% filter(gwas_id %in% gwas_ids)
  if (style=="fullname_author_year") {
    df.gwas_out_fmt <- df.gwas_database.select %>% 
      mutate(fmt=paste0(trait_name_full, " (", author, ", ", year, ")"))
  } else if (style=="fullname") {
    df.gwas_out_fmt <- df.gwas_database.select %>% 
      mutate(fmt=paste0(trait_name_full))
  } else if (style=="abrv_author_year") {
      df.gwas_out_fmt <- df.gwas_database.select %>% 
        mutate(fmt=paste0(trait_name_abrv, " (", author, ", ", year, ")"))
  } else if (style=="abrv_year") {
    df.gwas_out_fmt <- df.gwas_database.select %>% 
      mutate(fmt=paste0(trait_name_abrv, " (", year, ")"))
  } else if (style=="abrv") {
    df.gwas_out_fmt <- df.gwas_database.select %>% 
      mutate(fmt=paste0(trait_name_abrv))
  } else {
    stop("Internal function error")
  }
  ### Join renamed with input
  df.all <- tibble(gwas_id=gwas_ids) # this ensures that the return value has the same length as the input
  df.all <- df.all %>% left_join(df.gwas_out_fmt, by="gwas_id")
  ### Add fmt column for gwas_id non-matches
  ### df.all <- df.all %>% mutate(fmt=replace(fmt, is.na(fmt), gwas_id)) # GIVES WARNING ---> Warning message: In x[list] <- values :number of items to replace is not a multiple of replacement length
  df.all <- df.all %>% mutate(fmt=if_else(is.na(fmt), gwas_id, fmt)) 
  # ^ gwas_id non-matches have NA values in the 'fmt' column
  # ^ we replace na with the input gwas_id name
  if (return_as_df) {
    gwas_out_fmt <- df.all
  } else(
    gwas_out_fmt <- df.all %>% pull(fmt) # returns vector
  )
  return(gwas_out_fmt)
}

### EXAMPLE CALL
# style <- "fullname_author_year"
# gwas_ids <- c("CROHNS_Jostins2012", "LIPIDS_HDL_Teslovich2010") # must pass
# gwas_ids.fmt <- utils.rename_gwas(gwas_ids, style, return_as_df=F)
# gwas_ids.fmt




# ======================================================================= #
# =============================== RENAME FUNCTION ================================= #
# ======================================================================= #
# THIS FUNCTION WORKS


utils.rename_gwas.drop_no_matches <- function(gwas_ids, style, check_all_matches=T, return_as_df=F) {
  ### INPUT
  # gwas_ids:           a vector of gwas_ids's
  # style:              style.
  # check_all_matches:  if true, the function will raise an exception if the input gwas_ids are not all found in the database.
  #                     if false, non-matches will be dropped silently.
  ### OUTPUT
  # a vector of formatted gwas names. The length is equal to the number of gwas_id that matched the database.
  
  ### DEV
  # style <- "fullname_author_year"
  # ## gwas_ids <- c("CROHNS_Jostins2012", "LIPIDS_HDL_Teslovich2010", "asdasd") # must fail
  # gwas_ids <- c("CROHNS_Jostins2012", "LIPIDS_HDL_Teslovich2010") # must pass
  # return_as_df <- FALSE
  # 
  ### Check input
  allowed_styles <- c("fullname_author_year",
                      "fullname",
                      "abrv_author_year",
                      "abrv_year",
                      "abrv" # OBS: abrv may not be unique!
  )
  if (!style %in% allowed_styles) {
    stop("Got wrong style. Allowed styles are: [%s]", paste(allowed_styles, collapse=", "))
  }
  if (style == "abrv") {
    print("Warning: trait name abreviations are not unique.")
  }
  ### Read data
  df.gwas_database <- read_tsv(.data.rename_gwas)
  if (any(duplicated(df.gwas_database$gwas_id))) {
    stop("internal renaming data frame contains duplicates. Fix the .csv file")
  }
  ### Check that all gwas_ids are present in "data base"
  bool.found <- gwas_ids %in% df.gwas_database$gwas_id
  if (check_all_matches && !all(bool.found)) {
    gwas_ids_not_found <- gwas_ids[!bool.found]
    stop(sprintf("Some gwas_ids not found: [%s]", paste(gwas_ids_not_found, collapse=",")))
  }
  ### Format output conditional on style
  df.gwas_database.select <- df.gwas_database %>% filter(gwas_id %in% gwas_ids)
  if (style=="fullname_author_year") {
    df.gwas_out_fmt <- df.gwas_database.select %>% 
      mutate(fmt=paste0(trait_name_full, " (", author, ", ", year, ")"))
  } else if (style=="fullname") {
    df.gwas_out_fmt <- df.gwas_database.select %>% 
      mutate(fmt=paste0(trait_name_full))
  } else if (style=="abrv_author_year") {
    df.gwas_out_fmt <- df.gwas_database.select %>% 
      mutate(fmt=paste0(trait_abrv, " (", author, ", ", year, ")"))
  } else if (style=="abrv_year") {
    df.gwas_out_fmt <- df.gwas_database.select %>% 
      mutate(fmt=paste0(trait_abrv, " (", year, ")"))
  } else if (style=="abrv") {
    df.gwas_out_fmt <- df.gwas_database.select %>% 
      mutate(fmt=paste0(trait_abrv))
  } else {
    stop("Internal function error")
  }
  if (return_as_df) {
    gwas_out_fmt <- df.gwas_out_fmt
  } else(
    gwas_out_fmt <- df.gwas_out_fmt %>% pull(fmt) # returns vector
  )
  return(gwas_out_fmt)
}
