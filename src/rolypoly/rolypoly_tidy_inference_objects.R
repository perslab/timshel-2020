############### SYNOPSIS ###################
# Analyzing multiple RolyPoly run

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

############################################



# ======================================================================================================= #
#=============================================== USAGE ================================================== #
# ======================================================================================================= #



# ======================================================================================================= #
# ============================================ OptParse ================================================= #
# ======================================================================================================= #

library(optparse)
option_list <- list( 
  make_option("--input_file", type="character", default=NULL,
              help = "File path to RData object file")
)
opt <- parse_args(OptionParser(option_list=option_list))

file.rp <- opt$input_file

# file.rp <- "out.rolypoly_objs.body_BMI_Locke2015.squared_tss_10kb.final.RData"
#file.rp <- "out.rolypoly_objs-v1/out.rolypoly_objs.disease_CARDIOVASCULAR.squared_tss_10kb.final.RData"
# file.rp <- "/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v3.univariate-nboot1000/rolypoly_objs.body_BMI_Yengo2018.tss.10kb.pos_only.all_genes.inference.RData"

# name.gwas <- stringr::str_split(basename(file.rp), "\\.")[[1]][2] # --> name.gwas = "body_BMI_Locke2015"
# print(name.gwas)

### Define output file
file.out.inference <- basename(file.rp) # writes the same name to the current working directory.

# ======================================================================================================= #
# ==========================================       CONFIG      =========================================== #
# ======================================================================================================= #

# wd <- "/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/"
# setwd(wd)

# ======================================================================================================= #
# ==========================================       SETUP      =========================================== #
# ======================================================================================================= #

library(tidyverse)

dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library


# ======================================================================================================= #
# =========================================== LOAD ROLY POLY ============================================ #
# ======================================================================================================= #

print("Loading data...")
tmp.loadtime <- system.time(load(file.rp)) # list.rp_inference, list.run_parameters
print(tmp.loadtime)


# ======================================================================================================= #
# ============================================= CLEAN OBJECT ============================================ #
# ======================================================================================================= #



### Clean 
print("Cleaning objects...")
for (name.expr.data in names(list.rp_inference)) {
  for (name.anno in names(list.rp_inference[[name.expr.data]])) {
    list.rp_inference[[name.expr.data]][[name.anno]][["rp"]][["full_results"]][["model"]] <- NULL # delete "lm model". Takes up 500 MB
  }
}

### Save
print(sprintf("Saving object: %s", file.out.inference))
tmp.runtime <- system.time(save(list.rp_inference, list.run_parameters, file=file.out.inference))
print(tmp.runtime)

print("SCRIPT DONE!")

# ======================================================================================================= #
# ============================================ SIMPLE STATS ============================================= #
# ======================================================================================================= #

### new data structure
# lapply(lapply(list.rp_wrapper, '[[', 'runtime'), '[[', 'elapsed')

