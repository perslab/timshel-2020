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

# Rscript /projects/timshel/sc-genetics/sc-genetics/src/RP-meta/rolypoly_tidy_gwas_linked_objects.R --input_file /projects/timshel/sc-genetics/sc-genetics/out/out.rolypoly_objs-v2/rolypoly_objs.blood_EOSINOPHIL_COUNT.tss.10kb.square.all_genes.gwas_linked.RData --dir_out /scratch/sc-genetics/tmp_gwas_linked"

# ======================================================================================================= #
# ============================================ OptParse ================================================= #
# ======================================================================================================= #

library(optparse)
option_list <- list( 
  make_option("--input_file", type="character", default=NULL,
              help = "File path to RData object file"),
  make_option("--dir_out", type="character", default=NULL,
              help = "File path output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

file.rp <- opt$input_file
dir_out <- opt$dir_out

dir.create(dir_out,showWarnings=F, recursive=T) # make new directory. cmd does nothing if the directory already exists.

# file.rp <- "out.rolypoly_objs.body_BMI_Locke2015.squared_tss_10kb.final.RData"
#file.rp <- "out.rolypoly_objs-v1/out.rolypoly_objs.disease_CARDIOVASCULAR.squared_tss_10kb.final.RData"
# file.rp <- "/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v3.univariate-nboot1000/rolypoly_objs.body_BMI_Yengo2018.tss.10kb.pos_only.all_genes.inference.RData"

# name.gwas <- stringr::str_split(basename(file.rp), "\\.")[[1]][2] # --> name.gwas = "body_BMI_Locke2015"
# print(name.gwas)

### Define output file
# file.out.inference <- basename(file.rp) # writes the same name to the current working directory.
file.out.gwas_linked <- file.path(dir_out, basename(file.rp)) # writes the same name to the current working directory.

# ======================================================================================================= #
# ==========================================       CONFIG      =========================================== #
# ======================================================================================================= #

# wd <- "/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/"
# setwd(wd)

# ======================================================================================================= #
# ==========================================       SETUP      =========================================== #
# ======================================================================================================= #

library(tidyverse)

# dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib"
# source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library

# ======================================================================================================= #
# =========================================== LOAD ROLY POLY ============================================ #
# ======================================================================================================= #

print("Loading data...")
# tmp.loadtime <- system.time(load(file.rp)) # list.rp_inference, list.run_parameters
tmp.loadtime <- system.time(load(file.rp)) # list.gwas_linked, list.run_parameters
print(tmp.loadtime)


# ======================================================================================================= #
# ============================================= CLEAN OBJECT ============================================ #
# ======================================================================================================= #

# list.gwas_linked.orig <- list.gwas_linked

### Clean 
print("Cleaning objects...")
list.gwas_linked[["rp"]] <- list.gwas_linked[["rp"]][[1]]
list.gwas_linked[["runtime"]] <- list.gwas_linked[["runtime"]][[1]]

### Save
print(sprintf("Saving object: %s", file.out.gwas_linked))
tmp.runtime <- system.time(save(list.gwas_linked, list.run_parameters, file=file.out.gwas_linked))
print(tmp.runtime)

print("SCRIPT DONE!")
