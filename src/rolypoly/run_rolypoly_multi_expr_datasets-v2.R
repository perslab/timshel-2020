############### SYNOPSIS ###################
# Running RolyPoly

### v2 change log: many things! Here is a selection:
# NB: v2 is still a *multiple regression* over all cell-types.
# NB: "loading pre-computed GWAS data" does NOT work (--gwas_linked_file)
# maf as covariate to model add_poly = T, n_degree = 1
# run_rolypoly.multi_attempt() function does both GWAS 'pre-computation' and per-data set calculation
# added rolypoly_add_ld_corrected_gwas_block_scores()
# change output file and variables names (to e.g. list.gwas_linked, list.rp_inference.run)


### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

############################################



# ======================================================================================================= #
#=============================================== USAGE ================================================== #
# ======================================================================================================= #

### body_WHR_Shungin2015
# time Rscript run_rolypoly_multi_expr_datasets.R --gwas_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_WHR_Shungin2015.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_WHR_Shungin2015 --run_name squared_tss_10kb --n_cores 10 |& tee log.body_WHR_Shungin2015.out.txt
### cov_EDU_YEARS
# time Rscript run_rolypoly_multi_expr_datasets.R --gwas_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/cov_EDU_YEARS.gwassumstats.rolypoly_fmt.tab.gz --gwas_name cov_EDU_YEARS --run_name squared_tss_10kb --n_cores 10 |& tee log.cov_EDU_YEARS.out.txt
### body_WHRadjBMIz
# time Rscript run_rolypoly_multi_expr_datasets.R --gwas_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_WHRadjBMIz.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_WHRadjBMIz --run_name squared_tss_10kb --n_cores 10 |& tee log.body_WHRadjBMIz.out.txt
### body_BMIz [MISSING]
# time Rscript run_rolypoly_multi_expr_datasets.R --gwas_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMIz.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_BMIz --run_name squared_tss_10kb --n_cores 10 |& tee log.body_BMIz.out.txt
### disease_T2D.gwassumstats
# time Rscript run_rolypoly_multi_expr_datasets.R --gwas_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/disease_T2D.gwassumstats.rolypoly_fmt.tab.gz --gwas_name disease_T2D --run_name squared_tss_10kb --n_cores 10 |& tee log.disease_T2D.out.txt
### body_HEIGHTz
# time Rscript run_rolypoly_multi_expr_datasets.R --gwas_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_HEIGHTz.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_HEIGHTz --run_name squared_tss_10kb --n_cores 10 |& tee log.body_HEIGHTz.out.txt
### mental_SCZ_Ripke2014
# time Rscript run_rolypoly_multi_expr_datasets.R --gwas_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/mental_SCZ_Ripke2014.gwassumstats.rolypoly_fmt.tab.gz --gwas_name mental_SCZ_Ripke2014 --run_name squared_tss_10kb --n_cores 10 |& tee log.mental_SCZ_Ripke2014.out.txt
### lipids_TC.Willer2013
# time Rscript run_rolypoly_multi_expr_datasets.R --gwas_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/lipids_TC.Willer2013.gwassumstats.rolypoly_fmt.tab.gz --gwas_name lipids_TC.Willer2013 --run_name squared_tss_10kb --n_cores 10 |& tee log.lipids_TC.Willer2013.out.txt

### body_BMI_Locke2015
# time Rscript run_rolypoly_multi_expr_datasets.R --gwas_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_BMI_Locke2015 --run_name squared_tss_10kb --n_cores 10 |& tee log.body_BMI_Locke2015.out.txt

### BMI test run
# time Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz --gwas_name TEST_body_BMI_Locke2015 --run_name gene.10kb.squared.protein_coding_only --window_position gene --window_size_kb 10 --pos_transformation square --protein_coding_only --n_cores 10 --test_run |& tee log.TEST_body_BMI_Locke2015.out.txt


# ======================================================================================================= #
# ============================================ OptParse ================================================= #
# ======================================================================================================= #

library(optparse)

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")

option_list <- list( 
  make_option("--gwas_file", type="character", default=NULL,
              help = "File path to GWAS file"),
  make_option("--gwas_name", type="character", default=NULL,
              help = "Prefix for output files"),
  make_option("--run_name", type="character", default="UNNAMED_RUN",
              help = "choose a name for the run, e.g. 'run_squared' if using the squared function for pos_transformation"),
  make_option("--outdir", type="character", default=NULL,
              help = "Full pathname for the output directory (will be created if it does not exists). Default is current working directory."),
  make_option("--gwas_linked_file", type="character", default=NULL,
              help = "[OPTIONAL] Full pathname to a pre-computed 'gwas linked' RolyPoly object. If argument is passed, then the gwas linked file will be loaded. OBS: This option overwrites loading a gwas_file with the same name ('caching' function)"),
  make_option("--window_position", type="character", default="tss",
              help = "Set the position of the window. Choose from 'tss' (transcription start site) and 'gene' (gene start and end)"),
  make_option("--window_size_kb", type="integer", default=10,
              help = "Window size in kilobases"),
  make_option("--pos_transformation", type="character", default="square",
              help = "Choose from 'square','abs','pos_only','none'"),
  make_option("--protein_coding_only", action="store_true", default=FALSE, 
              help="Run with protein coding genes only"),
  make_option("--test_run", action="store_true", default=FALSE, 
              help="Run in 'test mode': run on a subset of the GWAS data"),
  make_option("--n_cores", type="integer", default=10,
              help = "Number of cores to use for parallelization [default %default]")
  # make_option("--do_log_odds_beta", action="store_true", default=FALSE, 
  #             help="Log-transform GWAS beta column. Should only be used for case-control GWA studies"),
  
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

### get arguments
GWAS_FILE <- opt$gwas_file
GWAS_NAME <- opt$gwas_name
RUN_NAME <- opt$run_name
OUTDIR <- opt$outdir
GWAS_LINKED_FILE <- opt$gwas_linked_file
FLAG_TEST_RUN <- opt$test_run 
N_CORES <- opt$n_cores
N_BOOTSTRAP_ITERS <- 500 # number of bootstrap iterations

### get parameters
PARAM.WINDOW_POSITION <- opt$window_position
PARAM.WINDOW_SIZE_KB <- opt$window_size_kb
PARAM.POS_TRANSFORMATION <- opt$pos_transformation
PARAM.PROTEIN_CODING_ONLY <- opt$protein_coding_only

### Validate parameters 
if (is.null(GWAS_FILE)) { # if no argument given
  stop("Error: no gwas_file given")
}

if (is.null(GWAS_NAME)) { # if no argument given
  stop("Error: no gwas_name given")
}

if (is.null(OUTDIR)) {
  OUTDIR <- getwd() # get current working directory
  message(sprintf("No --outdir argument given. Will use current working directory: %s", OUTDIR))
} else {
  dir.create(OUTDIR,showWarnings=F, recursive=T) # make new directory. cmd does nothing if the directory already exists.
}

if (!is.null(GWAS_LINKED_FILE)) {
  stopifnot(file.exists(GWAS_LINKED_FILE)) # check that file exists
}


if (!PARAM.WINDOW_POSITION %in% c("tss", "gene")) { #
  stop(sprintf("Error: wrong argument for window_position given: %s", PARAM.WINDOW_POSITION))
}


print(sprintf("========== RUNNING GWAS = %s (%s) ==========", GWAS_NAME, GWAS_FILE))

# ======================================================================================================= #
# ==========================================       CONFIG      =========================================== #
# ======================================================================================================= #

### save arguments to list
list.run_parameters <- list(PARAM.WINDOW_POSITION = PARAM.WINDOW_POSITION,
                            PARAM.WINDOW_SIZE_KB = PARAM.WINDOW_SIZE_KB,
                            PARAM.POS_TRANSFORMATION = PARAM.POS_TRANSFORMATION,
                            PARAM.PROTEIN_CODING_ONLY = PARAM.PROTEIN_CODING_ONLY,
                            N_BOOTSTRAP_ITERS=N_BOOTSTRAP_ITERS,
                            GWAS_FILE=GWAS_FILE,
                            GWAS_NAME=GWAS_NAME,
                            RUN_NAME=RUN_NAME,
                            GWAS_LINKED_FILE=GWAS_LINKED_FILE,
                            FLAG_TEST_RUN=FLAG_TEST_RUN,
                            N_CORES=N_CORES)

print("============================ RUN PARAMETERS ===========================")
print(list.run_parameters)


# ======================================================================================================= #
# =============================================== DEBUGGING ============================================= #
# ======================================================================================================= #

# GWAS_FILE <- "/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz"
# GWAS_NAME <- "body_BMI_Locke2015"
# RUN_NAME <- "gene.10kb.squared.protein_coding_only"# or "gene.10kb.squared.all_genes"
# N_BOOTSTRAP_ITERS <- 500 # number of bootstrap iterations
# N_CORES <- 3 # number of cores for parallel
# FLAG_TEST_RUN = FALSE

# PARAM.WINDOW_POSITION = "gene"
# PARAM.WINDOW_SIZE_KB = 10
# #PARAM.POS_TRANSFORMATION = "square"
# PARAM.POS_TRANSFORMATION = "none"
# PARAM.PROTEIN_CODING_ONLY = FALSE

# ======================================================================================================= #
# ==========================================       SETUP      =========================================== #
# ======================================================================================================= #


suppressMessages(library(rolypoly))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))

dir.sc_genetics.data <- "/projects/timshel/sc-genetics/sc-genetics/data" # no trailing slash
dir.ldfiles <- file.path(dir.sc_genetics.data, "rolypoly/EUR_LD_FILTERED_NONAN_R") # Linkage disequilibrium files

dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library


# ======================================================================================================= #
# ============================================ OUTPUT FILES  ============================================ #
# ======================================================================================================= #

file.out.gwas_linked <- sprintf("%s/rolypoly_objs.%s.%s.gwas_linked.RData", OUTDIR, GWAS_NAME, RUN_NAME)
file.out.inference <- sprintf("%s/rolypoly_objs.%s.%s.inference.RData", OUTDIR, GWAS_NAME, RUN_NAME)

# ======================================================================================================= #
# ========================================= LOAD EXISTING DATA  ========================================= #
# ======================================================================================================= #

flag.loaded_gwas_linked_rdata <- FALSE # default value
flag.loaded_inference_rdata <- FALSE # default value

# *OBS* notice the order of this if-statement determines the priority for loading existing data: first cmd-line argument; then by name matching
if (!is.null(GWAS_LINKED_FILE)) { 
  print(sprintf("*Will load existing *gwas_linked* data file specified in the commend line argument: %s", GWAS_LINKED_FILE))
  load(GWAS_LINKED_FILE) # list.gwas_linked, list.run_parameters
    # list.gwas_linked ---> list("rp"=list(rp.gwas_linked), "runtime"=list(runtime.gwas_linked))
  rp.gwas_linked <- list.gwas_linked[["rp"]] # *OBS*: extract pre-computed object
  flag.loaded_gwas_linked_rdata <- TRUE
  print("Done...")
} else if (file.exists(file.out.gwas_linked)) {
  print(sprintf("*Will load existing *gwas_linked* data file WITH the same GWAS_NAME and RUN_NAME*: %s", file.out.gwas_linked))
  load(file.out.gwas_linked) # list.gwas_linked, list.run_parameters
    # list.gwas_linked ---> list("rp"=list(rp.gwas_linked), "runtime"=list(runtime.gwas_linked))
  rp.gwas_linked <- list.gwas_linked[["rp"]] # *OBS*: extract pre-computed object
  flag.loaded_gwas_linked_rdata <- TRUE
  print("Done...")
}


if (file.exists(file.out.inference)) {
  print(sprintf("*Will load existing inference data file WITH the same GWAS_NAME and RUN_NAME*: %s", file.out.inference))
  load(file.out.inference) # list.rp_inference, list.run_parameters
  flag.loaded_inference_rdata <- TRUE
  print("Done...")
}



# ======================================================================================================= #
# ============================================== GWAS DATA ============================================== #
# ======================================================================================================= #

# file.gwas <- file.path(dir.sc_genetics.data, "gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz") # BMI
# file.gwas <- file.path(dir.sc_genetics.data, "gwas_sumstats/mental_SCZ_Ripke2014.gwassumstats.rolypoly_fmt.tab.gz") # SCZ

df.gwas.rp <- read_gwas(GWAS_FILE, exlcude.HLA=T, do.log_odds=F)


if (FLAG_TEST_RUN) {
  print("RUNNING IN GWAS TEST MODE!")
  df.gwas.rp <- df.gwas.rp %>% slice(1:1000)
}

# ======================================================================================================= #
# ====================================== EXPRESSION DATA FILES  ========================================= #
# ======================================================================================================= #


### t-test data
# file.expr <- file.path(dir.sc_genetics.data, "expression/arc_lira/arc_lira.celltype_expr.ttest.hsapiens_orthologs.csv.gz") # arc_lira
# file.expr <- file.path(dir.sc_genetics.data, "expression/maca/maca.per_tissue.celltype_expr.ttest.hsapiens_orthologs.csv.gz") # maca per_tissue

### average expression
# file.expr <- file.path(dir.sc_genetics.data, "expression/maca/maca.per_tissue.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz") # maca per_tissue
# file.expr <- file.path(dir.sc_genetics.data, "expression/maca/maca.per_tissue_celltype.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz") # maca per_tissue_celltype


#files <- list.files(file.path(dir.sc_genetics.data, "expression"), recursive=T, pattern="*avg_expr.hsapiens_orthologs*") # matching on avg_expr
files <- list.files(file.path(dir.sc_genetics.data, "expression"), recursive=T, pattern="*.hsapiens_orthologs*") # matching on all data sets
files.paths <- file.path(dir.sc_genetics.data, "expression", files)
names(files.paths) <- files # named vector

print(files.paths)
# example:
# arc_lira/arc_lira.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz 
# "/projects/timshel/sc-genetics/sc-genetics/data/expression/arc_lira/arc_lira.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz" 
# gtex/gtex.sub_tissue.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz 
# "/projects/timshel/sc-genetics/sc-genetics/data/expression/gtex/gtex.sub_tissue.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz" 

# ======================================================================================================= #
# ========================================== GENE ANNOTATION =========================================== #
# ======================================================================================================= #

# file.gene_annot <- file.path(dir.sc_genetics.data, "gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz") # all 1-1 human-mouse ortholog
file.gene_annot <- file.path(dir.sc_genetics.data, "gene_annotations/gene_annotation.hsapiens_all_genes.GRCh37.ens_v91.txt.gz") # all human genes

df.gene_annot <- read_gene_annotation(file.gene_annot, protein_coding_only=PARAM.PROTEIN_CODING_ONLY, delim="\t")
# df.gene_annot %>% head()

# ======================================================================================================= #
# ====================================== LOAD EXPRESSION DATA  ========================================= #
# ======================================================================================================= #


wrapper.read_expression_data <- function(file.expr) {
  if (grepl("ttest.hsapiens_orthologs", file.expr)) {
    method_specific_expr="tstat"
  } else if (grepl("avg_expr.hsapiens_orthologs", file.expr)) {
    method_specific_expr="average"
  }
  df.expr <- read_expression_data(file.expr=file.expr,
                                  method_specific_expr=method_specific_expr, 
                                  pos_transformation=PARAM.POS_TRANSFORMATION, 
                                  genes.filter=df.gene_annot$ensembl_gene_id,  
                                  col_genes="gene", 
                                  delim=",") 
  return(df.expr)
}

### read data
list.df_expr <- lapply(files.paths, wrapper.read_expression_data) # returns named list
# names(list.df_expr)
# [1] "arc_lira/arc_lira.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz"            
# [2] "gtex/gtex.sub_tissue.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz"         
# [3] "gtex/gtex.tissue.gene_tpm.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz"    
# [4] "maca/maca.per_celltype.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz"       
# [5] "maca/maca.per_tissue_celltype.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz"
# [6] "maca/maca.per_tissue.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz"   

### List the union of all the genes
# all_genes <- unique(unlist(sapply(list.df_expr, rownames)))

# ======================================================================================================= #
# ========================================== BLOCK ANNOTATION =========================================== #
# ======================================================================================================= #




# COLUMNS in df.gene_annot --> ensembl_gene_id chromosome_name start_position end_position   gene_biotype
if (PARAM.WINDOW_POSITION == "gene") {
  window_boundary_end <- "end_position"
} else if (PARAM.WINDOW_POSITION == "tss") {
  window_boundary_end <- "start_position" # assumuing
} else {
  
}

# REQUIRED COLUMNS IN block_annotation DATA FRAME: c('chrom', 'start', 'end', 'label')
df.block_annotation <- df.gene_annot %>% transmute(
  label=ensembl_gene_id,
  chrom=chromosome_name,
  start=start_position-(PARAM.WINDOW_SIZE_KB * 1000), # TSS-PARAM.WINDOW_SIZE_KB
  end=UQ(rlang::sym(window_boundary_end))+(PARAM.WINDOW_SIZE_KB * 1000) # TSS+PARAM.WINDOW_SIZE_KB
)



# ======================================================================================================= #
# ======================================== ROLY FUNCTION ================================================ #
# ======================================================================================================= #

run_rolypoly.multi_attempt <- function(rp.precomputed,
                                       gwas_data,
                                       block_annotation,
                                       block_data,
                                       ld_folder,
                                       n_attempts) {
  N_MAX_ATTEMPTS <- 10 # constant. max number of attempts
  n_attempts <- n_attempts + 1 # incrementing attempts
  
  rp <- tryCatch(
    {
      # 'tryCatch()' will return the last evaluated expression  in case the "try" part was completed successfull
      message(sprintf("Attempt = %s | *Try block* | running attempt...", n_attempts))
      rolypoly_roll(rolypoly = rp.precomputed, # RolyPoly object precomputed (linked GWAS and gene annotation data)
                    gwas_data = gwas_data,
                    block_annotation = block_annotation,
                    block_data = block_data,
                    ld_folder = ld_folder,
                    add_poly = T, # it is ok for this to be set to TRUE, even though we do not link the gwas and gene annotations. (add_poly is only used when loading the gwas data, i.e. when gwas_data is not NULL)
                    n_degree = 1, # same as for add_poly
                    bootstrap_iters = N_BOOTSTRAP_ITERS, # use at least 200
                    gwas_link_parallel = T, # it is ok for this to be set to TRUE, even though we do not link the gwas and gene annotations.
                    bootstrap_parallel = T)
      # The return value of the above expression is the actual valuethat will be returned in case there is no condition (e.g. warning or error). 
      # You don't need to state the return value via `return()`.
    },
    error=function(cond) {
      message(sprintf("Attempt = %s | *Error block* | An error happend! Here's the original error message:", n_attempts))
      message(cond)
      message() # gives extra space.
      if (n_attempts < N_MAX_ATTEMPTS) {
        message(sprintf("Attempt = %s | *Error block* | Will run function again...", n_attempts))
        ### Calling the function with the SAME arguments as it was called the first time.
        run_rolypoly.multi_attempt(rp.precomputed = rp.precomputed,
                                  gwas_data = gwas_data, 
                                  block_annotation = block_annotation,
                                  block_data = block_data,
                                  ld_folder = ld_folder,
                                  n_attempts = n_attempts) # *OBS*: calling with updated variable
      } else {
        message(sprintf("Attempt = %s | *Error block* | REACHED MAX ATTEMPTS - will return 'NULL'", n_attempts))
        return(NULL) # return value if everything keeps failing
      }
    } # end error handling
  ) # end tryCatch
  return(rp) # return rolypoly
}


# ======================================================================================================= #
# ======================================== START PARALLEL ================================================ #
# ======================================================================================================= #

print("Starting registerDoParallel...")
library(doParallel)
registerDoParallel(N_CORES)
getDoParWorkers()


# ======================================================================================================= #
# ======================================== RUN FIRST RUN ================================================ #
# ======================================================================================================= #

if (!flag.loaded_gwas_linked_rdata) {
  print("Making first RolyPoly call to link GWAS data...")
  
  ### Run RolyPoly
  runtime.gwas_linked <- system.time(
    rp.gwas_linked <- run_rolypoly.multi_attempt(rp.precomputed = NULL, # no object parsed, so it will be created.
                                               gwas_data = as.data.frame(df.gwas.rp), 
                                               block_annotation = as.data.frame(df.block_annotation),
                                               block_data = as.data.frame(list.df_expr[[1]]), # first expression data set
                                               ld_folder = dir.ldfiles,
                                               n_attempts = 0) # *OBS*: important to make first call with n_attempts=0
  )
  print("Done with RolyPoly first run. Runtime:")
  print(runtime.gwas_linked)
  print("Printing any warnings from first run:")
  print(warnings())
  
  ### Gene scores
  print("Calculating gene scores...")
  rp.gwas_linked <- rolypoly_add_ld_corrected_gwas_block_scores(rolypoly=rp.gwas_linked, fast_calculation=T) # returns RolyPoly object | this takes time - it is calculated per gene
    # sets the slot rp$data$<BLOCK/GENE>$corrected_gwas_block_score_pval
  ### You might get the warning:
  # In CompQuadForm::imhof(raw_score, lambda = lambda, epsabs = 1e-15,  ... :
  #                          Note that Qq + abserr is positive.
  print("Done calculating gene scores")
  
  ### Combine into list object
  list.gwas_linked <- list("rp"=list(rp.gwas_linked), "runtime"=list(runtime.gwas_linked))
  
  ### Save
  print("Saving...")
  save(list.gwas_linked, list.run_parameters, file=file.out.gwas_linked)
}


### Error
# starting bootstrap iteration: 492
# Attempt = 1 | *Error block* | An error happend! Here's the original error message:
# task 4 failed - "variable lengths differ (found for 'vectorized_rolypoly_data$x')"
# Attempt = 1 | *Error block* | Will run function again...
# Attempt = 2 | *Try block* | running attempt...
# adding gwas
# filtering out SNPs with MAF < 1%
# adding block annotations
# beginning processs of linking gwas a

# Attempt = 1 | *Error block* | An error happend! Here's the original error message:
# task 394 failed - "variable lengths differ (found for 'vectorized_rolypoly_data$x')"
# Attempt = 1 | *Error block* | Will run function again...
# Attempt = 2 | *Try block* | running attempt...
# adding gwas
# filtering out SNPs with MAF < 1%
# adding block annotations


# ======================================================================================================= #
# ======================================= WRAPPER FUNCTION ============================================== #
# ======================================================================================================= #

wrapper.run_rolypoly_over_expr_datasets <- function(idx.list, list, rp.precomputed) {
  ### INPUT
  # idx.list     a integer specifying the index of the list to run the function on.
  # list         a list that contains the index provided by "idx.list". (Preferably the list is named)
  
  n.elements <- length(list)
  if (is.null(names(list))) { # list is unnamed
    name.list <- "<UNNAMED LIST ELEMENT>"
  } else {
    name.list <- names(list)[idx.list] # get list element name
  }
  print(sprintf("========== #%s/#%s | Status = START RUNNING | wrapper.run_rolypoly_over_expr_datasets | list element = %s ==========", idx.list, n.elements, name.list)) # prints "character(0)" if list has no names (i.e. names(list) is NULL)
  df.expr <- list[[idx.list]] # get expression data frame
  runtime <- system.time(
    rp <- run_rolypoly.multi_attempt(rp.precomputed = rp.precomputed,
                                     gwas_data = NULL, # set to NULL to avoid GWAS-linking computation
                                     block_annotation = NULL, # set to NULL to avoid GWAS-linking computation
                                     block_data = as.data.frame(df.expr), # *VARIABLE*
                                     ld_folder = NULL, # set to NULL to avoid GWAS-linking computation
                                     n_attempts = 0) # *OBS*: important to make first call with n_attempts=0
  )
  
  ### Deleting data slots from RolyPoly object
  # We do this to reduce later load times and save space (storage and memory).
  rp$data <- "REMOVED" # gene level information (GWAS-linked) [takes up 98% of the size of the object]
  # rp$blocks <- "REMOVED" # gene==block annotation
  # rp$raw_block_data <- "REMOVED" # gene==block expression data
  
  # system.time(): returns a names numerical vector of class "proc_time"
  res <- c("rp"=list(rp), "runtime"=list(runtime), "name"=name.list) # elements must be enclosed in list() to keep their native data structure/'integrity' (unless they are 'atomic')
  print(sprintf("========== #%s/#%s | Status = DONE | wrapper.run_rolypoly_over_expr_datasets | list element = %s ==========", idx.list, n.elements, name.list)) # prints "character(0)" if list has no names (i.e. names(list) is NULL)
  return(res)
}


# ======================================================================================================= #
# =================================== FILTER EXPRESSION DATA SETS TO RUN ================================ #
# ======================================================================================================= #

list.df_expr.run <- list.df_expr # default value
if (flag.loaded_inference_rdata) { # if we have loaded inference data
  # ASSUMPTION: list.rp_inference and list.df_expr share the same name for inference on expression datasets that have analyzed
  # > names(list.df_expr)
  # [1] "arc_lira/arc_lira.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz"            
  # [2] "gtex/gtex.sub_tissue.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz"  
  
  names.expr_datasets_already_run <- names(list.df_expr) %in% list.rp_inference # return TRUE for dataset names already analyzed
  list.df_expr.run <- list.df_expr[!names.expr_datasets_already_run]
}

# ======================================================================================================= #
# ========================================== RUN WRAPPER ================================================ #
# ======================================================================================================= #

print("Running wrapper.run_rolypoly_over_expr_datasets with the following (new) expression data sets:")
print(names(list.df_expr.run))

### Run wrapper
list.rp_inference.run <- lapply(seq_along(list.df_expr.run), wrapper.run_rolypoly_over_expr_datasets, list=list.df_expr.run, rp.precomputed=rp.gwas_linked) # "list" is here the extra (named) argument to wrapper.run_rolypoly_over_expr_datasets()
names(list.rp_inference.run) <- names(list.df_expr.run) # set names of list
print("Done with RolyPoly wrapper loop")
print("Printing any warnings:")
print(warnings())

# ======================================================================================================= #
# ========================================== ADDING NEW DATA TO OLD ===================================== #
# ======================================================================================================= #

if (flag.loaded_inference_rdata) { # if we have loaded inference data
  print("Adding new run data to previous inference data...")
  list.rp_inference <- modifyList(list.rp_inference, list.rp_inference.run) # combining lists using modifyList(). 
  # c() operator should also do the trick, but this is nicer if there - for some reason - are elements with the same name in both lists.
  # modifyList(x,y): elements from y will be added to the list; any elements in x which is also in y will be replaced by the y elements..
} else {
  list.rp_inference <- list.rp_inference.run
}

# ======================================================================================================= #
# ============================================ EXPORT ================================================== #
# ======================================================================================================= #

print("Saving inference object...")
save(list.rp_inference, list.run_parameters, file=file.out.inference)

# ======================================================================================================= #
# ============================================== FINISH ================================================== #
# ======================================================================================================= #

print("SCRIPT DONE!")

# ======================================================================================================= #
# ============================================== PLOTS ================================================== #
# ======================================================================================================= #

### Plot the log10(p) to rank tissues by the strength of their association
# p.pval_ranking <- plot_rolypoly_annotation_ranking(rp)
# p.pval_ranking


### Plot gamma estimate with 95% confidence intervals
### gamma: influence of cell-type-specific gene expression on the variance of GWAS effect sizes
# p.gamma_est <- plot_rolypoly_annotation_estimates(rp)
# p.gamma_est


# ======================================================================================================= #
# ============================================== MISC ================================================== #
# ======================================================================================================= #

# rp$full_results$parameters %>% sort
# rp$bootstrap_results %>% arrange(-bt_value) %>% head



# ======================================================================================================= #
# ==========================================   LEFT OVERS   ============================================ #
# ======================================================================================================= #
