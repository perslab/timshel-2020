############### SYNOPSIS ###################
# Running RolyPoly

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
# time Rscript run_rolypoly_multi_expr_datasets.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_WHR_Shungin2015.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_WHR_Shungin2015 --run_name squared_tss_10kb --n_cores 10 |& tee log.body_WHR_Shungin2015.out.txt
### cov_EDU_YEARS
# time Rscript run_rolypoly_multi_expr_datasets.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/cov_EDU_YEARS.gwassumstats.rolypoly_fmt.tab.gz --gwas_name cov_EDU_YEARS --run_name squared_tss_10kb --n_cores 10 |& tee log.cov_EDU_YEARS.out.txt
### body_WHRadjBMIz
# time Rscript run_rolypoly_multi_expr_datasets.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_WHRadjBMIz.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_WHRadjBMIz --run_name squared_tss_10kb --n_cores 10 |& tee log.body_WHRadjBMIz.out.txt
### body_BMIz [MISSING]
# time Rscript run_rolypoly_multi_expr_datasets.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMIz.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_BMIz --run_name squared_tss_10kb --n_cores 10 |& tee log.body_BMIz.out.txt
### disease_T2D.gwassumstats
# time Rscript run_rolypoly_multi_expr_datasets.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/disease_T2D.gwassumstats.rolypoly_fmt.tab.gz --gwas_name disease_T2D --run_name squared_tss_10kb --n_cores 10 |& tee log.disease_T2D.out.txt
### body_HEIGHTz
# time Rscript run_rolypoly_multi_expr_datasets.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_HEIGHTz.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_HEIGHTz --run_name squared_tss_10kb --n_cores 10 |& tee log.body_HEIGHTz.out.txt
### mental_SCZ_Ripke2014
# time Rscript run_rolypoly_multi_expr_datasets.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/mental_SCZ_Ripke2014.gwassumstats.rolypoly_fmt.tab.gz --gwas_name mental_SCZ_Ripke2014 --run_name squared_tss_10kb --n_cores 10 |& tee log.mental_SCZ_Ripke2014.out.txt
### lipids_TC.Willer2013
# time Rscript run_rolypoly_multi_expr_datasets.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/lipids_TC.Willer2013.gwassumstats.rolypoly_fmt.tab.gz --gwas_name lipids_TC.Willer2013 --run_name squared_tss_10kb --n_cores 10 |& tee log.lipids_TC.Willer2013.out.txt

### body_BMI_Locke2015
# time Rscript run_rolypoly_multi_expr_datasets.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_BMI_Locke2015 --run_name squared_tss_10kb --n_cores 10 |& tee log.body_BMI_Locke2015.out.txt

### BMI test run
# time Rscript run_rolypoly_multi_expr_datasets.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz --gwas_name TEST_body_BMI_Locke2015 --run_name squared_tss_10kb --n_cores 10 --test_run |& tee log.TEST_body_BMI_Locke2015.out.txt


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

if (is.null(opt$gwas_file)) { # if no argument given
  stop("Error: no gwas_file given")
}

if (is.null(opt$gwas_name)) { # if no argument given
  stop("Error: no gwas_name given")
}

### get arguments
GWAS_FILE <- opt$gwas_file
GWAS_NAME <- opt$gwas_name
RUN_NAME <- opt$run_name
FLAG_TEST_RUN <- opt$test_run 
N_CORES <- opt$n_cores


print(sprintf("========== RUNNING GWAS = %s (%s) ==========", GWAS_NAME, GWAS_FILE))


# LOAD_EXISTING_FIRST_RUN <- TRUE # *OBS*
LOAD_EXISTING_FIRST_RUN <- FALSE # *OBS*

# ======================================================================================================= #
# ==========================================       CONFIG      =========================================== #
# ======================================================================================================= #

POS_TRANSFORMATION <- "square" # for expression data
WINDOW_SIZE <- 10e3 # 10kb
N_BOOTSTRAP_ITERS <- 500 # number of bootstrap iterations

### save arguments to list
list.run_parameters <- list(N_BOOTSTRAP_ITERS=N_BOOTSTRAP_ITERS,
                            POS_TRANSFORMATION=POS_TRANSFORMATION,
                            WINDOW_SIZE=WINDOW_SIZE,
                            GWAS_FILE=GWAS_FILE,
                            GWAS_NAME=GWAS_NAME,
                            RUN_NAME=RUN_NAME,
                            FLAG_TEST_RUN=FLAG_TEST_RUN,
                            N_CORES=N_CORES)

print(list.run_parameters)

# GWAS_FILE <- "DUMMY"
# GWAS_NAME <- "SCZ"
# RUN_NAME <- "run_squared"
# N_CORES <- 3 # number of cores for parallel
# FLAG_TEST_RUN = FALSE


# ======================================================================================================= #
# ==========================================       SETUP      =========================================== #
# ======================================================================================================= #

library(rolypoly)
library(tidyverse)
library(ggplot2)

dir.sc_genetics.data <- "/projects/timshel/sc-genetics/sc-genetics/data" # no trailing slash
dir.ldfiles <- file.path(dir.sc_genetics.data, "rolypoly/EUR_LD_FILTERED_NONAN_R") # Linkage disequilibrium files

dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library


# ======================================================================================================= #
# ============================================ OUTPUT FILES  ============================================ #
# ======================================================================================================= #


file.out.final <- sprintf("out.rolypoly_objs.%s.%s.final.RData", GWAS_NAME, RUN_NAME)
file.out.first_run <- sprintf("out.rolypoly_objs.%s.%s.first_run.RData", GWAS_NAME, RUN_NAME)

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


files <- list.files(file.path(dir.sc_genetics.data, "expression"), recursive=T, pattern="*avg_expr.hsapiens_orthologs*")
files.paths <- file.path(dir.sc_genetics.data, "expression", files)
names(files.paths) <- files

print(files.paths)


# ======================================================================================================= #
# ========================================== GENE ANNOTATION =========================================== #
# ======================================================================================================= #

# file.gene_annot <- file.path(dir.sc_genetics.data, "gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz") # all 1-1 human-mouse ortholog
file.gene_annot <- file.path(dir.sc_genetics.data, "gene_annotations/gene_annotation.hsapiens_all_genes.GRCh37.ens_v91.txt.gz") # all human genes

df.gene_annot <- read_gene_annotation(file.gene_annot, delim="\t")
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
                                  pos_transformation=POS_TRANSFORMATION, 
                                  genes.filter=df.gene_annot$ensembl_gene_id,  
                                  col_genes="gene", 
                                  delim=",") 
  return(df.expr)
}

### read data
list.df_expr <- lapply(files.paths, wrapper.read_expression_data)


### List the union of all the genes
# all_genes <- unique(unlist(sapply(list.df_expr, rownames)))

# ======================================================================================================= #
# ========================================== BLOCK ANNOTATION =========================================== #
# ======================================================================================================= #

# REQUIRED COLUMNS IN DATA FRAME: c('chrom', 'start', 'end', 'label')

df.block_annotation <- df.gene_annot %>% transmute(
  label=ensembl_gene_id,
  chrom=chromosome_name,
  start=start_position-WINDOW_SIZE, # TSS-WINDOW_SIZE
  end=start_position+WINDOW_SIZE # TSS+WINDOW_SIZE
)

# df.block_annotation %>% head


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

if (LOAD_EXISTING_FIRST_RUN) {
  print(sprintf("*Will load existing first_run data file*: %s", file.out.first_run))
  load(file=file.out.first_run)
} else {
  print("Making first RolyPoly call...")
  runtime.first_run <- system.time(
    rp.first_run <- run_rolypoly.multi_attempt(rp.precomputed = NULL, # no object parsed, so it will be created.
                                               gwas_data = as.data.frame(df.gwas.rp), 
                                               block_annotation = as.data.frame(df.block_annotation),
                                               block_data = as.data.frame(list.df_expr[[1]]), # first expression data set
                                               ld_folder = dir.ldfiles,
                                               n_attempts = 0) # *OBS*: important to make first call with n_attempts=0
  )
  print("Done with RolyPoly first run. Runtime:")
  print(runtime.first_run)
  print("Printing any warnings from first run:")
  print(warnings())
  
  ### Save
  print("Saving...")
  save(rp.first_run, runtime.first_run, 
       list.run_parameters, file=file.out.first_run)
}

# Warning messages:
#   1: In calculate_block_values(rolypoly$raw_block_data, rolypoly$full_results$parameters) :
#   NA pameters found, take results with grain of salt
# 2: In calculate_annotation_block_heritability(rolypoly$raw_block_data,  :
#   NA pameters found, take results with grain of salt


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
# ========================================== RUN WRAPPER ================================================ #
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
  # system.time(): returns a names numerical vector of class "proc_time"
  list.res <- c("rp"=list(rp), "runtime"=list(runtime), "name"=name.list) # elements must be enclosed in list() to keep their native data structure/'integrity'
  # print(sprintf("========== #%s/#%s | Status = DONE | wrapper.run_rolypoly_over_expr_datasets | list element = %s ==========", idx.list, n.elements, name.list)) # prints "character(0)" if list has no names (i.e. names(list) is NULL)
  # print(sprintf("Will save temporary snapshot of variables to RData file: %s", file.out.final))
  # save(rp.first_run, runtime.first_run, 
  #      list.run_parameters,
  #      list.rp_wrapper, file=file.out.final)
  # print("Done saving...")
  return(list.res)
}

print("Running wrapper.run_rolypoly_over_expr_datasets with the following files:")
print(files.paths)

### Run wrapper
list.rp_wrapper <- lapply(seq_along(list.df_expr), wrapper.run_rolypoly_over_expr_datasets, list=list.df_expr, rp.precomputed=rp.first_run) # "list" is here the extra (named) argument to wrapper.run_rolypoly_over_expr_datasets()
names(list.rp_wrapper) <- names(list.df_expr) # set name of list
print("Done with RolyPoly wrapper loop")
print("Printing any warnings:")
print(warnings())

# ======================================================================================================= #
# ============================================ EXPORT ================================================== #
# ======================================================================================================= #


save(rp.first_run, runtime.first_run, 
     list.run_parameters,
     list.rp_wrapper, file=file.out.final)

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
