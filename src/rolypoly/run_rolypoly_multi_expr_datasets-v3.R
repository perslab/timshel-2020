############### SYNOPSIS ###################
# Running RolyPoly

### v3 change log: 
# * univariate regression: change FROM *multiple regression* over all cell-types TO *multiple "univariate"* regressions. I.e. each cell-type is run separately to avoid colinarity issues.
# * Restructure code: split into 1) linked_gwas [rolypoly_load_gwas()] snp_annotation arguments; 2) rolypoly_perform_inference. This was done to allow support of loading snp_annotations.
# * added support for WGCNA 'kme' 'expression data'
# * fixed issue --gwas_linked file loading
# * added --expr_data_list to read specific expression datasets
# + LIKELY MORE CHANGES

### v3 change log (2019-01-30 update to accomodate new data)
# * df.gwas.rp <- read_gwas(file.gwas=GWAS_FILE, .. ) --> update colnames

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
# time Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz --gwas_name TEST_body_BMI_Locke2015 --run_name gene.10kb.squared.protein_coding_only --window_position gene --window_size_kb 10 --pos_transformation square --protein_coding_only --n_cores 10 --test_run |& tee log.TEST_body_BMI_Locke2015.out.txt


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
              help = "[OPTIONAL] Choose a name for the run. Used to construct output filename. E.g. 'run_squared' if using the squared function for pos_transformation"),
  make_option("--outdir", type="character", default=NULL,
              help = "Full pathname for the output directory (will be created if it does not exists). Default is current working directory."),
  make_option("--gwas_linked_file", type="character", default=NULL,
              help = "[OPTIONAL] Full pathname to a pre-computed 'gwas linked' RData file. If argument is passed, then the gwas linked file will be loaded. OBS: This option overwrites the option for loading a gwas_file ('caching' function)"),
  # make_option("--inference_file", type="character", default=NULL,
  #            help = "[OPTIONAL] Full pathname to a 'inference' RolyPoly RData file. If argument is passed, then the inference file will be loaded."),
  make_option("--expr_data_list", type="character", default=NULL,
              help = "[OPTIONAL] Full pathname to a tab delimited file specifying what expression datasets to analyse. File format: no header; use '#' for outcommenting lines; Col1=expr. dataset name; Col2=full path to expr. dataset"),
  make_option("--window_position", type="character", default="tss",
              help = "Set the position of the window. Choose from 'tss' (transcription start site) and 'gene' (gene start and end)"),
  make_option("--window_size_kb", type="integer", default=10,
              help = "Window size in kilobases"),
  make_option("--pos_transformation", type="character", default="square",
              help = "Choose from 'square','abs','pos_only','none'"),
  make_option("--protein_coding_only", action="store_true", default=FALSE, 
              help="Run with protein coding genes only"),
  make_option("--n_bootstrap", type="integer", default=100, 
              help="Number of bootstrap iterations to run. Diego used approx 500-1000. Use at minimum 200 unless testing"),
  make_option("--gwas_linking_only", action="store_true", default=FALSE, 
              help="Only run GWAS linking step"),
  make_option("--ldsc_input_mode", action="store_true", default=FALSE, 
              help="Run in 'ldsc input mode': set the 'se' column in the GWAS data to a constant value. (Typically the input 'se' column will be all NA's"),
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
EXPR_DATA_LIST <- opt$expr_data_list
GWAS_LINKING_ONLY <- opt$gwas_linking_only
LDSC_INPUT_MODE <- opt$ldsc_input_mode # *OBS* tmp
FLAG_TEST_RUN <- opt$test_run 
N_CORES <- opt$n_cores
N_BOOTSTRAP_ITERS <- opt$n_bootstrap

### get parameters
PARAM.WINDOW_POSITION <- opt$window_position
PARAM.WINDOW_SIZE_KB <- opt$window_size_kb
PARAM.POS_TRANSFORMATION <- opt$pos_transformation
PARAM.PROTEIN_CODING_ONLY <- opt$protein_coding_only

### Validate parameters 
if (is.null(GWAS_FILE) & is.null(GWAS_LINKED_FILE)) { # if no arguments given
  stop("Error: no --gwas_file or --gwas_linked_file argument provided. Please provide one of them.")
} else if (!is.null(GWAS_FILE) & !is.null(GWAS_LINKED_FILE)) { # both arguments given
    stop("Error: arguments provided to both --gwas_file and --gwas_linked_file. The script requires that only one of these arguments are provided.")
} else if (!is.null(GWAS_FILE)) { # if only GWAS_FILE file given
  stopifnot(file.exists(GWAS_FILE)) # check that file exists
} else if (!is.null(GWAS_LINKED_FILE)) { # if only GWAS_LINKED_FILE given
  stopifnot(file.exists(GWAS_LINKED_FILE)) # check that file exists
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

if (!is.null(GWAS_LINKED_FILE) & FLAG_TEST_RUN) {
  stop("Error: test mode does not work with --gwas_linked_file")
}

if (!is.null(EXPR_DATA_LIST)) {
  stopifnot(file.exists(EXPR_DATA_LIST)) # check that file exists
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
                            EXPR_DATA_LIST=EXPR_DATA_LIST,
                            FLAG_TEST_RUN=FLAG_TEST_RUN,
                            N_CORES=N_CORES)

print("============================ RUN PARAMETERS ===========================")
print(list.run_parameters)


# ======================================================================================================= #
# =============================================== DEBUGGING ============================================= #
# ======================================================================================================= #

# GWAS_FILE <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz"
# GWAS_NAME <- "body_BMI_Locke2015"
# RUN_NAME <- "TMPTMPTMP_debug.gene.10kb.squared.protein_coding_only"# or "gene.10kb.squared.all_genes"
# OUTDIR <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2-tmp_test"
# GWAS_LINKED_FILE <- "/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.rolypoly_objs-v2.pos_only/rolypoly_objs.body_BMI_Locke2015.gene.10kb.squared.protein_coding_only.gwas_linked.RData"
# EXPR_DATA_LIST <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/list_expr_data.test.txt"
# N_BOOTSTRAP_ITERS <- 500 # number of bootstrap iterations
# N_CORES <- 3 # number of cores for parallel
# FLAG_TEST_RUN = FALSE
# 
# PARAM.WINDOW_POSITION = "gene"
# PARAM.WINDOW_SIZE_KB = 10
# PARAM.POS_TRANSFORMATION = "square"
# ###PARAM.POS_TRANSFORMATION = "none"
# PARAM.PROTEIN_CODING_ONLY = TRUE

# ======================================================================================================= #
# ==========================================       SETUP      =========================================== #
# ======================================================================================================= #


suppressMessages(library(rolypoly))
suppressMessages(library(tidyverse))
suppressMessages(library(doParallel))

dir.sc_genetics.data <- "/projects/timshel/sc-genetics/sc-genetics/data" # no trailing slash
dir.ldfiles <- file.path(dir.sc_genetics.data, "rolypoly/EUR_LD_FILTERED_NONAN_R") # Linkage disequilibrium files

dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library

### Check that correct version of RolyPoly is loaded
if (!exists('is_rolypoly_pascaltimshel_fork', mode='function')) {
  # 'is_rolypoly_pascaltimshel_fork' function only exists in pascaltimshel fork.
  # exists() returns TRUE only if function is found.
  stop("You have not loaded the pascaltimshel forked version of rolypoly. Univariate functionality is only supported in pascaltimshel forked version. Will exit...")
}


# ======================================================================================================= #
# ============================================ OUTPUT FILES  ============================================ #
# ======================================================================================================= #

file.out.gwas_linked <- sprintf("%s/rolypoly_objs.%s.%s.gwas_linked.RData", OUTDIR, GWAS_NAME, RUN_NAME)
file.out.inference <- sprintf("%s/rolypoly_objs.%s.%s.inference.RData", OUTDIR, GWAS_NAME, RUN_NAME)
file.out.inference.tmp <- sprintf("%s/rolypoly_objs.%s.%s.inference.tmp.RData", OUTDIR, GWAS_NAME, RUN_NAME)
file.out.inference.tmp_red <- sprintf("%s/rolypoly_objs.%s.%s.inference.tmp_red.RData", OUTDIR, GWAS_NAME, RUN_NAME)


# ======================================================================================================= #
# ========================================= LOAD EXISTING DATA  ========================================= #
# ======================================================================================================= #

flag.loaded_gwas_linked_rdata <- FALSE # default value
flag.loaded_inference_rdata <- FALSE # default value

# *OBS* notice the order of this if-statement determines the priority for loading existing data: first cmd-line argument; then by name matching
if (!is.null(GWAS_LINKED_FILE)) { 
  print(sprintf("*Will load existing *gwas_linked* data file specified in the commend line argument: %s", GWAS_LINKED_FILE))
  load(GWAS_LINKED_FILE) # list.gwas_linked, list.run_parameters
    # list.gwas_linked ---> list("rp"=rp.gwas_linked, "runtime"=runtime.gwas_linked)
  print(sprintf("Size list.gwas_linked: %s", format(object.size(list.gwas_linked),units="MB")))
  rp.gwas_linked <- list.gwas_linked[["rp"]] # *OBS*: extract pre-computed object
  flag.loaded_gwas_linked_rdata <- TRUE
  print("Done...")
} else if (file.exists(file.out.gwas_linked)) {
  print(sprintf("*Will load existing *gwas_linked* data file WITH the same GWAS_NAME and RUN_NAME*: %s", file.out.gwas_linked))
  load(file.out.gwas_linked) # list.gwas_linked, list.run_parameters
    # list.gwas_linked ---> list("rp"=list(rp.gwas_linked), "runtime"=list(runtime.gwas_linked))
  print(sprintf("Size list.gwas_linked: %s", format(object.size(list.gwas_linked),units="MB")))
  rp.gwas_linked <- list.gwas_linked[["rp"]] # *OBS*: extract pre-computed object
  flag.loaded_gwas_linked_rdata <- TRUE
  print("Done...")
}


if (file.exists(file.out.inference)) {
  print(sprintf("*Will load existing inference data file WITH the same GWAS_NAME and RUN_NAME*: %s", file.out.inference))
  load(file.out.inference) # list.rp_inference, list.run_parameters
  print(sprintf("Size list.rp_inference: %s", format(object.size(list.rp_inference),units="MB")))
  flag.loaded_inference_rdata <- TRUE
  print("Done...")
}



# ======================================================================================================= #
# ============================================== GWAS DATA ============================================== #
# ======================================================================================================= #

# file.gwas <- file.path(dir.sc_genetics.data, "gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz") # BMI
# file.gwas <- file.path(dir.sc_genetics.data, "gwas_sumstats/mental_SCZ_Ripke2014.gwassumstats.rolypoly_fmt.tab.gz") # SCZ

if (!flag.loaded_gwas_linked_rdata) { # only load GWAS data if not already loaded precomputed object.

  # df.gwas.rp <- read_gwas(GWAS_FILE, exlcude.HLA=T, do.log_odds=F) # columns c('chrom', 'pos', 'rsid', 'beta', 'se', 'maf')
    # ^ DEFAULT rolypoly format GWAS. 
  
  df.gwas.rp <- read_gwas(file.gwas=GWAS_FILE, 
                        exlcude.HLA=T, 
                        do.log_odds=F,
                        delim="\t", # file delimiter
                        col_chrom="chr",
                        col_pos="pos",
                        col_rsid="SNP", # NEW
                        col_beta="BETA", # NEW
                        col_se="SE", # NEW
                        col_maf="snp_maf" 
                        ) # 2019-01-30 format from LDSC munged GWAS (data/gwas_sumstats_ldsc/timshel-collection/)

  if (LDSC_INPUT_MODE) {
    print("RUNNING IN LDSC_INPUT_MODE. Will se values for 'se' column")
    df.gwas.rp$se <- 0.001
  }
  
  if (FLAG_TEST_RUN) {
    print("RUNNING IN GWAS TEST MODE!")
    df.gwas.rp <- df.gwas.rp %>% slice(1:1000)
  }

  
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


if (is.null(EXPR_DATA_LIST)) {
  stop("2019-01-30 update: --expr_data_list argument must be passed.")
  
  #files <- list.files(file.path(dir.sc_genetics.data, "expression"), recursive=T, pattern="*avg_expr.hsapiens_orthologs*") # matching on avg_expr
  files <- list.files(file.path(dir.sc_genetics.data, "expression"), recursive=T, pattern="*.hsapiens_orthologs*") # matching on all data sets
  files.paths <- file.path(dir.sc_genetics.data, "expression", files)
  names(files.paths) <- files # named vector
  
  # example:
  # arc_lira/arc_lira.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz 
  # "/projects/timshel/sc-genetics/sc-genetics/data/expression/arc_lira/arc_lira.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz" 
  # gtex/gtex.sub_tissue.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz 
  # "/projects/timshel/sc-genetics/sc-genetics/data/expression/gtex/gtex.sub_tissue.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz" 
} else {
  df.expr_data_list <- read_tsv(EXPR_DATA_LIST, col_names=FALSE, comment="#") # no header. Use "#" for commenting out lines.
  # print(sprintf("Received n=%s expression data from file %s", nrow(df.expr_data_list), EXPR_DATA_LIST))
  #print(df.expr_data_list)
  files.paths <- df.expr_data_list %>% pull(2) # extract second column as vector
  names(files.paths) <- df.expr_data_list %>% pull(1) # set names 
}

print(sprintf("Will read the following n=%s expression data files:", length(files.paths)))
print(files.paths)

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
    scale_genes=FALSE
  } else if (grepl("avg_expr.hsapiens_orthologs", file.expr)) {
    scale_genes=TRUE
  } else if (grepl("kme.hsapiens_orthologs", file.expr)) {
    scale_genes=FALSE
  } else if (grepl("depict_152tissues", file.expr)) { # *OBS*: semi-hack
    scale_genes=FALSE
  } else {
    scale_genes=FALSE
    # stop(sprintf("Error: detected unexpected data processing pattern stamp in expression data set file %s. Accepted patterns are 'tstat', 'avg_expr', 'kme'.", file.expr))
  }
  df.expr <- read_expression_data(file.expr=file.expr,
                                  scale_genes=scale_genes,
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
# ======================================== RUN GWAS_LINKING ================================================ #
# ======================================================================================================= #


if (!flag.loaded_gwas_linked_rdata) {
  
  ### Start parallel
  print("Starting registerDoParallel...")
  # cl <- parallel::makeCluster(N_CORES, type="FORK")
  # registerDoParallel(cl) # The registerDoParallel function is used to register the parallel backend with the foreach package.
  registerDoParallel(N_CORES)
  print(sprintf("Number of registered workers (getDoParWorkers()): %s", getDoParWorkers()))
  
  
  
  print("Making first RolyPoly call to link GWAS data...")
  
  ### Run RolyPoly
  runtime.gwas_linked <- system.time(
    rp.gwas_linked <- run_rolypoly_precomputation.multi_attempt(gwas_data = as.data.frame(df.gwas.rp),
                                                                snp_annotations = NULL, # TODO: implement argument to add baseline model as SNP covariates (load via GWAS data)
                                                                block_annotation = as.data.frame(df.block_annotation),
                                                                ld_folder = dir.ldfiles,
                                                                n_attempts=0) # *OBS*: important to make first call with n_attempts=0
  )
  print("Done with RolyPoly GWAS linking. Runtime:")
  print(runtime.gwas_linked)
  print("Printing any warnings from first run:")
  print(warnings())
  
  ### Gene scores
  print("Calculating gene scores...")
  # rp.gwas_linked <- rolypoly_add_ld_corrected_gwas_block_scores(rolypoly=rp.gwas_linked, fast_calculation=T) # returns RolyPoly object | this takes time - it is calculated per gene
    # sets the slot rp$data$<BLOCK/GENE>$corrected_gwas_block_score_pval
  ### You might get the warning:
  # In CompQuadForm::imhof(raw_score, lambda = lambda, epsabs = 1e-15,  ... :
  #                          Note that Qq + abserr is positive.
  print("Done calculating gene scores")
  
  ### Combine into list object
  list.gwas_linked <- list("rp"=rp.gwas_linked, "runtime"=runtime.gwas_linked)
  
  ### Save
  print("Saving...")
  save(list.gwas_linked, list.run_parameters, file=file.out.gwas_linked)
  
  ### stop cluster
  stopImplicitCluster() # stop cluster [stopCluster() should not be used or?]
}


if (GWAS_LINKING_ONLY) {
  print("GWAS_LINKING_ONLY enabled: will quit gracefully now.")
  quit(save = "no", status = 0) # status: the (numerical) error status to be returned to the operating system, where relevant.  Conventionally ‘0’ indicates successful completion.
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
# =========================== CACHING INFERENCE: FILTER EXPRESSION DATA SETS TO RUN ===================== #
# ======================================================================================================= #

list.df_expr.run <- list.df_expr # default value
if (flag.loaded_inference_rdata) { # if we have loaded inference data
  # ASSUMPTION: list.rp_inference and list.df_expr share the same name for inference on expression datasets that have analyzed
  # > names(list.df_expr)
  # [1] "arc_lira/arc_lira.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz"            
  # [2] "gtex/gtex.sub_tissue.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz"  
  
  names.expr_datasets_already_run <- names(list.df_expr) %in% names(list.rp_inference) # return TRUE for dataset names already analyzed
  print(sprintf("*CACHING INFERENCE* flag.loaded_inference_rdata = %s | Found n=%s matching expression data set names with existing RolyPoly inference results.", flag.loaded_inference_rdata, sum(names.expr_datasets_already_run)))
  print("Will not run RolyPoly inference on the following datasets:")
  print(names(list.df_expr)[names.expr_datasets_already_run])
  list.df_expr.run <- list.df_expr[!names.expr_datasets_already_run]
  
  if (length(list.df_expr.run) == 0) { ### check if there are any datasets left to run
    print("*CACHING INFERENCE* | Detected that all datasets in list.df_expr have already been analyzed!")
    print("*CACHING INFERENCE* | Will exit gracefully and not save any new files.")
    quit(save = "no", status = 0) # status: the (numerical) error status to be returned to the operating system, where relevant.  Conventionally ‘0’ indicates successful completion.
  }
}



# ======================================================================================================= #
# ========================================== RUN WRAPPER ================================================ #
# ======================================================================================================= #

print("Running wrapper.run_rolypoly_over_expr_datasets with the following (new) expression data sets:")
print(names(list.df_expr.run))


### Run wrapper - parallel anno
tmp.runtime <- system.time(list.rp_inference.run <- wrapper.run_rolypoly_univariate_over_expr_datasets(list.df_expr = list.df_expr.run,
                                                                                                       rp.precomputed = rp.gwas_linked,
                                                                                                       do.parallel_annotations=TRUE,
                                                                                                       do.parallel_bootstrap=TRUE))
### Run wrapper - NO parallel anno
# tmp.runtime <- system.time(list.rp_inference.run <- wrapper.run_rolypoly_univariate_over_expr_datasets(list.df_expr = list.df_expr.run, 
#                                                                                                        rp.precomputed = rp.gwas_linked, 
#                                                                                                        do.parallel_annotations=FALSE,
#                                                                                                        do.parallel_bootstrap=FALSE))

print("Run wrapper - parallel anno")
print(tmp.runtime)

print("Done with RolyPoly wrapper loop")
print("Printing any warnings:")
print(warnings())

# ======================================================================================================= #
# ================================ CACHING INFERENCE: ADDING NEW DATA TO OLD ============================ #
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

print("FINAL EXPORT | Saving final inference object...")
print(sprintf("Size list.rp_inference: %s", format(object.size(list.rp_inference),units="MB")))
runtime <- system.time(save(list.rp_inference, list.run_parameters, file=file.out.inference))
print(sprintf("Time saving inference object: %s", file.out.inference))
print(runtime)

# ======================================================================================================= #
# ============================================== FINISH ================================================== #
# ======================================================================================================= #

print("SCRIPT DONE!")



# ======================================================================================================= #
# ==========================================   LEFT OVERS   ============================================ #
# ======================================================================================================= #
