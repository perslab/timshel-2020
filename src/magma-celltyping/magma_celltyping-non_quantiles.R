############### SYNOPSIS ###################
# MAGMA Celltyping prototype - using NON-QUANTILES expression models

### *OBS*
# calculate_celltype_associations() and calculate_conditional_celltype_associations() are run with 'SPECIES = human' because the input 'quantiles' have human gene symbols

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:



### USAGE: 
# time Rscript magma_celltyping-non_quantiles.R BMI_Yengo2018

# ======================================================================= #
# ============================  CMD LINE ARGS  ========================== #
# ======================================================================= #

GWAS_NAME <- commandArgs(trailingOnly=TRUE)[1]
print(sprintf("RUNNING GWAS_NAME=%s", GWAS_NAME))

# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

# wd <- "/projects/timshel/sc-genetics/sc-genetics/src/nb-magma_celltyping"
# setwd(wd)


### Install packages
# devtools::install_github(repo = 'NathanSkene/MAGMA_Celltyping') # main package
# devtools::install_github(repo = 'NathanSkene/EWCE')

library(MAGMA.Celltyping) # loads EWCE package
library(tidyverse)


# ======================================================================= #
# ============================  CONSTANTS  =============================== #
# ======================================================================= #

genome_ref_path = "/tools/magma/1.06/g1000_eur/g1000_eur" # should include the 'prefix' plink files

# Load the mouse to human 1:1 orthologs
data(ortholog_data_Mouse_Human)  # loaded from the EWCT package [loads ortholog_data_Mouse_Human]


# ======================================================================= #
#=============================  GWAS DATA  =============================== #
# ======================================================================= #

### SELF FORMATED (SCRATCH)
gwas_sumstats_path <- file.path("/scratch/tmp-magma_gwas/", sprintf("%s.txt", GWAS_NAME)) # not a gziped file

### LDSC PATH
# gwas_sumstats_path <- file.path("/projects/timshel/sc-genetics/sc-genetics/data/gwas_magma_ldsc", sprintf("%s.magma_fmt.txt", GWAS_NAME)) # not a gziped file

### GWAS MANUAL
# gwas_sumstats_path <- "/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/BMI_HEIGHT_Yengo2018/bmi_yengo2018_magma_fmt.txt" # WORKS

# gwas_sumstats_path <- "/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/BMI_HEIGHT_Yengo2018/Meta-analysis_Locke_et_al+UKBiobank_2018.txt"
# ERROR - reading p-value file: non-numeric or non-integer value for sample size variable N on line 447477
# 	line: 6	6176144	rs12210959	T	C	0.7384	-0.0014	0.0022	 5.1e-01	6e+05


# ======================================================================= #
#==========================  MAGMA - ANNOTATE  ============================= #
# ======================================================================= #

## ------------------------------------------------------------------------
#### Format GWAS data (i.e. column headers etc)
# source("magma_celltyping_modified_functions.R")
# col_headers = format_sumstats_for_magma.PT(gwas_sumstats_path) # ---> NEEDS FIXING. DOES NOT WORK YET. (It is also slow.)
## ------------------------------------------------------------------------

source("magma_celltyping_modified_functions.R")
genesOutPath = map.snps.to.genes.PT(gwas_sumstats_path,genome_ref_path=genome_ref_path)


# ======================================================================= #
#===========================  LOAD EXPRESSION DATA  ============================ #
# ======================================================================= #

### read gene mapping for human ENSG --> gene symbol.
DIR.gene_mapping <- "/projects/timshel/git/timshel-lib/bioinformatics/gene_mapping/gene_mapping.GRCh38.ens_v90/GRCh38.ens_v90.ensembl2gene_name_version.txt.gz"
df.gene_mapping <- read_tsv(DIR.gene_mapping)

### read expression data - MOUSEBRAIN
# DIR.expr_data <- "/projects/timshel/sc-genetics/sc-genetics/data/expression/mousebrain"
# filenames <- list.files(path=DIR.expr_data,  pattern="*.hsapiens_orthologs.csv.gz") # has *HUMAN* ENSEMBL IDs
# list.df_expr <- lapply(file.path(DIR.expr_data, filenames), read_csv)
# filenames_shorten <- stringr::str_match(filenames, "celltype_expr\\.(.*)\\.hsapiens_orthologs.csv.gz")[,2] # e.g. "celltype_expr.specificity_quantiles.hsapiens_orthologs.csv.gz" --> "specificity_quantiles"
# names(list.df_expr) <- paste0("MOUSEBRAIN_EXT_", filenames_shorten)
# names(list.df_expr) # e.g. "MOUSEBRAIN_EXT_avg_expr"

### read expression data - NOVO BULK
# DIR.expr_data <- "/projects/timshel/sc-genetics/sc-genetics/data/expression/novo_bulk"
# filenames <- list.files(path=DIR.expr_data,  pattern="*.hsapiens_orthologs.csv.gz") # has *HUMAN* ENSEMBL IDs
# list.df_expr <- lapply(file.path(DIR.expr_data, filenames), read_csv)
# filenames_shorten <- stringr::str_match(filenames, "celltype_expr\\.(.*)\\.hsapiens_orthologs.csv.gz")[,2] # e.g. "celltype_expr.specificity_quantiles.hsapiens_orthologs.csv.gz" --> "specificity_quantiles"
# names(list.df_expr) <- paste0("novo_bulk_", filenames_shorten)
# names(list.df_expr) # e.g. xxx

### read expression data - WGCNA
DIR.expr_data <- "/projects/timshel/sc-genetics/sc-genetics/data/expression/wgcna/kme_combined"
filenames <- list.files(path=DIR.expr_data,  pattern="*.hsapiens_orthologs.csv.gz") # has *HUMAN* ENSEMBL IDs
list.df_expr <- lapply(file.path(DIR.expr_data, filenames), read_csv)
filenames_shorten <- stringr::str_match(filenames, "(.*)\\.hsapiens_orthologs.csv.gz")[,2] # e.g. "maca_tissue_cell_type.kme.kme_na_replaced.hsapiens_orthologs.csv.gz"
names(list.df_expr) <- paste0("kme.", filenames_shorten)
names(list.df_expr) # e.g. xxx


### read expression data - ROSA
# DIR.expr_data <- "/projects/timshel/sc-genetics/sc-genetics/data/expression/rosa"
# filenames <- list.files(path=DIR.expr_data,  pattern="*.hsapiens_orthologs.csv.gz") # has *HUMAN* ENSEMBL IDs
# list.df_expr <- lapply(file.path(DIR.expr_data, filenames), read_csv)
# filenames_shorten <- stringr::str_match(filenames, "celltype_expr\\.(.*)\\.hsapiens_orthologs.csv.gz")[,2] # e.g. "celltype_expr.specificity_quantiles.hsapiens_orthologs.csv.gz" --> "specificity_quantiles"
# names(list.df_expr) <- paste0("rosa_flu_hypo_", filenames_shorten)
# names(list.df_expr) # e.g. xxx


### [***NOT FINISHED***] read expression data - HYPOTHALAMUS z_score_pos_log 
# DIR.expr_data <- "/projects/timshel/sc-genetics/sc-genetics/data/expression/hypothalamus"
# filenames <- list.files(path=DIR.expr_data,  pattern="XXXXXXXz_score_pos_log.hsapiens_orthologs.csv.gz") # has *HUMAN* ENSEMBL IDs
# list.df_expr <- lapply(file.path(DIR.expr_data, filenames), read_csv)
# filenames_shorten <- stringr::str_match(filenames, "celltype_expr\\.(.*)\\.hsapiens_orthologs.csv.gz")[,2] # e.g. "celltype_expr.specificity_quantiles.hsapiens_orthologs.csv.gz" --> "specificity_quantiles"
# names(list.df_expr) <- paste0("MOUSEBRAIN_EXT_", filenames_shorten)
# names(list.df_expr) # e.g. "MOUSEBRAIN_EXT_avg_expr"



### Map human ENSG --> human gene symbols
list.df_expr <- lapply(list.df_expr, function(df) {
  df <- df %>% mutate(gene_symbol = df.gene_mapping$gene_name_optimal[match(df$gene, df.gene_mapping$ensembl_gene_id)]) %>% select(gene_symbol, everything()) # map genes
  print(sprintf("Genes that could not be mapped: %s", sum(is.na(df$gene_symbol))))
  df <- df %>% filter(!is.na(gene_symbol)) # discard non-mapped genes
  df <- df %>% select(-gene) %>% column_to_rownames(var="gene_symbol") # drop ENSG gene names column and SET ROWNAMES
})
### Testing function
# df <- list.df_expr[[1]]
# df <- df %>% mutate(gene_symbol = df.gene_mapping$gene_name_optimal[match(df$gene, df.gene_mapping$ensembl_gene_id)]) %>% select(gene_symbol, everything())
# df.x <- df %>% filter(!is.na(gene_symbol))
# sum(is.na(df$gene_symbol)) # --> 109 genes could not be mapped
# SANITY CHECK: tail(df %>% select(gene_symbol, gene)) DOLPP1 --> ENSG00000167130 [Google says it OK!]


for (DATA_SET in names(list.df_expr)) {
  
  # ======================================================================= #
  # ==========================  PREP EXPR DATA  =========================== #
  # ======================================================================= #
  
  
  # str(ctd) 
  # list[<data_set/annot_level>] --> list of 3 elements (or 4 elements if 'quantiles' has been set using prepare.quantile.groups())
  # list[<data_set/annot_level>]["mean_exp"] # data.frame. genes x cell_types (has rownames with mouse genes)
  # list[<data_set/annot_level>]["specificity"] # matrix. genes x cell_types (has rownames with mouse genes)
  # list[<data_set/annot_level>]["annot"] # character. cell labels
  # list[<data_set/annot_level>]["quantiles"] # matrix. genes x cell_types (has rownames with mouse genes)
  
  df.data_in <- list.df_expr[[DATA_SET]] # data.frame
  
  ### All the 'calculate_celltype_associations' function requires for the ctd object, is that it has an '$quantiles' element
  ctd = list(list(  # OBS: Remember that ctd is list of lists.
    "mean_exp"="DUMMY",
    "annot"="DUMMY",
    "specificity"="DUMMY",
    "quantiles"=as.matrix(df.data_in)
  ))
  
  # ======================================================================= #
  #==========================  MAGMA - MAIN  ============================= #
  # ======================================================================= #
  
  ctAssocs = calculate_celltype_associations(ctd,gwas_sumstats_path,genome_ref_path=genome_ref_path, specificity_species="human") # *OBS*: SPECIES = human
  # this function calls create_gene_covar_file() which MAPS genes from mouse to human if 'specificity_species="mouse"'. (Also maps gene symbol to Entrez ID)
  
  ### Extract results
  df.results <- ctAssocs[[1]][["results"]] # create data frame
  df.results %>% arrange(P) %>% head()
  df.results %>% arrange(P) %>% write_csv(sprintf("out.magma_celltyping.%s.results.%s.csv", GWAS_NAME, DATA_SET))
  
  ## ------------------------------------------------------------------------
  plot_celltype_associations(ctAssocs, savePDF=F)
  ggsave(sprintf("out.magma_celltyping.%s.plot.%s.pdf", GWAS_NAME, DATA_SET), w=14, h=40)
  
  # ======================================================================= #
  #=========================  MAGMA - CONDITIONAL  ======================== #
  # ======================================================================= #
  
  # Error in calculate_conditional_celltype_associations(ctd, gwas_sumstats_path,  :
  # No celltypes reach significance with Q<0.05
  
  tryCatch({
    ## ------------------------------------------------------------------------
    ### Conditional analysis
    ctCondAssocs = calculate_conditional_celltype_associations(ctd,gwas_sumstats_path,genome_ref_path=genome_ref_path,controlTopNcells=3,specificity_species="human") # *OBS*: SPECIES = human
    
    ### Try this
    df.results_cond <- ctCondAssocs[[1]][["results"]] # create data frame
    df.results_cond %>% arrange(P) %>% head()
    df.results_cond %>% distinct(CONTROL_label)
    
    plot_celltype_associations(ctCondAssocs, savePDF=F)
    ggsave(sprintf("out.magma_celltyping.%s.plot_cond.%s.pdf", GWAS_NAME, DATA_SET), w=20, h=40)
  },error=function(cond) {
    print(sprintf("There was an error running calculate_conditional_celltype_associations: %s", cond))
  })
  
  
  # ======================================================================= #
  #=========================  SAVE DATA SESSION =========================== #
  # ======================================================================= #
  
  save.image(file=sprintf("rsession_celltyping.%s.%s.RData", GWAS_NAME, DATA_SET))
  
  # load(file=sprintf("rsession_celltyping.%s.%s.RData", GWAS_NAME, DATA_SET))
  
}  

print("SCRIPT DONE")



