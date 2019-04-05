############### SYNOPSIS ###################
# MAGMA Celltyping prototype

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


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
#===========================  LOAD EWCT DATA  ============================ #
# ======================================================================= #

DATA_SETS <- c("allKI_lvl1_skene", "allKI_lvl2_skene", "MACA", "mousebrain") # "hypothalamus_log_avg", "skene_KI",
for (DATA_SET in DATA_SETS) {
  ### Load the celltype data
  # DATA_SET <- "MACA"
  # DATA_SET <- "mousebrain"
  
  if (DATA_SET == "skene_KI") {
    data(ctd) # loaded from the EWCT package
  } else {
    load(sprintf("data-ctd/CellTypeData.%s.RData", DATA_SET)) # loads ctd file
  }
  
  # str(ctd) 
  # list[<data_set/annot_level>] --> list of 3 elements (or 4 elements if 'quantiles' has been set using prepare.quantile.groups())
  # list[<data_set/annot_level>]["mean_exp"] # data.frame. genes x cell_types (has rownames)
  # list[<data_set/annot_level>]["specificity"] # matrix. genes x cell_types (has rownames)
  # list[<data_set/annot_level>]["annot"] # character. cell labels
  # list[<data_set/annot_level>]["quantiles"] # matrix. genes x cell_types (has rownames)
  # x <- ctd[[1]][["mean_exp"]]
  # x <- ctd[[1]][["specificity"]]
  # x <- ctd[[1]][["annot"]]
  
  # ======================================================================= #
  # ==========================  PREP EXPR DATA  =========================== #
  # ======================================================================= #
  
  ctd = prepare.quantile.groups(ctd,specificity_species="mouse",numberOfBins=40)
  
  ## ------------------------------------------------------------------------
  # Examine how the quantile groups look
  print(ctd[[1]]$quantiles[c("Gfap","Dlg4","Aif1"),])
  print(table(ctd[[1]]$quantiles[,1]))
  
  
  # ======================================================================= #
  #==========================  MAGMA - MAIN  ============================= #
  # ======================================================================= #
  
  ctAssocs = calculate_celltype_associations(ctd,gwas_sumstats_path,genome_ref_path=genome_ref_path, specificity_species="mouse")
  
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
    ctCondAssocs = calculate_conditional_celltype_associations(ctd,gwas_sumstats_path,genome_ref_path=genome_ref_path,controlTopNcells=3)
    
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



