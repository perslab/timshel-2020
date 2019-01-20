############### SYNOPSIS ###################
# Process WGCNA kme files: MOUSEBRAIN AND MACA.
# THIS SCRIPT IS KIND OF TEMPORARY, SO MOVE CONTENT ELSEWHERE

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

library(tidyverse)

### Source custom scripts
dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib/"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library


wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-wgcna_modules/"
setwd(wd)


# ======================================================================= #
# ==========================  read kME files  ============================ #
# ======================================================================= #

### MACA
# DATA_SET_NAME <- "maca_tissue_cell_type.kme" # ONLY USED FOR OUTPUT FILES
# DIR.expr_data <- "/raid5/projects/jonatan/tmp-maca/tables"
# filenames <- list.files(path=DIR.expr_data,  pattern="maca_tissue_cell_type_kME_(.*)_kME.csv") # maca_tissue_cell_type_kME_Tongue_keratinocyte_kME.csv
# list.df_expr <- lapply(file.path(DIR.expr_data, filenames), read_csv)
# filenames_shorten <- stringr::str_match(filenames, "maca_tissue_cell_type_kME_(.*)_kME.csv")[,2]
# names(list.df_expr) <- filenames_shorten
# names(list.df_expr) # e.g. xxx

### MOUSEBRAIN neurons
DATA_SET_NAME <- "mousebrain_Neurons_ClusterName.kme"  # ONLY USED FOR OUTPUT FILES
DIR.expr_data <- "/raid5/projects/jonatan/tmp-mousebrain/tables"
filenames <- list.files(path=DIR.expr_data,  pattern="mousebrain_Neurons_ClusterName_2_(.*)_kME.csv") # mousebrain_Vascular_ClusterName_2_ENMFB_kME.csv
list.df_expr <- lapply(file.path(DIR.expr_data, filenames), read_csv)
filenames_shorten <- stringr::str_match(filenames, "mousebrain_Neurons_ClusterName_2_(.*)_kME.csv")[,2]
names(list.df_expr) <- filenames_shorten
names(list.df_expr) # e.g. xxx

# ======================================================================= #
# ==========================  Merge, Map Orthologs and Set NA to Zero  ============================ #
# ======================================================================= #

### merge | *OBS* there will be MANY 'NA' values after joning
df.kme <- list.df_expr %>% purrr::reduce(full_join, by = "genes")
# mousebrain neurons --> 18647 x 3833
# maca tissue cell-type --> 17695 x 3100

### Ortholog mapping
df.human <- mouse_to_human_ortholog_gene_expression_mapping(df.kme %>% column_to_rownames(var="genes"), type_mouse_gene_ids="ensembl") #  input data frame *MUST have rownames* with *mouse genes*.

### Set NA values to zero. This is a fair thing to do: the genes are not part of the pathway and should not count
### MAGMA will throw away covariates if too many of them have NA values
df.human[is.na(df.human)] <- 0
### ^ dplyrr::replace_na(0) is not helpful here because it requires a named list.
any(is.na(df.human)) # ---> False


# ======================================================================= #
# ======================= WRITE OUT ========================== #
# ======================================================================= #

file.out.kme.human <- sprintf("%s.kme_na_replaced.hsapiens_orthologs.csv.gz", DATA_SET_NAME)
write_csv(df.human %>% rownames_to_column(var="gene"), path=file.out.kme.human)

# ======================================================================= #
# ======================= SEM models ========================== #
# ======================================================================= #

source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library
# ======= Skene ======
df.specificity.quantiles <- wrapper_specificity_quantiles(df.human)
file.out.kme.human <- sprintf("%s.specificity_quantiles.hsapiens_orthologs.csv.gz", DATA_SET_NAME)
write_csv(df.specificity.quantiles %>% rownames_to_column(var="gene"), path=file.out.kme.human)




