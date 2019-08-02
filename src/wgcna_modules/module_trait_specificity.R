############### SYNOPSIS ###################
# Module "trait specificity": LDSC genetic prioritization for selected modules across multiple traits


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

source(here("src/lib/load_functions.R")) # load sc-genetics library

setwd(here("src/wgcna_modules"))

# ======================================================================= #
# ============================ PARAMETERS =============================== #
# ======================================================================= #

prefix_genomic_annot <- "wgcna.mousebrain-190213.fdr_sign_celltypes.continuous"

# ======================================================================= #
# ============================ FUNCTIONS =============================== #
# ======================================================================= #


load_ldsc_cts_results_wgcna <- function(file.ldsc_cts) {
  df.ldsc_cts <- read_tsv(file.ldsc_cts)
  mat_tmp_split_str <- stringr::str_split_fixed(df.ldsc_cts$Name,pattern="\\__",n=Inf)
  df.ldsc_cts <- df.ldsc_cts %>% mutate(module_id=mat_tmp_split_str[,ncol(mat_tmp_split_str)]) # add module ID from 'Name' string, e.g. maca_tissue_cell_type.slateblue4
  df.ldsc_cts <- df.ldsc_cts %>% rename(BETA=Coefficient, SE=Coefficient_std_error, P=Coefficient_P_value)
  return(df.ldsc_cts)
}

# ======================================================================= #
# ======================== LOAD LDSC CTS RESULTS ======================== #
# ======================================================================= #


### Load - MULTI GWAS
dir.data <- "/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/"
filenames <- list.files(path=dir.data,  pattern=sprintf("%s.(.*).cell_type_results.txt", prefix_genomic_annot))
filenames <- filenames[!grepl(pattern="__CONDITIONAL__", filenames, perl=T)] # exclude any conditional results
list.dfs <- lapply(file.path(dir.data, filenames), load_ldsc_cts_results_wgcna)
names(list.dfs) <- stringr::str_match(filenames, pattern=sprintf("%s__(.*).cell_type_results.txt", prefix_genomic_annot))[,2] # ALT: filenames
names(list.dfs)
df.ldsc_cts <- list.dfs %>% bind_rows(.id="gwas")

# ======================================================================= #
# ================================ EXPORT =============================== #
# ======================================================================= #

df.ldsc_cts <- df.ldsc_cts %>% filter(module_id %in% c("lavenderblush", "lightpink3")) # filter modules

### format file [*not* compatible with multi-GWAS plotting]
df.ldsc_cts.export <- df.ldsc_cts %>% select(gwas, P, module_id) %>% spread(key=module_id, value=P) # format file
df.ldsc_cts.export <- df.ldsc_cts.export %>% arrange(lavenderblush) %>% mutate(gwas=factor(gwas, levels=gwas)) # order factor
# file.out <- here("results/", sprintf("prioritization_modules.%s.multi_gwas.csv.gz", prefix_genomic_annot))
file.out <- here("results/", "prioritization_modules--fdr_sign_celltypes.multi_gwas.csv.gz")

# df.ldsc_cts.export %>% write_csv(file.out)


