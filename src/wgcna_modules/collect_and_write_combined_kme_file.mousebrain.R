############### SYNOPSIS ###################
# Combine kME files for Mousebrain from each cell-type
# E.g. the kME file for the "PER3" cell-type is /projects/jonatan/tmp-mousebrain/tables/mousebrain_Vascular_ClusterName_2_PER3_kME.csv

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# =============================== SETUP ================================ #
# ======================================================================= #

library(tidyverse)

wd <- "/projects/timshel/sc-genetics/sc-genetics/src/wgcna_modules/"
setwd(wd)

# ======================================================================= #
# =============================== READ kME FILES ================================ #
# ======================================================================= #

### Read
dir.data <- "/projects/jonatan/tmp-mousebrain/tables"
filenames <- list.files(path=dir.data,  pattern="mousebrain_(.*)_ClusterName_2_(.*)_kME.csv")
list.df <- lapply(file.path(dir.data, filenames), data.table::fread, nThread=24, data.table=F, showProgress=T)
filenames_shorten <- stringr::str_match(filenames, "mousebrain_.*_ClusterName_2_(.*)_kME.csv")[,2]
names(list.df) <- filenames_shorten

### Add cell-type to columnnames to make them unique
add_cell_type_to_columnames <- function(name.list, list) {
  df <- list[[name.list]]
  colnames(df)[-1] <- paste0(name.list, ".", colnames(df)[-1])
  return(df)
  }
list.df.mapped <- lapply(names(list.df), add_cell_type_to_columnames, list.df)

# list.df.mapped[[1]]

### Join
df.join <- list.df.mapped %>% purrr::reduce(full_join, by = "genes")


# ======================================================================= #
# =============================== WRITE ================================ #
# ======================================================================= #

#data.table::fwrite(df.join, file="/projects/timshel/sc-genetics/sc-genetics/out/module_correlation_graph/mousebrain.all_modules_kME.csv", sep=",", nThread=24, showProgress=T)

# ======================================================================= #
# ======================== Calculate correlation ======================== #
# ======================================================================= #

system.time(df.cor.pearson <- df.join %>% select(-genes) %>% cor(use = "pairwise.complete.obs", method="pearson") %>% as.data.frame())
system.time(df.cor.spearman <- df.join %>% select(-genes) %>% cor(use = "pairwise.complete.obs", method="spearman") %>% as.data.frame())
file.out.pearson <- "/projects/timshel/sc-genetics/sc-genetics/out/module_correlation_graph/mousebrain.all_modules_kME.correlation.pearson.csv"
file.out.spearman <- "/projects/timshel/sc-genetics/sc-genetics/out/module_correlation_graph/mousebrain.all_modules_kME.correlation.spearman.csv"
write_csv(df.cor.pearson, path=file.out.pearson)
write_csv(df.cor.spearman, path=file.out.spearman)


# ======================================================================= #
# =============================== XXXXXX ================================ #
# ======================================================================= #


# ======================================================================= #
# =============================== XXXXXX ================================ #
# ======================================================================= #






