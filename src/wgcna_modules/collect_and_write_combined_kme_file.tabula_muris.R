############### SYNOPSIS ###################
# Combine kME files for Tabula Muris from each cell-type
# E.g. "/projects/jonatan/tmp-maca/tables/maca_tissue_cell_type_kME_Trachea_stromal cell_kME.csv"

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
dir.data <- "/projects/jonatan/tmp-maca/tables"
filenames <- list.files(path=dir.data,  pattern="maca_tissue_cell_type_kME_(.*)_kME.csv")
list.df <- lapply(file.path(dir.data, filenames), data.table::fread, nThread=24, data.table=F, showProgress=T)
filenames_shorten <- stringr::str_match(filenames, "maca_tissue_cell_type_kME_(.*)_kME.csv")[,2]
names(list.df) <- filenames_shorten

### Add cell-type to columnnames to make them unique
### ---> SKIP THIS STEP FOR tabula muris BECAUSE MODULE IDs are UNIQUE (1 WGCNA run)
# add_cell_type_to_columnames <- function(name.list, list) {
#   df <- list[[name.list]]
#   colnames(df)[-1] <- paste0(name.list, ".", colnames(df)[-1])
#   return(df)
# }
# list.df.mapped <- lapply(names(list.df), add_cell_type_to_columnames, list.df)

list.df.mapped <- list.df # stupid copy but what the heck...
### Join
df.join <- list.df.mapped %>% purrr::reduce(full_join, by = "genes")


# ======================================================================= #
# =============================== WRITE ================================ #
# ======================================================================= #

# data.table::fwrite(df.join, file="/projects/timshel/sc-genetics/sc-genetics/out/module_correlation_graph/tabula_muris.all_modules_kME.csv", sep=",", nThread=24, showProgress=T)

# ======================================================================= #
# ======================== Calculate correlation ======================== #
# ======================================================================= #

print("Doing correlation pearson")
system.time(df.cor.pearson <- df.join %>% select(-genes) %>% cor(use = "pairwise.complete.obs", method="pearson") %>% as.data.frame())
print("Doing correlation spearman")
system.time(df.cor.spearman <- df.join %>% select(-genes) %>% cor(use = "pairwise.complete.obs", method="spearman") %>% as.data.frame())
file.out.pearson <- "/projects/timshel/sc-genetics/sc-genetics/out/module_correlation_graph/tabula_muris.all_modules_kME.correlation.pearson.csv"
file.out.spearman <- "/projects/timshel/sc-genetics/sc-genetics/out/module_correlation_graph/tabula_muris.all_modules_kME.correlation.spearman.csv"
write_csv(df.cor.pearson, path=file.out.pearson)
write_csv(df.cor.spearman, path=file.out.spearman)


# [1] "Doing correlation pearson"
#     user   system  elapsed
# 1318.126    3.507 1321.642 ---> 21 min
# [1] "Doing correlation spearman"
#     user   system  elapsed
# 6867.744   53.667 6921.549 ---> 2 h

# ======================================================================= #
# =============================== XXXXXX ================================ #
# ======================================================================= #


# ======================================================================= #
# =============================== XXXXXX ================================ #
# ======================================================================= #






