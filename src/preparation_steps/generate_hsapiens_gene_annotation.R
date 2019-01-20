############### SYNOPSIS ###################
# Get biomaRt annotation

### DESCRIPTION
# ...

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:
# Ensembl biomaRt guide: http://www.ensembl.org/info/data/biomart/biomart_r_package.html

############################################

### Installation
# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
library(biomaRt)
library(tidyverse)


# rm(list=ls())

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src"
setwd(wd)


# ======================================================================= #
# ==============================  biomaRt  ============================== #
# ======================================================================= #
### Set version variable
#ENSEMBL_VERSION = NULL # for newest version
ENSEMBL_VERSION = 91 # Ensembl release 91 = December 2017 | for specific version
GENOME_VERSION = 37 # GRCh37=hg19; GRCh38=hg38

### Connect to the BioMart database and dataset hosted by Ensembl
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=ENSEMBL_VERSION, GRCh=GENOME_VERSION, verbose=T)


### Inspect attributes and filters
df.ensembl.attr <- listAttributes(ensembl) # ~1280 attributes
df.ensembl.filters <- listFilters(ensembl) # ~265 filers (e.g. entrezgene)


# ======================================================================= #
# ========================= PAGE = GENES FEATURES  ====================== #
# ======================================================================= #

### Get all Ensembl Gene IDs - JUST TO SEE HOW MANY GENES THERE ARE IN TOTAL
df.BM.all.ids <- getBM(attributes = c("ensembl_gene_id"), mart=ensembl)
nrow(df.BM.all.ids) # 66203 for ensembl_v84
str(df.BM.all.ids)

### Features page
df.BM.feature <- getBM(attributes = c("ensembl_gene_id", 
                                      "chromosome_name",
                                      "start_position", # Gene start (bp)
                                      "end_position", # Gene end (bp)
                                      "gene_biotype"
), mart=ensembl)
str(df.BM.feature)
nrow(df.BM.feature) 

length(unique(df.BM.feature$ensembl_gene_id)) == nrow(df.BM.feature) # --> TRUE, all ensembl_gene_id are unique

# ======================================================================= #
# ===============================  EXPORT  ============================== #
# ======================================================================= #


df.BM.feature %>% count(chromosome_name) %>% print(n=50)
# ---> does contain X and Y chromosomes
df.BM.feature %>% distinct(chromosome_name) %>% nrow() # ---> 265 different chromosomes

# ======================================================================= #
# ===============================  EXPORT  ============================== #
# ======================================================================= #

### Write table
file.ensmbl_annotation <- "gene_annotation.hsapiens_all_genes.GRCh37.ens_v91.txt.gz"
write_tsv(df.BM.feature, path=file.ensmbl_annotation)



