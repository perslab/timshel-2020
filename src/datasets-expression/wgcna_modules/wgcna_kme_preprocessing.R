############### SYNOPSIS ###################
# Process WGCNA kme files: ortholog mapping

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# ==============================  USAGE  =============================== #
# ======================================================================= #

# time Rscript wgcna_kme_preprocessing.R campbell_agrp /projects/tp/tmp-bmi-brain/data/wgcna/campbell_arc/agrp_kmemedule.csv
# time Rscript wgcna_kme_preprocessing.R campbell_pomc /projects/tp/tmp-bmi-brain/data/wgcna/campbell_arc/pomc_kmemedule.csv

### AGRP and POMC (the same output)
# [1] "Number of genes not mapped from MGI to Ensembl ID: 118 out of 5000 genes (2.36 pct)"
# [1] "Number of genes not mapped to human ortholog: 483 out of 4882 genes (9.89 pct)"


# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

library(tidyverse)
library(here)

library(ggplot2)
library(GGally)


### Source custom scripts
dir.project_src <- "/projects/timshel/sc-arc_lira/src"
dir.pers_lab_sc_lib <- "/projects/timshel/git/perslab-sc-library"
source(sprintf("%s/seurat_functions/load_functions.R", dir.pers_lab_sc_lib)) # load Pers lab/Timshel single-cell library

source(here("src/lib/load_functions.R")) # load sc-genetics library


# ======================================================================= #
# ============================  PARAMS  ============================== #
# ======================================================================= #

wd <- "/projects/timshel/sc-genetics/sc-genetics/src/GE-wgcna_modules/"
setwd(wd)

data_set_name <- commandArgs(trailingOnly=TRUE)[1]
file.data <- commandArgs(trailingOnly=TRUE)[2]
print(sprintf("RUNNING file.data=%s", file.data))


# file.data <- "/projects/tp/tmp-bmi-brain/data/wgcna/campbell_arc/agrp_kmemedule.csv"
# data_set_name <- "campbell_agrp"
# file.data <- "/projects/tp/tmp-bmi-brain/data/wgcna/campbell_arc/pomc_kmemedule.csv"
# data_set_name <- "campbell_pomc"

# ======================================================================= #
# ============================  INITIALIZE  ============================== #
# ======================================================================= #

### Load data

df <- read_csv(file.data)
df <- df %>% as.data.frame() %>% column_to_rownames(var="X1")

# ======================================================================= #
# ================================  PARAMs  ============================= #
# ======================================================================= #


### Output files
file.out.kme <- sprintf("%s.celltype_expr.kme.csv.gz", data_set_name)
file.out.kme.human <- sprintf("%s.celltype_expr.kme.hsapiens_orthologs.csv.gz", data_set_name)


# ======================================================================= #
# ==========================  Map Orthologs  ============================ #
# ======================================================================= #

### Ortholog
df.human <- mouse_to_human_ortholog_gene_expression_mapping(df) #  input data frame *MUST have rownames* with *MGI symbols*.
write_csv(df.human %>% rownames_to_column(var="gene"), path=file.out.kme.human) # write_* automatically deals with .gz, .bz2 or .xz output filenames

# ======================================================================= #
# =====================  [OPTIONAL - look at kme values]  ============== #
# ======================================================================= #

### correlation
p <- ggcorr(df, geom="circle", label = TRUE, label_size = 1.5, label_round = 2, label_alpha = TRUE, hjust = 1, size = 2, layout.exp = 5)
p <- p + labs(title=sprintf("module correlation: %s", data_set_name))
ggplot2::ggsave(filename=sprintf("plot.module_correlation.%s.pdf", data_set_name), plot=p, width=10, height=10)

### density
p <- ggplot(df %>% gather(key="module", value="gene_kme_value"), aes(x=gene_kme_value, fill=module)) + geom_density(alpha=0.4)
p <- p + labs(title=sprintf("dataset: %s", data_set_name))
ggplot2::ggsave(filename=sprintf("plot.kme_distribution.%s.pdf", data_set_name), plot=p, width=10, height=10)


# ======================================================================= #
# =======================      FINISH       ========================== #
# ======================================================================= #

print("Script done!")

