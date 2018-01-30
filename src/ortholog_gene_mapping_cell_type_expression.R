############### SYNOPSIS ###################
# Transform expression : map genes and XXX

### DESCRIPTION
# INPUT: avg expression file + mouse MGI-Ensembl mapping + ortholog mapping file.
# 1. Load avg expression file.
# 4. Map MOUSE GENE NAME to ENSEMBL: Map expression matrices MGI gene names to mouse Ensembl IDs
#     1. use timshel-lib from gene name to ensembl: "Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz"
# 5. SUBSET genes on 1-1 orthologs: in the expression file, keep only Ensembl IDs listed in "1-1 ortholog file"
# 6. REPLACE mouse ID with human ID.

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

############################################
library(Seurat)
library(tidyverse)


rm(list=ls())

wd <- "/Users/djw472/Dropbox/0_Projects/p_sc_genetics/analysis/src"
setwd(wd)


# ======================================================================= #
# ==============================  READ DATA  ============================== #
# ======================================================================= #

file.expr <- "../data/expression/arc_lira-celltype_expr_ttest.csv"
# df.expr <- read_csv(file.expr)
df.expr <- read.csv(file.expr)


file.gene_name_mapping <- "../data/gene_annotations/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz"
df.gene_name_mapping <- read_delim(file.gene_name_mapping, delim="\t")

file.ortholog_mapping <- "../data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz"
df.ortholog_mapping <- read_delim(file.ortholog_mapping, delim="\t")

# ======================================================================= #
# =========  EXPRESSION DATA: NUMERICAL TRANSFORMATION  ================= #
# ======================================================================= #

# df.expr.trans <- df.expr^2

# ======================================================================= #
# ==========================  Map: MGI to EnsemblID  ==================== #
# ======================================================================= #

### Get gene names
genes.name <- rownames(df.expr)

### Genes not found
bool.not_found <- !(genes.name %in% df.gene_name_mapping$gene_name_optimal) 
sum(bool.not_found) # n genes could not be mapped --> 607
genes.name[bool.not_found]
# e.g Fam150a not found, but this is because it is called ALKAL1. Fam150a is listed as a synonym in Emsembl now.

### Match
idx.match <- match(genes.name, df.gene_name_mapping$gene_name_optimal) # find mapping idx. 'nomatch' is set to NA
head(idx.match)
genes.ensembl <- df.gene_name_mapping$ensembl_gene_id[idx.match] # vector with same length as "genes.name", but we NA values for genes not found.
head(genes.ensembl)
# df.tmp <- data.frame(genes.name, genes.ensembl)

### Filter out genes not mapped, and replace gene names with IDs.
df.expr.ens <- df.expr[!is.na(genes.ensembl),] # subset to mapped genes
rownames(df.expr.ens) <- genes.ensembl[!is.na(genes.ensembl)] # replace gene names with IDs

head(df.expr.ens)



# ======================================================================= #
# ========================  Subset on orthologs  ======================== #
# ======================================================================= #

### Find mouse genes with human orthologs
idx.match <- match(rownames(df.expr.ens), df.ortholog_mapping$mmusculus_homolog_ensembl_gene)
genes.human.orthologs <- df.ortholog_mapping$ensembl_gene_id[idx.match] # NA values for genes with no ortholog.

### Genes not found
bool.not_found <- !(rownames(df.expr.ens) %in% df.ortholog_mapping$mmusculus_homolog_ensembl_gene) 
sum(bool.not_found) # genes with no ortholog --> 4853


### Subset genes
df.expr.ens.human <- df.expr.ens[!is.na(genes.human.orthologs),] # subset to ortholog genes
rownames(df.expr.ens.human) <- genes.human.orthologs[!is.na(genes.human.orthologs)] # replace mouse ensembl ID with human

head(df.expr.ens.human)

# ======================================================================= #
# ===============================  EXPORT  ============================== #
# ======================================================================= #

### Write table
file.out <- "arc_lira-celltype_expr_ttest.hsapiens_orthologs.csv.gz"
file.out.gz <- gzfile(file.out, 'w')
write.table(df.expr.ens.human, file=file.out.gz, col.names=T, row.names=T, quote=F, sep=",")
close(file.out.gz)





