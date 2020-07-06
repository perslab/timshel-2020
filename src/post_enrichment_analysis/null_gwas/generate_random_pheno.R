############### SYNOPSIS ###################
# Generate random phenotype for "null GWAS".

### DESCRIPTION
# ...

### OUTPUT: 
# ....

### REMARKS:
# ....


############################################

library(here)

# ==================

set.seed(1)

# ==================

### Get number of individuals
file.plink.fam <- here("src","post_enrichment_analysis", "null_gwas", "phenotypes", "ALL.chr_merged.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.fam")

df.plink.fam <- read.table(file.plink.fam, header=F)


### Set params
n.individuals <- nrow(df.plink.fam)
n.null_gwas <- 1000

### Generate matrix
pheno <- matrix(replicate(n.null_gwas, rnorm(n.individuals)), ncol=n.null_gwas)

# ==================
### Create data.frame
df.pheno <- data.frame(familyID=df.plink.fam[,1], within_familyID=df.plink.fam[,2], pheno)
df.pheno.format <- format(df.pheno, digits=1, nsmall=2)

### Check
stopifnot(df.pheno[1,3]+0.6264538<= 1e-6) # a small check for the first element, to ensure that the matrix does not change.
# ==================


### Write file

file.out <- here("src","post_enrichment_analysis", "null_gwas", "phenotypes",sprintf("null_gwas_guassian_phenotypes.n%s.phe", n.null_gwas))
write.table(df.pheno.format, file=file.out, sep="\t", quote=F, col.names=F, row.names=F)

# ==================

### .pphe file description (used for --pheno)
# A text file with no header line, and one line per sample with the following P+2 fields (where P is the requested number of permutations):
# 1. Family ID
# 2. Within-family ID
# 3-(P+2). Permuted phenotypes





