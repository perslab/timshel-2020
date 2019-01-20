### SYNOPSIS: fit MAGMA ZSTATs to SEM models

library(tidyverse)

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/nb-sem/"
setwd(wd)

# ======================================================================= #
# =============================== FUNCTIONS ================================ #
# ======================================================================= #

### Map:  ENTREZ --> ENSEMBL
map_entrez2ensembl <- function(df.magma) {
  file.mapping <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_id_mapping.hsapiens.ensembl_entrez.txt.gz"
  df.mapping <- read_tsv(file.mapping)
  
  genes_mapped <- df.mapping$ensembl_gene_id[match(df.magma$GENE, df.mapping$entrezgene)]
  print(sprintf("Number of genes mapped: %s",sum(!is.na(genes_mapped))))
  print(sprintf("Number of genes not mapped: %s",sum(is.na(genes_mapped)))) # number of not mapped genes
  df.magma.clean <- df.magma %>% 
    mutate(ensembl_gene_id=genes_mapped) %>%
    filter(!is.na(ensembl_gene_id))
  return(df.magma.clean)
}

# ========================================================================== #
# =============================== READ MAGMA GWAS  ======================== #
# ========================================================================== #

### MAGMA GWAS
file.magma_gwas <- "/scratch/tmp-magma_gwas/BMI_Yengo2018.txt.10UP.1.5DOWN.genes.out"
### SNIPPIT
## GENE       CHR      START       STOP  NSNPS  NPARAM       N        ZSTAT            P
## 375790       1     945503     992999      1       1  605446            0          0.5
## 401934       1    1005626    1019687      6       1  662343    -0.036852       0.5147
df.magma_gwas.raw <- map_entrez2ensembl(read_table(file.magma_gwas)) # magma annotation window: 10kb up, 1.5kb down
# [1] "Number of genes mapped: 17467"
# [1] "Number of genes not mapped: 158"

### MAGMA GWAS RESIDUALS (genesize.density)
file.magma_gwas.rsid <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/BMI_Yengo2018.baseline_resid.genesize.density.genes.out"
### SNIPPIT (not that first line is outcommented)
## CONDITIONED_INTERNAL = gene size, gene density, log(gene size), log(gene density) 
# GENE       CHR      START       STOP  NSNPS  NPARAM       N        ZSTAT  ZFITTED_BASE  ZRESID_BASE
# 375790       1     945503     992999      1       1  605446      -1.0936   -0.00018461      -1.0934
# 54991        1    1015698    1061736     24       5  649301      -1.0953   -0.00018461      -1.0951
df.magma_gwas.resid <- map_entrez2ensembl(read_table(file.magma_gwas.rsid, skip=1))


### Select
df.magma_gwas <- df.magma_gwas.resid

# ============================================================================ #
# ======================= READ SEM DATA ================ #
# ============================================================================ #

### read SEM expression data
DIR.sem_data <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/mousebrain"
filenames <- list.files(path=DIR.sem_data,  pattern="*.hsapiens_orthologs.csv.gz") # has *HUMAN* ENSEMBL IDs
list.df_sem <- lapply(file.path(DIR.sem_data, filenames), read_csv)
filenames_shorten <- stringr::str_match(filenames, "celltype_expr\\.(.*)\\.hsapiens_orthologs.csv.gz")[,2] # e.g. "celltype_expr.specificity_quantiles.hsapiens_orthologs.csv.gz" --> "specificity_quantiles"
names(list.df_sem) <- filenames_shorten
names(list.df_sem)


# ============================================================================ #
# ======================= TESTING FIT ================ #
# ============================================================================ #

df.sem <- list.df_sem[["specificity_quantiles"]]
# top: HBGLU3
# bottom: VLMC2

### combine data
df.model <- df.sem %>% 
  select(ensembl_gene_id=gene, top=HBGLU3, bottom=VLMC2) %>% 
  left_join(df.magma_gwas %>% select(ensembl_gene_id, ZSTAT), by="ensembl_gene_id")

df.model
### point
ggplot(df.model, aes(x=top, y=ZSTAT)) + geom_point()
### violin
ggplot(df.model, aes(x=top, y=ZSTAT, group=top)) + geom_violin()
ggplot(df.model, aes(x=bottom, y=ZSTAT, group=bottom)) + geom_violin()
### loess
ggplot(df.model, aes(x=top, y=ZSTAT)) + geom_smooth(method = "loess", size = 1.5)
ggplot(df.model, aes(x=top, y=ZSTAT)) + geom_smooth(method = "loess", size = 1.5) + geom_point(alpha=0.2)
ggplot(df.model, aes(x=bottom, y=ZSTAT)) + geom_smooth(method = "loess", size = 1.5)



m <- lm(ZSTAT ~ top, data=df.model)
summary(m)
m <- lm(ZSTAT ~ bottom, data=df.model)
summary(m)

# ============================================================================ #
# ======================= ZSTAT RESID ================ #
# ============================================================================ #
# gene size, gene density, log(gene size), log(gene density) 

# ### TOP
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.038681   0.023874    1.62    0.105    
# top         0.012648   0.001177   10.74   <2e-16 *** ---> less significant than unconditioned ZSTAT
# Residual standard error: 1.899 on 13771 degrees of freedom
# (1172 observations deleted due to missingness)
# Multiple R-squared:  0.008313,	Adjusted R-squared:  0.008241 
# F-statistic: 115.4 on 1 and 13771 DF,  p-value: < 2.2e-16
# 
# 
# ### Bottom
# (Intercept) 0.189421   0.026063   7.268 3.85e-13 ***
#   bottom      0.002240   0.001207   1.856   0.0634 .   ---> slightly more significant than unconditioned ZSTAT
# Residual standard error: 1.907 on 13771 degrees of freedom
# (1172 observations deleted due to missingness)
# Multiple R-squared:  0.0002501,	Adjusted R-squared:  0.0001775 
# F-statistic: 3.445 on 1 and 13771 DF,  p-value: 0.06345


# ============================================================================ #
# ======================= ZSTAT not conditioned on any internal covariates ================ #
# ============================================================================ #

### TOP
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 1.519098   0.024230   62.69   <2e-16 ***
#   top         0.017897   0.001195   14.98   <2e-16 ***
#   Residual standard error: 1.929 on 13787 degrees of freedom
# (1171 observations deleted due to missingness)
# Multiple R-squared:  0.01602,	Adjusted R-squared:  0.01595 
# F-statistic: 224.4 on 1 and 13787 DF,  p-value: < 2.2e-16

### Bottom
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.81025    0.02656  68.147   <2e-16 ***
#   bottom      -0.00144    0.00123  -1.171    0.242    
# Residual standard error: 1.944 on 13787 degrees of freedom
# (1171 observations deleted due to missingness)
# Multiple R-squared:  9.938e-05,	Adjusted R-squared:  2.686e-05 
# F-statistic:  1.37 on 1 and 13787 DF,  p-value: 0.2418

# ============================================================================ #
# ======================= XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXS ================ #
# ============================================================================ #


# ============================================================================ #
# ======================= XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXS ================ #
# ============================================================================ #


# ============================================================================ #
# ======================= XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXS ================ #
# ============================================================================ #
