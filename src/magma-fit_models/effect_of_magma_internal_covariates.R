### SYNOPSIS: test effect of MAGMA 'correction' (internal covariates) on the significance of the cell-type prioritization.
### OBS: this investigation draws its conclusions based on MAGMA analysis of the 'mousebrain.sem.specificity_quantiles' SEM data.
### It is possible that other SEM data would give different results.

library(tidyverse)

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/"
setwd(wd)

# ======================================================================= #
# =============================== MAGMA COMMANDS ================================ #
# ======================================================================= #


# ### MAGMA correct_genesize_density
# magma \
# --gene-results /scratch/tmp-magma_gwas/BMI_Yengo2018.txt.10UP.1.5DOWN.genes.raw \
# --gene-covar /raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/mousebrain.sem.specificity_quantiles.entrez.txt \
# --model \
# correct=include,genesize,density \
# direction-covar=greater \
# --settings abbreviate=0 gene-info \
# --out /raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/magma.specificity_quantiles.correct_genesize_density

# ### MAGMA correct_density
# magma \
# --gene-results /scratch/tmp-magma_gwas/BMI_Yengo2018.txt.10UP.1.5DOWN.genes.raw \
# --gene-covar /raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/mousebrain.sem.specificity_quantiles.entrez.txt \
# --model \
# correct=include,density \
# direction-covar=greater \
# --settings abbreviate=0 gene-info \
# --out /raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/magma.specificity_quantiles.correct_density


# ### MAGMA correct_genesize
# magma \
# --gene-results /scratch/tmp-magma_gwas/BMI_Yengo2018.txt.10UP.1.5DOWN.genes.raw \
# --gene-covar /raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/mousebrain.sem.specificity_quantiles.entrez.txt \
# --model \
# correct=include,genesize \
# direction-covar=greater \
# --settings abbreviate=0 gene-info \
# --out /raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/magma.specificity_quantiles.correct_genesize

# ### MAGMA correct_none
# magma \
# --gene-results /scratch/tmp-magma_gwas/BMI_Yengo2018.txt.10UP.1.5DOWN.genes.raw \
# --gene-covar /raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/mousebrain.sem.specificity_quantiles.entrez.txt \
# --model \
# correct=none \
# direction-covar=greater \
# --settings abbreviate=0 gene-info \
# --out /raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/magma.specificity_quantiles.correct_none
# 
# 
# ### MAGMA correct_all
# magma \
# --gene-results /scratch/tmp-magma_gwas/BMI_Yengo2018.txt.10UP.1.5DOWN.genes.raw \
# --gene-covar /raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/mousebrain.sem.specificity_quantiles.entrez.txt \
# --model \
# correct=all \
# direction-covar=greater \
# --settings abbreviate=0 gene-info \
# --out /raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/magma.specificity_quantiles.correct_all

# ============================================================================ #
# ======================= READ MAGMA *.gsa.out (cell-type enrichment) ================ #
# ============================================================================ #

### SNIPPET *.gsa.out
# # MEAN_SAMPLE_SIZE = 68378
# # TOTAL_GENES = 13985
# # TEST_DIRECTION = one-sided, positive (set), one-sided, positive (covar)
# # CONDITIONED_INTERNAL = gene size, gene density, sample size, inverse mac, log(gene size), log(gene density), log(sample size), log(inverse mac)
# VARIABLE      TYPE  NGENES         BETA     BETA_STD           SE            P
# ENT9         COVAR   13985    0.0044786     0.060696   0.00089045   2.4932e-07
# ENT8         COVAR   13985    0.0031883     0.043565   0.00088087   0.00014816
# ENT6         COVAR   13985    0.0032336     0.044283   0.00087604   0.00011213
# ENT5         COVAR   13985    0.0029439     0.040209   0.00087223   0.00037008

### read *.gsa.out data
DIR.sem_data <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/data-effect_of_magma_internal_covariates/"
filenames <- list.files(path=DIR.sem_data,  pattern="*.gsa.out")
list.df_sem <- lapply(file.path(DIR.sem_data, filenames), read_table, comment = "#")
filenames_shorten <- stringr::str_match(filenames, "magma.specificity_quantiles\\.(.*)\\.gsa.out")[,2] # e.g. "celltype_expr.specificity_quantiles.hsapiens_orthologs.csv.gz" --> "specificity_quantiles"
names(list.df_sem) <- filenames_shorten
names(list.df_sem)


### Row bind
df <- bind_rows(list.df_sem, .id="correction")

### Process
df <- df %>% mutate(minus_log10P=-log10(P))


# ============================================================================ #
# ======================= PLOT cell-type enrichment ================ #
# ============================================================================ #

df.plot <- df %>% select(correction, cell_type=VARIABLE, minus_log10P) %>% spread(key=correction, value=minus_log10P)
df.plot

# cell_type correct_all correct_genesize_density correct_none
# <chr>          <dbl>                    <dbl>        <dbl>
#   1 ABC            0.367                    0.375       0.0301
# 2 ACBG           2.48                     2.54        2.02  
# 3 ACMB           1.17                     1.19        1.04  
# 4 ACNT1          1.23                     1.38        1.03  
# 5 ACNT2          2.46                     2.53        2.70  


ggplot(df.plot, aes(x=correct_all, y=correct_none)) + geom_point() + geom_abline() + labs(title="-log10_Pval(MAGMA cell-type priori)")
ggplot(df.plot, aes(x=correct_genesize_density, y=correct_all)) + geom_point() + geom_abline() + labs(title="-log10_Pval(MAGMA cell-type priori)")

### Effect of genesize
ggplot(df.plot, aes(x=correct_genesize_density, y=correct_density)) + geom_point() + geom_abline() + labs(title="-log10_Pval(MAGMA cell-type priori)") + theme_gray() # the default
ggplot(df.plot, aes(x=correct_density, y=correct_none)) + geom_point() + geom_abline() + labs(title="-log10_Pval(MAGMA cell-type priori)") + theme_gray() # the default

### Effect of density
ggplot(df.plot, aes(x=correct_genesize_density, y=correct_genesize)) + geom_point() + geom_abline() + labs(title="-log10_Pval(MAGMA cell-type priori)") + theme_gray() # the default
ggplot(df.plot, aes(x=correct_genesize, y=correct_none)) + geom_point() + geom_abline() + labs(title="-log10_Pval(MAGMA cell-type priori)") + theme_gray() # the default





# ============================================================================ #
# ======================= READ MAGMA *.gsa.genes.out (Gene P-values) ================ #
# ============================================================================ #


### SNIPPET *.gsa.genes.out
# # CONDITIONED_INTERNAL = gene size, gene density, sample size, inverse mac, log(gene size), log(gene density), log(sample size), log(inverse mac)
# GENE       CHR      START       STOP  NSNPS  NPARAM       N        ZSTAT  ZFITTED_BASE  ZRESID_BASE
# 375790       1     945503     992999      1       1  605446     -0.56891   -0.00016317     -0.56875
# 54991        1    1015698    1061736     24       5  649301     -0.81953   -0.00016317     -0.81937
# 254173       1    1099286    1134815     11       4  643729     -0.86689   -0.00016317     -0.86672

### read *.gsa.genes.out data
DIR.sem_data <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/data-effect_of_magma_internal_covariates/"
filenames <- list.files(path=DIR.sem_data,  pattern="*.gsa.genes.out")
list.df_sem <- lapply(file.path(DIR.sem_data, filenames), read_table, comment = "#")
filenames_shorten <- stringr::str_match(filenames, "magma.specificity_quantiles\\.(.*)\\.gsa.genes.out")[,2] # e.g. "celltype_expr.specificity_quantiles.hsapiens_orthologs.csv.gz" --> "specificity_quantiles"
names(list.df_sem) <- filenames_shorten
names(list.df_sem)


### Row bind
df <- bind_rows(list.df_sem, .id="correction")


# ============================================================================ #
# ======================= PLOT ZSTAT ================ #
# ============================================================================ #

df.plot <- df %>% select(correction, GENE, ZSTAT) %>% spread(key=correction, value=ZSTAT)
df.plot

ggplot(df.plot, aes(x=correct_all, y=correct_none)) + geom_point() + geom_abline() + labs(title="MAGMA ZSTAT")
ggplot(df.plot, aes(x=correct_genesize, y=correct_none)) + geom_point() + geom_abline() + labs(title="MAGMA ZSTAT")
ggplot(df.plot, aes(x=correct_density, y=correct_none)) + geom_point() + geom_abline() + labs(title="MAGMA ZSTAT")

