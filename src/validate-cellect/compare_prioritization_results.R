############### SYNOPSIS ###################
# Compare CELLECT snakemake pipeline results to "190306_fix results".
# Overall aim is to show that CELLECT produces the same results as PNT 'manual python workflow' WHEN GIVEN THE SAME ESmu INPUT
# (Some minor deviations are allowed because of how 'all genes' are automatically handled in the snakemake workflow)


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

# source(here("src/lib/load_functions.R")) # load sc-genetics library

# ======================================================================= #
# =============================== LOAD DATA ================================= #
# ======================================================================= #

#### thesis (before fix) | NOT IMPORTANT HERE
file.results <- here("src/validate-cellect/data/19_thesis_results_that_is_before_es_fix.prioritization_celltypes--mousebrain_all.multi_gwas.csv.gz")
df.before_fix <- read_csv(file.results) %>% filter(gwas=="BMI_UKBB_Loh2018")
df.before_fix <- df.before_fix %>% select(annotation, p.value)

### 190306_es_fix
file.results <- here("src/validate-cellect/data/celltypes.mousebrain_190306_es_fix.all__BMI_UKBB_Loh2018.cell_type_results.txt")
df.fix <- read_tsv(file.results)
df.fix <- df.fix %>% separate(col=Name, into=c("X1", "X2", "X3", "annotation", "X4"), sep="\\.")
df.fix <- df.fix %>% select(annotation, p.value="Coefficient_P_value")
df.fix

### snakemake run from 191110 with 190306_es_fix ESmu matrix as input.
# CP FROM ---> /scratch/timshel/BMIcelltypes/190306_mb_replication_tmp/
file.results <- here("src/validate-cellect/data/191110-cellect_snakemake_with_190306_es_fix_as_input.prioritization.csv")
df.fix_snake <- read_csv(file.results)
df.fix_snake <- df.fix_snake %>% filter(gwas=="BMI_UKBB_Loh2018")
df.fix_snake <- df.fix_snake %>% select(annotation, p.value=pvalue)

# ======================================================================= #
# =============================== COMPARE ================================= #
# ======================================================================= #


### snakemake vs manual
### THIS IS THE COMPARISON OF INTEREST
df.cmp <- full_join(df.fix_snake %>% rename(snakemake=p.value), 
                    df.fix %>% rename(PT_workflow=p.value), 
                    by="annotation")
ggplot(df.cmp, aes(snakemake, PT_workflow)) + geom_point() + geom_abline()
# ---> perfect correlation!


### thesis vs es_fix - 'NOT IMPORTANT'
### NB: we have already investigated this in GDOCs 'RESULTS _ 190316 CELL-TYPE PRIORITIZATION _ Fixed ES calculation'
df.cmp <- full_join(df.before_fix %>% rename(before_fix=p.value), df.fix %>% rename(es_fix=p.value), by="annotation")
df.cmp
ggplot(df.cmp, aes(before_fix, es_fix)) + geom_point() + geom_abline() + scale_x_log10() + scale_y_log10()


# ======================================================================= #
# =============================== XXXXX ================================= #
# ======================================================================= #

