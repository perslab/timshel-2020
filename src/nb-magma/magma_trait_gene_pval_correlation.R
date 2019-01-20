############### SYNOPSIS ###################
### Investigate correlation between MAGMA P-value for different traits


library(tidyverse)
library(ggpubr) # http://www.sthda.com/english/rpkgs/ggpubr/index.html


# ========================================================================== #
# ======================= SETUP ======================= #
# ========================================================================== #

wd <- "/projects/timshel/sc-genetics/sc-genetics/src/magma/"
setwd(wd)

# ========================================================================== #
# ======================= READ MAGMA GWAS PVALUES FILES ======================= #
# ========================================================================== #


### READ MAGMA PVALUES
DIR.gwas <- "/scratch/tmp-magma_gwas/"
filenames <- list.files(DIR.gwas, pattern="*genes.out")
list.magma <- lapply(file.path(DIR.gwas, filenames), read_table) # read_table() and read_table2() are designed to read the type of textual data where each column is separated by one (or more) columns of space.
names(list.magma) <- stringr::str_match(filenames, pattern="(.*).txt*")[,2] # e.g. "AN_PGC_Duncan2017.txt.10UP.1.5DOWN.genes.out" --> "AN_PGC_Duncan2017"
names(list.magma)
df.magma.raw <- bind_rows(list.magma, .id="gwas")

### subset
gwas_of_interest <- c("BMI_Yengo2018",
  "EA2_Okbay2016",
  "EA3_Lee2018",
  "ASD_iPSYCH_PGC_Grove2018",
  "SCZ_Ripke2014",
  "HEIGHT_Yengo2018",
  "blood_EOSINOPHIL_COUNT",
  "1KG_phase3_EUR_null_gwas_P1")
df.magma <- df.magma.raw %>% filter(gwas %in% gwas_of_interest)

# ========================================================================== #
# ======================= PLOT ======================= #
# ========================================================================== #

### reformat df for plotting (a small trick for getting "BMI" ZSTATs as the column to plot)
df.plot <- df.magma %>% 
  select(gwas, GENE, ZSTAT) %>% 
  spread(key=gwas, value=ZSTAT) %>%
  gather(key="gwas", value=ZSTAT, -BMI_Yengo2018, -GENE)

### PLOT | single-gwas
ggplot(df.plot %>% filter(gwas=="SCZ_Ripke2014"), aes(x=ZSTAT, y=BMI_Yengo2018)) + geom_point(alpha=0.3) + geom_smooth(method = "lm", se = FALSE) + ggpubr::stat_cor(method = "spearman")

### PLOT | multi-gwas
ggplot(df.plot, aes(x=ZSTAT, y=BMI_Yengo2018)) + geom_point(alpha=0.3) + geom_smooth(method = "lm", se = FALSE) + ggpubr::stat_cor(method = "spearman") + facet_wrap(~gwas,ncol=1)



