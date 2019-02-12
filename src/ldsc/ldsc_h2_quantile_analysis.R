############### SYNOPSIS ###################
### AIM: Perform h2 quantile analysis (h2 enrichment for each quantile of a continuous annotation)

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #

library(tidyverse)
library(here)

library(GGally)
library(corrr) # devtools::install_github("drsimonj/corrr")
library(gghighlight)


dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib/"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library


setwd(here("src/ldsc"))

# ======================================================================= #
# ================================ LOAD DATA ================================ #
# ======================================================================= #

# This output file has one row for each quantile (starting with lowest values) 
# column summarizing the heritability explained by each quantile, its enrichment and corresponding standard error and P value.

### Load - MULTI GWAS AND ANNOTATIONS
dir.data <- "/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc_h2/"
filenames <- list.files(path=dir.data,  pattern="(.*).results_quantile_qfixed_h2")
# filenames <- filenames[!grepl(pattern="__CONDITIONAL__", filenames, perl=T)] # exclude any conditional results
list.dfs <- lapply(file.path(dir.data, filenames), read_tsv)
list.dfs <- lapply(list.dfs, function(df) {df %>% mutate(q=paste0("q", seq(0,n()-1)))}) # add quantile
names(list.dfs) <- stringr::str_match(filenames, pattern="(.*).results_quantile_qfixed_h2")[,2] 
names(list.dfs)
df.ldsc <- list.dfs %>% bind_rows(.id="run_str_full")

### Split sting
df.ldsc <- df.ldsc %>% separate(col=run_str_full, into=c("run_name", "annotation", "gwas"), sep="__") # e.g. run_str = "wgcna.mousebrain-190111__wheat3__T2D_UKBB_DIAMANTE_Mahajan2018"
### Filter
# df.ldsc <- df.ldsc %>% filter(gwas %in% c("BMI_UKBB_Loh2018", "T2D_DIAMANTE_Mahajan2018"))

# ======================================================================= #
#============================= QUANTILE H2 PLOT ========================== #
# ======================================================================= #

library(RColorBrewer)
display.brewer.all()
display.brewer.pal(n=6, name="YlOrRd")


### BMI_UKBB_Loh2018 FDR significant annotations
filter.annotations <- c("TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5","MEGLU10","TEGLU17","MEGLU11","TEGLU4","DEGLU4","TEINH12")
# filter.annotations <- c("TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5")
df.plot_h2q <- df.ldsc %>% filter(gwas %in% c("BMI_UKBB_Loh2018"),
                              annotation %in% filter.annotations)
df.plot_h2q <- df.plot_h2q %>% select(annotation, q, prop_h2g, enr, enr_pval) # %>% gather(key="key", value="value", annotation)
p <- ggplot(df.plot_h2q, aes(x=annotation, y=enr, fill=q)) + 
# p <- ggplot(df.plot_h2q, aes(x=annotation, y=-log10(enr_pval), fill=q)) + 
  geom_col(position=position_dodge()) + 
  scale_fill_brewer(palette = "YlOrRd") + 
  labs(title=expression("Quantile "*"h"[2]*" enrichment"), y=expression("h"[2]*" enrichment"))
p
# ggsave(filename="out.ldsc_h2_quantile.mb_celltypes.BMI_UKBB_Loh2018.pdf", p, width=12, height=6)

### HEIGHT_Yengo2018 or T2D_DIAMANTE_Mahajan2018 --> negative control
filter.annotations <- c("TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5","MEGLU10","TEGLU17","MEGLU11","TEGLU4","DEGLU4","TEINH12")
df.plot_h2q <- df.ldsc %>% filter(gwas %in% c("HEIGHT_Yengo2018"),
                                  annotation %in% filter.annotations)
df.plot_h2q <- df.plot_h2q %>% select(annotation, q, prop_h2g, enr, enr_pval) # %>% gather(key="key", value="value", annotation)
p <- ggplot(df.plot_h2q, aes(x=annotation, y=enr, fill=q)) + 
  geom_col(position=position_dodge()) + 
  scale_fill_brewer(palette = "YlOrRd") + 
  labs(title=expression("Quantile "*"h"[2]*" enrichment"), y=expression("h"[2]*" enrichment"))
p
# ggsave(filename="out.ldsc_h2_quantile.mb_celltypes.HEIGHT_Yengo2018.pdf", p, width=12, height=6)

# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #



# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #
