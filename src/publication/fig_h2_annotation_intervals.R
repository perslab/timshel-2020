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

source(here("src/lib/load_functions.R")) # load sc-genetics library

setwd(here("src/ldsc"))

# ======================================================================= #
# ================================ LOAD DATA ================================ #
# ======================================================================= #

### Export (selected columns)
file.data <- here("results/h2_annotation_intervals-multi_gwas.csv.gz")
df.ldsc <- read_csv(file.data)

# ======================================================================= #
# ========== COMPUTE AVERAGE across all meta-analysis traits ============ #
# ======================================================================= #

gwas_ids.meta_analysis <- utils.get_gwas_ids_for_meta_analysis()
### Dplyr summarise

# ======================================================================= #
#============================= QUANTILE H2 PLOT ========================== #
# ======================================================================= #

library(RColorBrewer)
display.brewer.all()
display.brewer.pal(n=6, name="YlOrRd")

# x <- df.ldsc %>% distinct(annotation)

### SELECTED ANNOTATIONS
filter.annotations <- c("TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5","MEGLU10","TEGLU17","MEGLU11","TEGLU4","DEGLU4","TEINH12")
# filter.annotations <- c("TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5")
# filter.annotations <- c(filter.annotations, "Brain_Non-Myeloid.neuron", "Brain_Non-Myeloid.oligodendrocyte_precursor_cell")

### GWAS
# filter.gwas <- "BMI_UKBB_Loh2018"
# filter.gwas <- "HEIGHT_Yengo2018" # T2D_DIAMANTE_Mahajan2018
filter.gwas <- c("BMI_UKBB_Loh2018", "HEIGHT_Yengo2018") # SCZ_Pardinas2018

### Extract data
df.plot_h2q <- df.ldsc %>% filter(gwas %in% filter.gwas,
                                  annotation %in% filter.annotations)
df.plot_h2q <- df.plot_h2q %>% select(run_name, gwas, annotation, q, prop_h2g, enr, enr_se, enr_pval) # %>% gather(key="key", value="value", annotation)
df.plot_h2q

### Rename
### Rename GWAS [do this after filtering, so ensure uniqueness]
### IMPORANTLY recode_factor() allows us to rename the factor values but KEEP the ORDER defined by filter.gwas
newnames <- utils.rename_gwas(filter.gwas, style="fullname_author_year") # "fullname_author_year","fullname","abrv_author_year","abrv_year","abrv" 
rename_vector <- newnames; names(rename_vector) <- filter.gwas
df.plot_h2q <- df.plot_h2q %>% mutate(gwas_fmt = recode_factor(gwas, !!!rename_vector)) 

### PLOT
pd <- position_dodge(0.9)
p <- ggplot(df.plot_h2q, aes(x=annotation, y=enr, fill=q)) + 
# p <- ggplot(df.plot_h2q, aes(x=annotation, y=-log10(enr_pval), fill=q)) + 
  geom_col(position=pd) + 
  geom_errorbar(aes(ymin=enr-enr_se, ymax=enr+enr_se), position=pd, width = 0.01, colour="black") +
  geom_hline(yintercept = 1, linetype="dashed", color="gray") +
  scale_fill_brewer(palette = "YlOrRd",
                    breaks = levels(as.factor(df.plot_h2q$q)), 
                    labels = c(expression(I0["0"]),
                               expression(I1["(0-0.2]"]),
                               expression(I2["(0.2-0.4]"]),
                               expression(I3["(0.4-0.6]"]),
                               expression(I4["(0.6-0.8]"]),
                               expression(I5["(0.8-1]"]))) +
  theme(legend.text.align = 0) + # needed for aligning legend text to the left
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="", y=expression("h"[2]*" enrichment"), fill="ES interval") # title=expression("h"[2]*" enrichment") 
p <- p + facet_grid(gwas_fmt~run_name, space="free_x", scales="free_x", drop=T)
# p <- p + facet_wrap(~gwas_fmt, ncol=1)
p <- p + theme(panel.grid.major = element_blank(), # REF: https://stackoverflow.com/questions/14185754/remove-strip-background-keep-panel-border
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black"))
p
file.out <- sprintf("out.ldsc_h2_intervals.mb_fdr_celltypes.%s.pdf", paste(filter.gwas, collapse="-"))
# ggsave(filename=file.out, p, width=12, height=6)



# ======================================================================= #
# =========================== EXPORT to results ========================= #
# ======================================================================= #

### Export (selected columns)
#file.out <- here("results/h2_annotations_quantile-multi_gwas.csv.gz")
# file.out
# df.ldsc %>% select(-Category) %>% arrange(gwas) %>% write_csv(file.out)


# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #



# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #
