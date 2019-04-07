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

library(RColorBrewer)
# display.brewer.all()
# display.brewer.pal(n=6, name="YlOrRd")

source(here("src/lib/load_functions.R")) # load sc-genetics library

setwd(here("src/publication"))


# ======================================================================= #
# ======================== FUNCTION: BASE PLOT ===================== #
# ======================================================================= #

plot_h2_annotation_intervals <- function(df.plot_h2q) {
  pd <- position_dodge(0.9)
  p <- ggplot(df.plot_h2q, aes(x=annotation, y=enr, fill=q)) + 
    geom_col(position=pd) + 
    geom_errorbar(aes(ymin=enr-enr_se, ymax=enr+enr_se), position=pd, width = 0.01, colour="black") +
    geom_hline(yintercept = 1, linetype="dashed", color="gray") +
    scale_fill_brewer(palette = "YlOrRd",
                      breaks = levels(as.factor(df.plot_h2q$q)), # IMPORTANT: q is q0,...q5
                      labels = c(expression(I0["0"]),
                                 expression(I1["(0-0.2]"]),
                                 expression(I2["(0.2-0.4]"]),
                                 expression(I3["(0.4-0.6]"]),
                                 expression(I4["(0.6-0.8]"]),
                                 expression(I5["(0.8-1]"]))) +
    theme(legend.text.align = 0) + # needed for aligning legend text to the left
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    labs(x="", y=expression("h"[2]*" enrichment"), fill="ES interval") # title=expression("h"[2]*" enrichment") 
  return(p)
}

# ======================================================================= #
# ================================ LOAD DATA ================================ #
# ======================================================================= #

### Export (selected columns)
file.data <- here("results/h2_annotation_intervals-multi_gwas.csv.gz")
df.ldsc <- read_csv(file.data)

### Rename GWAS
newnames <- utils.rename_gwas(df.ldsc$gwas, style="fullname_author_year") # "fullname_author_year","fullname","abrv_author_year","abrv_year","abrv" 
rename_vector <- newnames; names(rename_vector) <- df.ldsc$gwas
df.ldsc <- df.ldsc %>% mutate(gwas_fmt = recode_factor(gwas, !!!rename_vector)) 

### Rename run_name
df.ldsc <- df.ldsc %>% mutate(run_name = case_when(
  run_name == "celltypes.mousebrain" ~ "Mousebrain",
  run_name == "celltypes.tabula_muris" ~ "Tabula Muris",
  TRUE ~ as.character(run_name))
)

# ======================================================================= #
# ================= INTERVALS H2 PLOT | MAIN [MB only] ================== #
# ======================================================================= #

### SELECTED ANNOTATIONS
filter.annotations <- c("TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5","MEGLU10","TEGLU17","MEGLU11","TEGLU4","DEGLU4","TEINH12")
### SELECTED GWAS
filter.gwas <- "BMI_UKBB_Loh2018"
### Extract data
df.plot_h2q <- df.ldsc %>% filter(gwas %in% filter.gwas, annotation %in% filter.annotations)

### PLOT
p <- plot_h2_annotation_intervals(df.plot_h2q)
p
file.out <- sprintf("figs/fig_h2_annotation_intervals.main.mb_fdr_celltypes.%s.pdf", paste(filter.gwas, collapse="-"))
ggsave(filename=file.out, p, width=12, height=6)



# ======================================================================= #
# ================== COMPUTATION FOR META-ANALYSIS ====================== #
# ======================================================================= #

### COMPUTE AVERAGE across all meta-analysis traits
gwas_ids.meta_analysis <- utils.get_gwas_ids_for_meta_analysis()
df.ldsc.meta_analysis <- df.ldsc %>% 
  filter(gwas %in% gwas_ids.meta_analysis) %>% # do computations only on meta-analysis gwas
  group_by(annotation, q) %>% 
  summarise(
    enr_se = sd(enr, na.rm=T), # OBS: we take the sd of enr. Calculate this before you 'update' enr
    enr = mean(enr, na.rm=T),
    run_name = unique(run_name),
    gwas = "META_ANALYSIS",
    gwas_fmt = sprintf("Meta-analysis (n=%s)", length(gwas_ids.meta_analysis))
  ) %>% 
  ungroup()

### Row concatenate
df.ldsc.meta_analysis <- bind_rows(df.ldsc.meta_analysis, df.ldsc)

# ======================================================================= #
# ================ INTERVALS H2 PLOT | SOM [BMI+HEIGHT/META] ================= #
# ======================================================================= #

### GWAS [*SWITCH*]
filter.gwas <- c("BMI_UKBB_Loh2018", "META_ANALYSIS") 
# filter.gwas <- c("BMI_UKBB_Loh2018", "HEIGHT_Yengo2018") # SCZ_Pardinas2018


### SELECTED ANNOTATIONS
filter.annotations <- c("TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5","MEGLU10","TEGLU17","MEGLU11","TEGLU4","DEGLU4","TEINH12")
filter.annotations <- c(filter.annotations, "Brain_Non-Myeloid.neuron", "Brain_Non-Myeloid.oligodendrocyte_precursor_cell")

### Extract data [OBS: df.ldsc.meta_analysis]
df.plot_h2q <- df.ldsc.meta_analysis %>% filter(gwas %in% filter.gwas, annotation %in% filter.annotations)
df.plot_h2q <- df.plot_h2q %>% mutate(annotation = utils.rename_annotations.tabula_muris(annotation, style="tissue - celltype"))
df.plot_h2q

### PLOT
p <- plot_h2_annotation_intervals(df.plot_h2q)
p <- p + facet_grid(gwas_fmt~run_name, space="free_x", scales="free_x", drop=T)
p <- p + theme(panel.grid.major = element_blank(), # REF: https://stackoverflow.com/questions/14185754/remove-strip-background-keep-panel-border
               panel.grid.minor = element_blank(),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"))
p
file.out <- sprintf("figs/fig_h2_annotation_intervals.som.fdr_celltypes.%s.pdf", paste(filter.gwas, collapse="-"))
# ggsave(filename=file.out, p, width=12, height=10)



# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #



# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #
