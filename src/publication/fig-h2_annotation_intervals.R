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
    labs(x="", y="Heritability enrichment", fill="ES interval") # title=expression("h"[2]*" enrichment") 
  return(p)
}

plot_h2_annotation_intervals.lollipop <- function(df.plot_h2q) {
  pd <- position_dodge(0.9)
  p <- ggplot(df.plot_h2q, aes(x=annotation, y=enr, group=q)) + 
    geom_hline(yintercept = 1, linetype="dashed", color="gray") +
    geom_linerange(aes(x=annotation, ymin=0, ymax=enr), position=pd, color="grey",  alpha=0.5) + # REF: https://stackoverflow.com/a/21922792/6639640
    geom_point(aes(color=q), size=3, position=pd) + 
    geom_errorbar(aes(ymin=enr-enr_se, ymax=enr+enr_se), position=pd, width = 0.01, colour="black") +
    scale_color_brewer(palette = "YlOrRd",
                      breaks = levels(as.factor(df.plot_h2q$q)), # IMPORTANT: q is q0,...q5
                      labels = c(expression(I0["0"]),
                                 expression(I1["(0-0.2]"]),
                                 expression(I2["(0.2-0.4]"]),
                                 expression(I3["(0.4-0.6]"]),
                                 expression(I4["(0.6-0.8]"]),
                                 expression(I5["(0.8-1]"]))) +
    theme(legend.text.align = 0) + # needed for aligning legend text to the left
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    labs(x="", y="Heritability enrichment", color="ES interval") # title=expression("h"[2]*" enrichment") 
    # WORKS --> labs(x="", y=expression(atop(h^{2}~enrichment, (proportion~h^{2}/annotation~size))), color="ES interval") # title=expression("h"[2]*" enrichment") 
  p
  return(p)
}


# ======================================================================= #
# ================================ LOAD DATA ================================ #
# ======================================================================= #

### Export (selected columns)
file.data <- here("results/h2_annotation_intervals.multi_gwas.csv.gz")
df.ldsc <- read_csv(file.data)

### Rename GWAS
tmp_gwas_vector <- df.ldsc$gwas # convenience selector. If you want a specific order of the result, this vector should contain (unique) ordered values
newnames <- utils.rename_gwas(tmp_gwas_vector, style="fullname_author_year") # "fullname_author_year","fullname","abrv_author_year","abrv_year","abrv" 
rename_vector <- newnames; names(rename_vector) <- tmp_gwas_vector
df.ldsc <- df.ldsc %>% mutate(gwas_fmt = recode_factor(gwas, !!!rename_vector)) # Use a named character vector to recode factors with unquote splicing. | REF: https://dplyr.tidyverse.org/reference/recode.html
# OLD OK DELETE
# newnames <- utils.rename_gwas(df.ldsc$gwas, style="fullname_author_year") # "fullname_author_year","fullname","abrv_author_year","abrv_year","abrv" 
# rename_vector <- newnames; names(rename_vector) <- df.ldsc$gwas
# df.ldsc <- df.ldsc %>% mutate(gwas_fmt = recode_factor(gwas, !!!rename_vector)) 

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
# filter.annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")
# Prioritized cell-types ORDERED BY Prop-h2 from fig_h2_annotations
filter.annotations <- c("DEGLU4","MEGLU10","DEGLU5","MEGLU11","TEINH12","MEGLU1","DEINH3","MEINH2","TEGLU23","TEGLU17","TEGLU4") 
### SELECTED GWAS
filter.gwas <- "BMI_UKBB_Loh2018"
### Extract data
df.plot_h2q <- df.ldsc %>% 
  filter(gwas %in% filter.gwas, annotation %in% filter.annotations) %>%
  mutate(annotation = factor(annotation, levels=filter.annotations)) # set order

### PLOT
p <- plot_h2_annotation_intervals.lollipop(df.plot_h2q)
p <- p + coord_flip()
p <- p + theme(axis.text.x = element_text(angle=0, hjust=0.5))
p
file.out <- sprintf("figs/fig_h2_annotation_intervals.main.mb_fdr_celltypes.%s.pdf", paste(filter.gwas, collapse="-"))
ggsave(filename=file.out, p, width=6, height=4)



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
# filter.gwas <- c("BMI_UKBB_Loh2018", "META_ANALYSIS") 
filter.gwas <- c("BMI_UKBB_Loh2018", "HEIGHT_UKBB_Loh2018", "WHRadjBMI_UKBB_Loh2018") # "RA_Okada2014", "MS_Patsopoulos2011", "SCZ_Pardinas2018"

### SELECTED ANNOTATIONS
filter.annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")
filter.annotations <- c(filter.annotations, "Brain_Non-Myeloid.neuron", "Brain_Non-Myeloid.oligodendrocyte_precursor_cell")

### Extract data [OBS: df.ldsc.meta_analysis]
df.plot_h2q <- df.ldsc.meta_analysis %>% filter(gwas %in% filter.gwas, annotation %in% filter.annotations)
df.plot_h2q <- df.plot_h2q %>% mutate(annotation = utils.rename_annotations.tabula_muris(annotation, style="tissue - celltype"))
df.plot_h2q

### PLOT
p <- plot_h2_annotation_intervals.lollipop(df.plot_h2q)
p <- p + facet_grid(gwas_fmt~run_name, 
                    space="free_x", scales="free_x", drop=T)
                    # switch="x") # If "x", the top labels will be displayed to the bottom. If "y", the right-hand side labels will be displayed to the left. Can also be set to "both".
p <- p + theme(panel.grid.major = element_blank(), # REF: https://stackoverflow.com/questions/14185754/remove-strip-background-keep-panel-border
               panel.grid.minor = element_blank(),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               strip.text.y = element_text(size=rel(0.5)) # gwas text size smaller for easier post edditing
               )
p
file.out <- sprintf("figs/fig_h2_annotation_intervals.som.fdr_celltypes.%s.pdf", paste(filter.gwas, collapse="-"))
ggsave(filename=file.out, p, width=12, height=8)



# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #



# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #
