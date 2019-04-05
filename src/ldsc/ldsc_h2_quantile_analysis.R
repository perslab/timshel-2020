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
#============================= TM Label formatter ========================== #
# ======================================================================= #
# # Self-defined label formatting function
# tm_annotation_formatter <- function(annotation_tissue_celltype) {
#   # INPUT: annotation_tissue_celltype, string e.g. 'Marrow.Slamf1-negative_multipotent_progenitor_cell'
#   # USAGE 1: df %>% mutate
#   # USAGE 2: ggplot() + scale_y_discrete(label=tm_annotation_formatter)
#   label <- stringr::str_split_fixed(annotation_tissue_celltype, pattern="\\.", n=Inf)[,2]
#   label <- stringr::str_replace_all(label, pattern="-", " - ") # add space to any hyphens
#   label <- stringr::str_replace_all(label, pattern="_", " ") # convert _ to space
#   # label <- stringr::str_to_sentence(label) # title case
#   substr(label, 1, 1) <- toupper(substr(label, 1, 1)) # capitalize first character
#   return(label)
# }

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
# filter.annotations <- c(filter.annotations, "Brain_Non-Myeloid.neuron", "Brain_Non-Myeloid.oligodendrocyte")



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
