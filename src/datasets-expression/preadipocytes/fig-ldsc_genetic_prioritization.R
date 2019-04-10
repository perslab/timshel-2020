############### SYNOPSIS ###################
# Plot LDSC results


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

library(patchwork)

source(here("src/lib/load_functions.R")) # load sc-genetics library

source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/datasets-expression/preadipocytes"))

# ======================================================================= #
# =============================== LOAD DATA ================================= #
# ======================================================================= #

file.data <- "out/genetic_prioritization.1808_branch_deciles.multi_gwas.csv"
df.ldsc.raw <- read_csv(file.data)

# ======================================================================= #
# ============================ PRE-PROCESS DATA ========================= #
# ======================================================================= #

### Copy
df.ldsc <- df.ldsc.raw

### Gather
df.ldsc <- df.ldsc %>% gather(key="gwas", value="p.value", -annotation)
df.ldsc


### Process 'annotation' names: split into branch and decile
df.ldsc <- df.ldsc %>% separate(col=annotation, into=c("branch", "pseudotime_percentile"), sep="_")

### Rename to deciles
# recode: https://dplyr.tidyverse.org/reference/recode.html | When named, the argument names should be the current values to be replaced, and the argument values should be the new (replacement) values.
# mamed vector: c(name=value)
# *_top10 are the group of cells with the 10% highest pseudotime.
# *_top100 are the group of cells with the 10% lowest pseudotime.
recode_vector.deciles <- c(top100=1,
                   top90=2,
                   top80=3,
                   top70=4,
                   top60=5,
                   top50=6,
                   top40=7,
                   top30=8,
                   top20=9,
                   top10=10) 
df.ldsc <- df.ldsc %>% mutate(pseudotime_decile = recode(pseudotime_percentile, !!!recode_vector.deciles))
df.ldsc


### Rename branches
recode_vector.branches <- c(ECM="ECM",
                           oxidative="Metabolic",
                           preadipocyte="Progenitor")
df.ldsc <- df.ldsc %>% mutate(branch = recode(branch, !!!recode_vector.branches))
df.ldsc

### Filter GWAS [this vector also determines the *order of the GWAS* in the plot]
filter.gwas <- c("WHRadjBMI_UKBB_Loh2018",
                 "WHR_Pulit2019",
                 "LIPIDS_LDL_Teslovich2010",
                 "BMI_UKBB_Loh2018",
                 "HEIGHT_UKBB_Loh2018")
# "WHRadjBMI_Pulit2019" --> no reason to use this. Loh2018 is more clean

### Get all meta-analysis traits
# gwas_ids.meta_analysis <- utils.get_gwas_ids_for_meta_analysis()
# gwas_ids.meta_analysis
df.ldsc <- df.ldsc %>% filter(gwas %in% filter.gwas)
### Order GWAS by filter.gwas
df.ldsc <- df.ldsc %>% mutate(gwas=factor(gwas, levels=filter.gwas)) # Order GWAS
df.ldsc$gwas %>% levels()

### Add GWAS formated
### IMPORANTLY recode_factor() allows us to rename the factor values but KEEP the ORDER defined by filter.gwas
tmp_gwas_vector <- filter.gwas # convenience selector. If you want a specific order of the result, this vector should contain (unique) ordered values
newnames <- utils.rename_gwas(tmp_gwas_vector, style="fullname_author_year") # "fullname_author_year","fullname","abrv_author_year","abrv_year","abrv" 
rename_vector <- newnames; names(rename_vector) <- tmp_gwas_vector
df.ldsc <- df.ldsc %>% mutate(gwas_fmt = recode_factor(gwas, !!!rename_vector)) # Use a named character vector to recode factors with unquote splicing. | REF: https://dplyr.tidyverse.org/reference/recode.html
df.ldsc$gwas_fmt %>% levels()

# ======================================================================= #
# =========================== EXPORT SELECTED TRAITS ==================== #
# ======================================================================= #

df.ldsc.export <- df.ldsc %>% select(-gwas) %>% spread(key=gwas_fmt, value=p.value)
df.ldsc.export
### Write out
file.out <- "out/table-genetic_prioritization_branch_pseudotime.csv"
df.ldsc.export %>% write_csv(file.out)

# ======================================================================= #
# =============================== PLOT FUNCTION ========================== #
# ======================================================================= #

plot_genetic_prioritization <- function(df.plot) {
  p <- ggplot(df.plot, aes(x=as.factor(pseudotime_decile), y=p.value_mlog10, fill=branch)) + 
    #geom_col(position=position_dodge(), width=0.7) + 
    geom_col(position = position_dodge2(width = 0.9, preserve = "single")) + # equal widths of bars | REF: https://stackoverflow.com/questions/38101512/the-same-width-of-the-bars-in-geom-barposition-dodge
    geom_hline(yintercept = -log10(0.05/20), linetype="dashed", color="gray") + 
    labs(x="Pseudotime", y=expression(-log[10](P)), fill="Branch")
  p
  ### Colors
  colormap.branches <- c(
    Progenitor="#ecdd83",
    Metabolic="#e27268",
    ECM="#93c8bc")
  p <- p + scale_fill_manual(values=colormap.branches)
  ### Facet spacing of baseline
  p <- p + facet_grid(gwas_fmt~block_branch, # rows~columns. If using rows=vars(x) use loose the level ordering of x
                      space="free_x", scales="free_x", drop=T) + 
    theme(strip.background=element_blank(), # remove facet strip bg REF: https://stackoverflow.com/questions/14185754/remove-strip-background-keep-panel-border
          strip.text.y = element_text(size=rel(0.2)))
  ### Adjust theme
  # p <- p + theme_minimal()
  p <- p + theme(axis.text.x = element_blank())
  return(p)
}

# ======================================================================= #
# ================================== PLOT ============================= #
# ======================================================================= #

### Copy
df.plot <- df.ldsc
### Add columns
df.plot <- df.plot %>% mutate(p.value_mlog10 = -log10(p.value))
df.plot <- df.plot %>% mutate(block_branch = if_else(branch == "Progenitor", "Beginning", "End"))

### PLOT SINGLE
trait <- "WHRadjBMI_UKBB_Loh2018"
p.main <- plot_genetic_prioritization(df.plot %>% filter(gwas==trait))
p.main

### PLOT MULTI
trait_exclude <- c("WHRadjBMI_UKBB_Loh2018")
p.multi <- plot_genetic_prioritization(df.plot %>% filter(!gwas %in% trait_exclude))
p.multi

# ======================================================================= #
# ============================== PATCHWORK ================================= #
# ======================================================================= #

### Combine plots
p.patch <- p.main + p.multi + plot_layout(ncol = 1, heights = c(2, 4))
p.patch
### export
file.out <- "out/fig-genetic_prioritization_branch_pseudotime.pdf"
ggsave(plot=p.patch, filename=file.out, width=8, height=8)

# ======================================================================= #
# ============================== PLOT OLD ================================= #
# ======================================================================= #
# 
# trait <- "WHRadjBMI_UKBB_Loh2018"
# trait <- "WHR_Pulit2019"
# trait <- "BMI_UKBB_Loh2018"
# trait <- "LIPIDS_HDL_Teslovich2010"
# df.ldsc_cts.export.mod$LIPIDS_HDL_Teslovich2010
# 
# df.plot <- df.ldsc_cts %>% filter(gwas==trait)
# ggplot(df.plot, aes(x=fct_reorder(annotation,-log10(p.value)), y=-log10(p.value))) + 
#   geom_col() + 
#   geom_hline(yintercept = -log10(0.05/nrow(df.plot))) + 
#   labs(x="annotation", title=trait) +
#   coord_flip()
# # theme(axis.title.x = element_text(angle=70))


