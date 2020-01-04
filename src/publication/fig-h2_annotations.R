############### SYNOPSIS ###################
# Heritability explained figure


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

source(here("src/lib/load_functions.R")) # load sc-genetics library

setwd(here("src/publication"))

# ======================================================================= #
# ============================ READ DATA =============================== #
# ======================================================================= #

file.data <- here("results/cellect_ldsc/heritability.csv")
df <- read_csv(file.data)

### Rename GWAS [do this after filtering, so ensure uniqueness]
tmp_gwas_vector <- df$gwas # convenience selector. If you want a specific order of the result, this vector should contain (unique) ordered values
newnames <- utils.rename_gwas(tmp_gwas_vector, style="abrv_author_year") # "fullname_author_year","fullname","abrv_author_year","abrv_year","abrv" 
rename_vector <- newnames; names(rename_vector) <- tmp_gwas_vector
df <- df %>% mutate(gwas_fmt = recode_factor(gwas, !!!rename_vector)) # Use a named character vector to recode factors with unquote splicing. | REF: https://dplyr.tidyverse.org/reference/recode.html

### SELECTED ANNOTATIONS
filter.annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")
# filter.annotations <- c(filter.annotations, "Brain_Non-Myeloid.neuron", "Brain_Non-Myeloid.oligodendrocyte_precursor_cell")

### Filter
df <- df %>% filter(
  (annotation=="Limb_Muscle.mesenchymal_stem_cell" & gwas=="HEIGHT_UKBB_Loh2018") |
  (annotation=="Liver.hepatocyte" & gwas=="LIPIDS_LDL_Teslovich2010") |
  (annotation=="Spleen.T_cell" & gwas=="RA_Okada2014") |
  (annotation=="Pancreas.type_B_pancreatic_cell" & gwas=="T2D_UKBB_Loh2018") |
  (annotation %in% filter.annotations & gwas=="BMI_UKBB_Loh2018")
  )

### Rename TM annotations
df <- df %>% mutate(annotation = utils.rename_annotations.tabula_muris(annotation, style="tissue - celltype"))


### Add 'block' colunm separating BMI and other GWAS
df <- df %>% mutate(gwas_block = factor(if_else(gwas=="BMI_UKBB_Loh2018", "A_block", "B_block"))) # ordered factor.
### Order entries
# df <- df %>% mutate(idx_sort=order(gwas_block, Prop._h2)) # ---> THIS DOES NOT SORT CORRECTLY AND I HAVE NO IDEA WHY.
# ---> df.x <- df %>% mutate(idx_sort=base::order(gwas_block)) %>% arrange(idx_sort) # | even this does not work
# ---> df.x <- df[with(df, order(gwas_block, Prop._h2)),] # | this sorts correctly
df <- df %>% arrange(gwas_block, Prop._h2) # We want "other traits" to display in the top of the figure.
df <- df %>% mutate(annotation=factor(annotation, levels=annotation)) 


### Add percentage
df <- df %>% mutate(Pct._h2=Prop._h2*100,
                    Pct._h2_std_error=Prop._h2_std_error*100)


# ======================================================================= #
# ============================ PLOT =============================== #
# ======================================================================= #



### y=pct_h2; color=prop_snps
p <- ggplot(df, aes(x=annotation, y=Pct._h2)) + 
  geom_segment(aes(x=annotation, xend=annotation, y=0, yend=Pct._h2), color="grey") +
  geom_point(aes(color=Prop._SNPs), size=10) + 
  geom_errorbar(aes(ymin=Pct._h2-Pct._h2_std_error, ymax=Pct._h2+Pct._h2_std_error), width = 0.01, colour="black") +
  geom_text(aes(label=paste0(round(Pct._h2, 1), "%")), size=2.5, color="white", fontface="bold") +
  geom_text(data=df %>% filter(gwas_block=="B_block"), aes(y=35, label=gwas_fmt), size=4, color="black", hjust=0) +
  coord_flip(
    clip = 'off' # allow GWAS text outside plot [if y=25]
  ) + 
  #scale_color_continuous(low='darkmagenta', high='darkred') +
  scale_color_viridis(option="magma", begin=0.3, end=0.8, direction=-1) +
  labs(x="", y="Proportion of heritability (%)", color="Annotation size") + # y=expression("%"~h^{2})
  theme_classic() + 
  theme(plot.margin =  unit(c(1,4,1,1), "cm")) # trbl
p <- p + guides(color=guide_colorbar(direction = "horizontal")) + # make colorbar horizontal
  theme(legend.position="top")
p
p <- p + facet_grid(rows=vars(gwas_block), 
             space="free_y", scales="free_y", drop=T)
p <- p + theme(panel.grid.major = element_blank(), # remove gridlines REF: https://stackoverflow.com/questions/14185754/remove-strip-background-keep-panel-border
               panel.grid.minor = element_blank(),
               strip.background = element_blank(), # remove facet strip bg
               strip.text.y = element_blank()
)
p


file.out <- "figs/fig_h2_annotations.pdf"
ggsave(file.out, plot=p, width=8, height=6)







# ======================================================================= #
# ========================== FIRST VERSIONS [KEEP] ======================== #
# ======================================================================= #


### y=pct_h2
p <- ggplot(df, aes(x=annotation, y=Pct._h2, label=paste0(round(Pct._h2, 1), "%"))) + 
  geom_segment(aes(x=annotation, xend=annotation, y=0, yend=Pct._h2), color="grey") +
  geom_point(aes(color=gwas), size=10) + 
  geom_text(size=3) +
  coord_flip() + 
  labs(x="", y=expression("%"~h^{2}~explained)) +
  guides(color=F)
p
file.out <- "figs/fig_h2_annotations.simple_with_labels.pdf"
# ggsave(file.out, width=10, height=5)



### y=pct_h2 + size=enrichment
ggplot(df, aes(x=annotation, y=Pct._h2)) + 
  geom_point(aes(size=Enrichment, color=gwas)) + 
  geom_segment(aes(x=annotation, xend=annotation, y=0, yend=Pct._h2), color="grey") +
  geom_segment(aes(x=annotation, xend=annotation, y=0, yend=Prop._SNPs), color="black") +
  coord_flip() + 
  labs(size=expression(h^{2}~enrichment), x="", y=expression("%"~h^{2}~explained)) +
  guides(color=F)
file.out <- "figs/fig_h2_annotations.complex.pdf"
# ggsave(file.out, width=10, height=5)



# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #


# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #


