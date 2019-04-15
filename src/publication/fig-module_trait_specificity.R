############### SYNOPSIS ###################
# Module trait specificity


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
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))

# ======================================================================= #
# ================================ READ DATA =============================== #
# ======================================================================= #

file.data <- here("results/prioritization_modules--fdr_sign_celltypes.multi_gwas.csv.gz")
df.ldsc_cts <-  read_csv(file.data)

# ======================================================================= #
# ============================= PRE-PROCESS DATA ========================= #
# ======================================================================= #

# ======================================================================= #
# ================================ PLOT =============================== #
# ======================================================================= #


### SELECTED modules
filter.module <- "lavenderblush" # "lightpink3"
df.plot <- df.ldsc_cts %>% select(gwas, p.value=filter.module) # select and rename

### SELECTED gwas
filter.gwas <- utils.get_gwas_ids_for_meta_analysis()

### SELECT and order GWAS factor
df.plot <- df.plot %>% filter(gwas %in% filter.gwas)
order.gwas <- df.plot %>% arrange(p.value) %>% pull(gwas)
df.plot <- df.plot %>% mutate(gwas = factor(gwas, levels=order.gwas))
### ADD formated gwas name
tmp_gwas_vector <- df.plot$gwas # convenience selector. If you want a specific order of the result, this vector should contain (unique) ordered values
newnames <- utils.rename_gwas(tmp_gwas_vector, style="fullname") # "fullname_author_year","fullname","abrv_author_year","abrv_year","abrv" 
rename_vector <- newnames; names(rename_vector) <- tmp_gwas_vector
df.plot <- df.plot %>% mutate(gwas_fmt = recode_factor(gwas, !!!rename_vector))

### ADD columns
df.plot <- df.plot %>% mutate(p.value.mlog10=-log10(p.value))

### PLOT
fdr_threshold <- 0.05/nrow(df.plot)
fdr_threshold

### Lollipop plot
p <- ggplot(df.plot, aes(x=gwas_fmt, y=p.value.mlog10)) + 
  geom_hline(yintercept=-log10(fdr_threshold), color="gray", linetype="dashed") +
  # geom_col(data=df.plot %>% filter(p.value<fdr_threshold), fill="red", position=position_dodge()) + 
  geom_segment(aes(x=gwas_fmt, xend=gwas_fmt, y=0, yend=p.value.mlog10), color="grey") +
  geom_point(color="gray", size=5) + 
  geom_point(data=df.plot %>% filter(p.value <= fdr_threshold), color="black", size=5) + 
  ### color by pvalue
  # geom_point(aes(color=p.value.mlog10), size=5) + 
  # scale_color_viridis_c(direction=-1) + 
  labs(x="", y=expression(-log[10](P[S-LDSC]))) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size=rel(0.9))) + 
  guides(color=FALSE) +
  coord_cartesian(clip="off") + 
  theme(plot.margin = unit(c(2,1,1,2), "cm")) # trbl
p

file.out <- "figs/fig_module_trait_specificity.pdf"
ggsave(file.out, plot=p, width=15, height=6)

### BARPLOT
# p <- ggplot(df.plot, aes(x=gwas_fmt, y=p.value.mlog10)) + 
#   geom_col(position=position_dodge()) + 
#   geom_col(data=df.plot %>% filter(p.value<fdr_threshold), fill="red", position=position_dodge()) + 
#   geom_hline(yintercept=-log10(fdr_threshold), color="gray", linetype="dashed") +
#   labs(x="", y=expression(-log[10](P[S-LDSC]))) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   theme(plot.margin = unit(c(1,1,1,2), "cm"))
# p

