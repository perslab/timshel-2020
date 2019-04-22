############### SYNOPSIS ###################
### AIM: Make Campbell LDSC plot

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

source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))

# ======================================================================= #
# ============================ CAMPBELL LVL2 LDSC ============================ #
# ======================================================================= #

dataset_prefix <- "campbell_lvl2"
filter.gwas <- "BMI_UKBB_Loh2018"

# ================== LOAD LDSC CTS RESULTS (multi GWAS) ================ #

### Read LDSC results
file.results <- here("results/prioritization_celltypes--campbell_lvl2.multi_gwas.csv.gz")
df.ldsc_cts <- read_csv(file.results)

# =========================== FILTER GWAS =========================== #
df.ldsc_cts <- df.ldsc_cts %>% filter(gwas == filter.gwas)

# =========================== ADD METADATA =========================== #
df.metadata <- get_metadata(dataset_prefix)
df.ldsc_cts <- df.ldsc_cts %>% left_join(df.metadata, by="annotation")

### Clean annotation names (some meta-data for Campbell was not correctly provided by the authors)
# recode: x=new_name; name(x)=old_name
vector_rename <- c("n03"="n03.Th/Sst",
                   "n07"="n07.Arx/Nr5a2",
                   "n08"="n08.Th/Slc6a3",
                   "n27"="n27.Tbx19",
                   "n29.Nr5a1-Adcyap1"="n29.Nr5a1/Bdnf")
df.ldsc_cts <- df.ldsc_cts %>% mutate(annotation=recode(annotation, !!!vector_rename))

### Make new 'clean name' annotation 
tmp.split_str <- str_split(df.ldsc_cts$annotation, pattern="\\.") # list
df.ldsc_cts <- df.ldsc_cts %>% mutate(annotation_clean=tmp.split_str %>% purrr::map_chr(1))
### Add text column
df.ldsc_cts <- df.ldsc_cts %>% mutate(text=tmp.split_str %>% purrr::map_chr(2))
df.ldsc_cts <- df.ldsc_cts %>% mutate(text=str_replace_all(text,"-", "/")) # replace "-" with forwardslash ("/")
df.ldsc_cts <- df.ldsc_cts %>% mutate(text=str_replace_all(text,"_", " "))

### Clean taxonomy 
df.ldsc_cts <- df.ldsc_cts %>% mutate(taxonomy_lvl2=str_replace_all(taxonomy_lvl2, "_", " "))

### Rename Neuron to Neurons to match with mousebrain
# df.ldsc_cts <- df.ldsc_cts %>% mutate(taxonomy_lvl2=case_when(
#   taxonomy_lvl2 == "Neuron" ~ "Neurons",
#   TRUE ~ as.character(taxonomy_lvl2)
# )
# )

# ======================================================================= #
# ================================ PLOT LDSC ================================ #
# ======================================================================= #

### Add pvalue
df.plot <- df.ldsc_cts %>% mutate(p.value.mlog10 = -log10(p.value))
df.plot

### Colormap
factor_levels <- sort(unique(df.plot$taxonomy_lvl2))
colormap <- colorspace::qualitative_hcl(n=length(factor_levels), palette = "Dark 3")
names(colormap) <- factor_levels

fdr_threshold <- 0.05/nrow(df.plot)
p.priori <- ggplot() +
  ### fdr line
  geom_hline(yintercept=-log10(fdr_threshold), linetype="dashed", color="darkgray") + 
  ### points
  geom_point(data=df.plot, aes(x=annotation_clean, y=p.value.mlog10, color=taxonomy_lvl2)) +
  ### text/description
  geom_text(data=df.plot, aes(x=annotation_clean, y=p.value.mlog10, label=text), size=rel(2.5), hjust=0, nudge_y=0.08, show.legend=F) +
  ### axes
  labs(x="", y=expression(-log[10](P[S-LDSC])), color="") +
  # coord
  coord_flip(clip = 'off') + # This keeps the labels from disappearing 
  ### guides
  #...
  ### color
  scale_color_manual(values=colormap) +
  ### theme
  theme_classic() + 
  theme(axis.text.y=element_text(size=rel(0.6)),
        axis.ticks.y=element_blank()) +
  # theme(axis.text.x=element_text(angle=65, hjust=1, size=rel(0.6))) + # no coord_flip
  theme(legend.position="top")
### Add margin to plot if displaying it directly (without pathwork)
p.priori <- p.priori + theme(plot.margin = unit(c(1,2,1,1), "cm")) # (t, r, b, l) widen margin
p.priori

file.out <- "figs/fig_celltypepriori.campbell.priori.pdf"
ggsave(plot=p.priori, filename=file.out, width=8, height=8)
