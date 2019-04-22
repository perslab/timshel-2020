############### SYNOPSIS ###################
### AIM: Plot distributions of ES metrics.

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

library(gghighlight)
library(ggrepel)
library(viridis)

library(ggridges)


source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))

# ======================================================================= #
# ================================ LOAD DATA ============================ #
# ======================================================================= #

load(here("out/es/mousebrain_all.es_obj.RData"))
sem_obj.mb <- sem_obj
load(here("out/es/tabula_muris.es_obj.RData"))
sem_obj.tm <- sem_obj


# ======================================================================= #
# ======================== WIKI: finding marker genes ===================== #
# ======================================================================= #

### Look at highly expressed genes
# df.x <- df.es.tm_alpha %>% filter(es_metric=="expr_mean")

### Marker genes
# https://www.proteinatlas.org/humanproteome/tissue
# https://www.proteinatlas.org/humanproteome/tissue/heart
# Liver.hepatocyte - Liver: APOA2
# Heart.smooth_muscle_cell / Heart.cardiac_muscle_cell
# Heart: MYH6, ACTC1 and TNNI3 (examples of members of the myosin, actin and troponin families solely expressed in heart muscle)


# ======================================================================= #
# ======================== EXTRACT AND PROCESS DATA ===================== #
# ======================================================================= #

### Get data
df.es.mb_agrp <- get_es.annotation_centric(sem_obj.mb, annotation="DEINH6")
df.es.mb_msn <- get_es.annotation_centric(sem_obj.mb, annotation="MSN3") # or MSN2, MSN4
df.es.tm_alpha <- get_es.annotation_centric(sem_obj.tm, annotation="Pancreas.pancreatic_A_cell")
df.es.tm_hepato <- get_es.annotation_centric(sem_obj.tm, annotation="Liver.hepatocyte")

list.df <- list(mb_agrp=df.es.mb_agrp,
                mb_msn=df.es.mb_msn,
                tm_alpha=df.es.tm_alpha,
                tm_hepato=df.es.tm_hepato)
df.es <- bind_rows(list.df, .id="annotation")

### Munge data
df.es <- df.es %>% mutate(dataset = if_else(str_starts(annotation,pattern="mb_"), "mb", "tm"))
df.es <- df.es %>% 
  filter(!es_metric %in% c("expr_mean")) %>%
  mutate(
  es_weight_pseudo_log=sign(es_weight)*log10(abs(es_weight+1)),
  es_weight_log=log10(es_weight),
  es_weight_log_select=case_when(
    es_metric=="tstat" ~ log10(es_weight),
    es_metric=="ges" ~ log10(es_weight),
    TRUE ~ es_weight
  )
)

### Rename + order ES metrics
### IMPORANTLY recode_factor() allows us to rename the factor values but KEEP the ORDER defined by filter.gwas
order.es_metric <- rev(c(DET="tstat", GES="ges", EP="specificity", NSI="si", ESmu="es_mu"))
tmp_vector <- order.es_metric # convenience selector. If you want a specific order of the result, this vector should contain (unique) ordered values
newnames <- names(order.es_metric)
rename_vector <- newnames; names(rename_vector) <- tmp_vector
df.es <- df.es %>% mutate(es_metric_fmt = recode_factor(es_metric, !!!rename_vector)) # Use a named character vector to recode factors with unquote splicing. | REF: https://dplyr.tidyverse.org/reference/recode.html

### Rename + order annotation
order.annotations <- rev(c("Hypothalamus\nArgp+ neurons"="mb_agrp", "Striatum\nMedium spiny neurons"="mb_msn", "Liver\nHepatocytes"="tm_hepato", "Pancreas\nAlpha-cells"="tm_alpha"))
tmp_vector <- order.annotations # convenience selector. If you want a specific order of the result, this vector should contain (unique) ordered values
newnames <- names(order.annotations)
rename_vector <- newnames; names(rename_vector) <- tmp_vector
df.es <- df.es %>% mutate(annotation_fmt = recode_factor(annotation, !!!rename_vector)) # Use a named character vector to recode factors with unquote splicing. | REF: https://dplyr.tidyverse.org/reference/recode.html



# ======================================================================= #
# ============================ GET ESw thresholds ========================= #
# ======================================================================= #

df.es_threshold.mb <- get_es_w_threshold(sem_obj.mb, threshold_pval=0.05)
df.es_threshold.tm <- get_es_w_threshold(sem_obj.tm, threshold_pval=0.05)
df.es_threshold <- bind_rows(df.es_threshold.mb, df.es_threshold.tm) # ok to bind because annotation are unique

# ======================================================================= #
# ========================== PLOT geom_density_ridges =================== #
# ======================================================================= #

### Pre-process
df.plot <- df.es
df.plot <- df.plot %>% filter(es_weight > 0)
df.plot.genes_select <- df.plot %>% filter(dataset=="mb" & gene_name %in% c("AGRP", "DRD2") |
                                           dataset=="tm" & gene_name %in% c("GCG", "APOA2")
                                           ) # "RBFOX3" "FOS" "JUN"


### Define gene color mapping
factor_levels <- unique(df.plot.genes_select$gene_name)
colormap <- colorspace::qualitative_hcl(n=length(factor_levels), palette = "Dark 3")
names(colormap) <- factor_levels

### LOG-scale
p <- ggplot(df.plot) + 
  geom_density_ridges(aes(x=es_weight, y=es_metric_fmt), color="white", alpha=0.6) +
  geom_point(data=df.plot.genes_select, aes(x=es_weight, y=as.numeric(as.factor(es_metric_fmt)),
                                            color=gene_name)) +
  scale_x_continuous(trans=scales::log10_trans(), # same as "log10"
                     breaks = c(1e-4, 1e-2, 1, 1e2, 1e4),
                     # breaks = scales::trans_breaks("log10", function(x){10^x}),
                     labels = scales::trans_format("identity", function(x){formatC(x, format="fg", big.mark=",")})
                     # ^ scales::trans_format: nice formated y-labels ala 1, 100, 1,000, 100,000 
                     # labels = scales::trans_format("log10", scales::math_format(10^.x))
                     ) + 
  labs(x=expression(Gene~ES[w]), y="ES metric", color="Gene") + 
  scale_color_manual(values=colormap) + 
  facet_wrap(~annotation_fmt, nrow=1) +
  theme_classic() + 
  # theme_ridges() + 
  theme(strip.background=element_blank()) + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank())

p

file.out <- "figs/fig_es_conceptual.ggridges.pdf"
ggsave(plot=p, filename=file.out, width=10, height=5)

# ============================ TODO ========================= #
### [!!] ADD SEGMENT FOR "ES GENE THRESHOLD"
### [OK] plot "ESmu plot" with highlighted cell-types
### [!!] make patchwork: ggridges + ESmu plot
### [!!] JUST MAKE QUICK COMPILATION AND REFINE LATER [E.g. colors of annotations]

# ======================================================================= #
# ========================== PLOT ESmu =================== #
# ======================================================================= #


### AGRP
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj.mb, 
                                                   genes_select=c("AGRP"), 
                                                   es_metric="es_mu")
p.mb_agrp <- plot_es.gene_centric.single_es_metric(df.es.gene, 
                                           annotations_highlight="DEINH6",
                                           annotations_colormap=NULL)
p.mb_agrp

### MSN
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj.mb, 
                                                   genes_select=c("DRD2"), 
                                                   es_metric="es_mu")
p.mb_msn <- plot_es.gene_centric.single_es_metric(df.es.gene, 
                                                   annotations_highlight="MSN3",
                                                   annotations_colormap=NULL)
p.mb_msn

### Alpha
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj.tm, 
                                                   genes_select=c("GCG"), 
                                                   es_metric="es_mu")
p.tm_alhpa <- plot_es.gene_centric.single_es_metric(df.es.gene, 
                                                  annotations_highlight="Pancreas.pancreatic_A_cell",
                                                  annotations_colormap=NULL)
p.tm_alhpa


### Liver
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj.tm, 
                                                   genes_select=c("APOA2"), 
                                                   es_metric="es_mu")
p.tm_hepato <- plot_es.gene_centric.single_es_metric(df.es.gene, 
                                                    annotations_highlight="Liver.hepatocyte",
                                                    annotations_colormap=NULL)
p.tm_hepato

# ======================================================================= #
# ================================ PATCHWORK ============================ #
# ======================================================================= #
# p.mb_agrp
# p.mb_msn
# p.tm_alhpa
# p.tm_hepato

p.mb_agrp + p.mb_msn + p.tm_alhpa + p.tm_hepato + plot_layout(nrow=1)

# ======================================================================= #
# ================================= XXXXX =============================== #
# ======================================================================= #

