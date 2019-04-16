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

# annotation <- "Pancreas.pancreatic A cell" # GCG
# gene_highlight <- "GCG"
# # "Pancreatic Alpha-cells (glucagon+)"

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
order.annotations <- rev(c("Hypothalamus\nArgp+ cells"="mb_agrp", "Striatum\nMedium Spiny Neurons"="mb_msn", "Liver\nHepatocytes"="tm_hepato", "Pancreas\nAlpha-cells"="tm_alpha"))
tmp_vector <- order.annotations # convenience selector. If you want a specific order of the result, this vector should contain (unique) ordered values
newnames <- names(order.annotations)
rename_vector <- newnames; names(rename_vector) <- tmp_vector
df.es <- df.es %>% mutate(annotation_fmt = recode_factor(annotation, !!!rename_vector)) # Use a named character vector to recode factors with unquote splicing. | REF: https://dplyr.tidyverse.org/reference/recode.html
df.es
# ======================================================================= #
# ============================ geom_density_ridges ========================= #
# ======================================================================= #


# ============================ TODO ========================= #
### write function to get 'ES gene' threshold for each ES metric (tstat, ges, si, specificity).
# --> this will be used for fill coloring of density.
### ADD SEGMENT FOR "ES GENE THRESHOLD"

### [OK] CHANGE AGRP SEGMENT TO A POINT
### [OK] ADD GCG AS FACET WRAP
###### LAYOUT
### [OK] Rename before faceting
### [*NOW*] plot "ESmu plot" with highlighted cell-types
### make patchwork

# ============================ facet_wrap ========================= #
df.plot <- df.es
df.plot <- df.plot %>% filter(es_weight > 0)
df.plot.genes_select <- df.plot %>% filter(dataset=="mb" & gene_name %in% c("AGRP", "DRD2") |
                                           dataset=="tm" & gene_name %in% c("GCG", "APOA2")
                                           )
# "RBFOX3" "FOS" "JUN"
p <- ggplot(df.plot) + 
  geom_density_ridges(aes(x=es_weight_log, y=es_metric_fmt), color="white", alpha=0.6) +
  geom_point(data=df.plot.genes_select, aes(x=es_weight_log, y=as.numeric(as.factor(es_metric_fmt)),
                                            color=gene_name)) +
  # scale_x_continuous(limits=c(-1,2), expand = c(0.01, 0)) + 
  labs(x=expression(Gene~ES[w]), y="ES metric", color="Gene") + 
  theme_ridges() + 
  facet_wrap(~annotation_fmt, nrow=1) +
  theme(strip.background=element_blank())
  
p



# ============================ MB plot ========================= #
df.plot <- df.es %>% filter(dataset == "mb")
df.plot <- df.plot %>% filter(es_weight > 0)
p <- ggplot(df.plot) + 
  geom_density_ridges(aes(x=es_weight_log, y=es_metric_fmt), color="white", alpha=0.6) +
  geom_segment(data=df.plot %>% filter(gene_name == "AGRP"), aes(x=es_weight_log, xend=es_weight_log, 
                                                                 y=as.numeric(as.factor(es_metric_fmt)), yend=as.numeric(as.factor(es_metric_fmt))+0.2), color="red") +
  # scale_x_continuous(limits=c(-1,2), expand = c(0.01, 0)) + 
  labs(x=expression(Gene~ES[w]), y="ES metric")
  # theme_ridges()
p


### THINGS TO USE
# - BUILD CUSTOM quantile_fun: how to match to the correct ESmetric? only values and props are passed
# - custom vertical line: https://stackoverflow.com/questions/48393348/add-custom-vertical-line-joyplots-ggridges
# - theme_ridges() \ theme_ridges(grid = FALSE, center = TRUE)



# ======================================================================= #
# ===================== TESTING geom_density_ridges ===================== #
# ======================================================================= #


p <- ggplot(df.plot, aes(x=es_weight_log_select, y=es_metric, fill=factor(..quantile..))) +
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 4, quantile_lines = TRUE
  ) +
  scale_fill_viridis(discrete = TRUE, name = "Quartiles")
p


ggplot(iris, aes(x = Sepal.Length, y = Species)) +
  geom_density_ridges(
    jittered_points = TRUE, quantile_lines = TRUE, scale = 0.9, alpha = 0.7,
    vline_size = 1, vline_color = "red", vline_alpha = 1,
    point_size = 0.4, point_alpha = 1,
    position = position_raincloud(adjust_vlines = TRUE)
  )


# default scales
ggplot(iris, aes(x=Sepal.Length, y=Species, fill = Species, color = Species)) +
  geom_density_ridges(
    aes(vline_color = Species, vline_linetype = Species),
    alpha = .4, quantile_lines = TRUE
  ) +
  theme_ridges()

# stat_density_ridges(mapping = NULL, data = NULL,
#                     geom = "density_ridges", position = "identity", na.rm = FALSE,
#                     show.legend = NA, inherit.aes = TRUE, bandwidth = NULL,
#                     from = NULL, to = NULL, jittered_points = FALSE,
#                     quantile_lines = FALSE, calc_ecdf = FALSE, quantiles = 4,
#                     quantile_fun = quantile, ...)



# ======================================================================= #
# ============================ PLOT 1: facet_wrap ========================= #
# ======================================================================= #

### 
df.plot <- df.es
# df.plot <- df.plot %>% filter(es_metric=="tstat")
# df.plot <- df.plot %>% filter(es_metric == "mean nES")
df.plot <- df.plot %>% filter(es_weight > 0)
p <- ggplot(df.plot, aes(x=es_weight)) + 
  geom_density() + 
  scale_x_log10() +
  facet_wrap(~es_metric, scales = "free")
p

# ======================================================================= #
# ============================ TMP COLOR density ========================= #
# ======================================================================= #


# =======================  ======================= #
### REF: https://stackoverflow.com/questions/3494593/shading-a-kernel-density-plot-between-two-points
set.seed(1)
draws <- rnorm(100)^2
q75 <- quantile(draws, .75)
q95 <- quantile(draws, .95)
dens <- density(draws)
dd <- with(dens,data.frame(x,y)) # 'dens' contains elements 'x' and 'y'
qplot(x,y,data=dd,geom="line")+
  geom_ribbon(data=subset(dd,x>q75 & x<q95),aes(ymax=y),ymin=0,
              fill="red",colour=NA,alpha=0.5)
