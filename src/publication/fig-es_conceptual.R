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


library(ggridges)
library(patchwork)

source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))



# ======================================================================= #
# ================================ FUNCTIONS ============================ #
# ======================================================================= #

get_df_es_from_multiple_es_obj <- function(list.data) {
  list.res <- list()
  for (i in seq_along(list.data)) {
    annotation <- list.data[[i]][["annotation"]]
    df.es.tmp <- get_es.annotation_centric(list.data[[i]][["data_obj"]], annotation=annotation)
    df.es.tmp <- df.es.tmp %>% mutate(dataset=list.data[[i]][["dataset"]]) # add dataset
    list.res[[annotation]] <- df.es.tmp
  }
  df.es <- bind_rows(list.res, .id="annotation")
}


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
# ================================ CONSTANTS ============================ #
# ======================================================================= #

### recode_factor() allows us to rename the factor values but KEEP the ORDER defined by the renaming vector
### Use a named character vector to recode factors with unquote splicing. | REF: https://dplyr.tidyverse.org/reference/recode.html
### x=new_names; names(x)=old_names

rename_vector.es_metric <- rev(c(
  "tstat"="DET", 
  "ges"="GES", 
  "specificity"="EP",
  "si"="NSI", 
  "es_mu"="ESmu"))

rename_vector.annotation_to_annotation_fmt <- rev(c(
  "DEINH6"="Hypothalamus\nArgp+ neurons",
  "MSN3"="Striatum\nMedium spiny neurons",
  "Pancreas.pancreatic_A_cell"="Pancreas\nAlpha-cells",
  "Liver.hepatocyte"="Liver\nHepatocytes"))

# we need this variable on df.es to filter df.es_threshold
filter.annotations <- names(rename_vector.annotation_to_annotation_fmt)

# ======================================================================= #
# ============================ DEFINE DATA TO PLOT ======================= #
# ======================================================================= #

list.data <- list()
list.data[[1]] <- list(
  "annotation"="DEINH6",
  "dataset"="mb",
  "data_obj"=sem_obj.mb)
list.data[[2]] <- list(
  "annotation"="MSN3",
  "dataset"="mb",
  "data_obj"=sem_obj.mb)
list.data[[3]] <- list(
  "annotation"="Pancreas.pancreatic_A_cell",
  "dataset"="tm",
  "data_obj"=sem_obj.tm)
list.data[[4]] <- list(
  "annotation"="Liver.hepatocyte",
  "dataset"="tm",
  "data_obj"=sem_obj.tm)


# ======================================================================= #
# ======================== EXTRACT AND PROCESS DATA ===================== #
# ======================================================================= #

df.es <- get_df_es_from_multiple_es_obj(list.data)

### Exclude expr_mean
df.es <- df.es %>% filter(!es_metric %in% c("expr_mean"))

### Add 'pretty' ES metrics and order factor
df.es <- df.es %>% mutate(es_metric_fmt = recode_factor(es_metric, !!!rename_vector.es_metric))

### Add 'pretty' annotation and order factor
df.es <- df.es %>% mutate(annotation_fmt = recode_factor(annotation, !!!rename_vector.annotation_to_annotation_fmt))

# ======================================================================= #
# ============================ GET ESw thresholds ========================= #
# ======================================================================= #

df.es_threshold.mb <- get_es_w_threshold(sem_obj.mb, threshold_pval=0.05)
df.es_threshold.tm <- get_es_w_threshold(sem_obj.tm, threshold_pval=0.05)
df.es_threshold <- bind_rows(df.es_threshold.mb, df.es_threshold.tm) # ok to bind because annotation are unique

### process
df.es_threshold <- df.es_threshold %>% gather(key="es_metric", value="es_w_threshold", -annotation)

### filter annotations [important]
df.es_threshold <- df.es_threshold %>% filter(annotation %in% filter.annotations)

### Add 'pretty' annotation names
df.es_threshold <- df.es_threshold %>% mutate(annotation_fmt = recode_factor(annotation, !!!rename_vector.annotation_to_annotation_fmt))

### Add 'pretty' ES metric names
### ***IMPORTANT*** we set levels manually to match df.es levels that includes ESmu [because we later use as.numeric(es_metric_fmt)]
### Alternative is to add 1 in the plot
df.es_threshold <- df.es_threshold %>% mutate(es_metric_fmt = recode_factor(es_metric, !!!rename_vector.es_metric))
# ^ could not get recode_factor() .missing=NA_character_ or .default=NA_character_ to add ESmu
df.es_threshold <- df.es_threshold %>% mutate(es_metric_fmt = factor(es_metric_fmt, levels=levels(df.es$es_metric_fmt)))


# ======================================================================= #
# ========================== PLOT geom_density_ridges =================== #
# ======================================================================= #

### Init
df.plot <- df.es

### Filter 1: we can only plot positive weights
df.plot <- df.plot %>% filter(es_weight > 0)

# Select genes to plot for each dataset
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
  geom_point(data=df.plot.genes_select, aes(x=es_weight, 
                                            y=as.numeric(as.factor(es_metric_fmt)),
                                            color=gene_name)) +
  geom_segment(data=df.es_threshold, aes(x=es_w_threshold, xend=es_w_threshold,
                                         y=as.numeric(as.factor(es_metric_fmt)),
                                         yend=as.numeric(as.factor(es_metric_fmt))+0.05
                                         ), color="black") + 
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
ggsave(plot=p, filename=file.out, width=12, height=4)

# ============================ TODO ========================= #
### [OK] ADD SEGMENT FOR "ES GENE THRESHOLD"
### [OK] plot "ESmu plot" with highlighted cell-types
### [!!] make patchwork: ggridges + ESmu plot
### [!!] JUST MAKE QUICK COMPILATION AND REFINE LATER [E.g. colors of annotations]

# ======================================================================= #
# ================================= PLOT ESmu =========================== #
# ======================================================================= #

source(here("src/lib/load_functions.R")) # load sc-genetics library

# ================================= MB PLOTS ============================= #
annotations_highlight.mb <- c("DEINH6", "MSN3") 
names(annotations_highlight.mb) <- recode(annotations_highlight.mb, !!!rename_vector.annotation_to_annotation_fmt)
annotations_colormap.mb <- c("DEINH6"="black", "MSN3"="black")

### AGRP
genes_select <- c("AGRP")
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj.mb, genes_select=genes_select, es_metric="es_mu")
p.mb_agrp <- plot_es.gene_centric.single_es_metric(df.es.gene, annotations_highlight=annotations_highlight.mb, annotations_colormap=annotations_colormap.mb)
p.mb_agrp

### MSN
genes_select <- c("DRD2")
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj.mb, genes_select=genes_select, es_metric="es_mu")
p.mb_msn <- plot_es.gene_centric.single_es_metric(df.es.gene, annotations_highlight=annotations_highlight.mb, annotations_colormap=annotations_colormap.mb)
p.mb_msn

# ================================= TM PLOTS ============================= #
annotations_highlight.tm <- c("Pancreas.pancreatic_A_cell","Liver.hepatocyte")
names(annotations_highlight.tm) <- recode(annotations_highlight.tm, !!!rename_vector.annotation_to_annotation_fmt)
annotations_colormap.tm <- c("Pancreas.pancreatic_A_cell"="black","Liver.hepatocyte"="black")

### Alpha
genes_select <- c("GCG")
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj.tm, genes_select=genes_select, es_metric="es_mu")
p.tm_alhpa <- plot_es.gene_centric.single_es_metric(df.es.gene, annotations_highlight=annotations_highlight.tm, annotations_colormap=annotations_colormap.tm)
p.tm_alhpa


### Liver
genes_select <- c("APOA2")
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj.tm, genes_select=genes_select, es_metric="es_mu")
p.tm_hepato <- plot_es.gene_centric.single_es_metric(df.es.gene, annotations_highlight=annotations_highlight.tm, annotations_colormap=annotations_colormap.tm)
p.tm_hepato

# ================================ PATCHWORK ES MU ============================ #

# p.patch_esmu <- p.mb_agrp + p.mb_msn + p.tm_alhpa + p.tm_hepato + plot_layout(nrow=1)
p.patch_esmu <- p.tm_hepato + p.tm_alhpa + p.mb_msn + p.mb_agrp + plot_layout(nrow=1) 
# ^ PLOTS SHOULD GO IN THE CORRECT ORDER as determined by rename_vector.annotation_to_annotation_fmt
p.patch_esmu <- p.patch_esmu & labs(x="", y="")
p.patch_esmu <- p.patch_esmu & scale_y_continuous(expand=c(0,0), limits=c(0,1))
p.patch_esmu <- p.patch_esmu & theme(axis.ticks.x = element_blank(), # wipe away all x-axis elements
                                       axis.text.x = element_blank(),
                                       axis.line.x = element_blank())
# p.patch_esmu

# ======================================================================= #
# ================================ PATCHWORK COMBINED ============================ #
# ======================================================================= #

p.patch.comb <- p + p.patch_esmu + plot_layout(ncol=1, heights=c(1,0.5))

file.out <- "figs/fig_es_conceptual.pdf"
ggsave(plot=p.patch.comb, filename=file.out, width=12, height=6)





# ======================================================================= #
# ================================= LEFTOVERS =============================== #
# ======================================================================= #

### Munging of df.es
# df.es <- df.es %>% 
#   filter(!es_metric %in% c("expr_mean")) %>%
#   mutate(
#   es_weight_pseudo_log=sign(es_weight)*log10(abs(es_weight+1)),
#   es_weight_log=log10(es_weight),
#   es_weight_log_select=case_when(
#     es_metric=="tstat" ~ log10(es_weight),
#     es_metric=="ges" ~ log10(es_weight),
#     TRUE ~ es_weight
#   )
# )



# ================================= OLD MANUAL WAY OF GETTING df.es - WORKS =============================== #
# rename_vector.data_key_to_annotation_fmt <- rev(c("mb_agrp"="Hypothalamus\nArgp+ neurons",
#                                                   "mb_msn"="Striatum\nMedium spiny neurons",
#                                                   "tm_hepato"="Liver\nHepatocytes",
#                                                   "tm_alpha"="Pancreas\nAlpha-cells"))
# 

# filter.annotations <- c("DEINH6","MSN3","Pancreas.pancreatic_A_cell","Liver.hepatocyte")
# rename_vector.data_key_to_annotation <- c("DEINH6"="mb_agrp",
#                                           "MSN3"="mb_msn",
#                                           "Pancreas.pancreatic_A_cell"="tm_alpha",
#                                           "Liver.hepatocyte"="tm_hepato")


# ### Get data for different annotations
# df.es.mb_agrp <- get_es.annotation_centric(sem_obj.mb, annotation="DEINH6")
# df.es.mb_msn <- get_es.annotation_centric(sem_obj.mb, annotation="MSN3") # or MSN2, MSN4
# df.es.tm_alpha <- get_es.annotation_centric(sem_obj.tm, annotation="Pancreas.pancreatic_A_cell")
# df.es.tm_hepato <- get_es.annotation_centric(sem_obj.tm, annotation="Liver.hepatocyte")
# 
# ### Combine
# list.df <- list(mb_agrp=df.es.mb_agrp,
#                 mb_msn=df.es.mb_msn,
#                 tm_alpha=df.es.tm_alpha,
#                 tm_hepato=df.es.tm_hepato)
# df.es <- bind_rows(list.df, .id="data_key")
# 
# ### Add dataset variable
# df.es <- df.es %>% mutate(dataset = if_else(str_starts(data_key,pattern="mb_"), "mb", "tm"))

### Add 'real' annotation name (order does not matter)
# df.es <- df.es %>% mutate(annotation = recode_factor(data_key, !!!rename_vector.data_key_to_annotation))


# # ======================================================================= #
# # =========================== PLOT ESmu - WITHOUT RENAMING [OK DELETE] =============== #
# # ======================================================================= #
# 
# # ================================= MB PLOTS ============================= #
# 
# ### AGRP
# df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj.mb, 
#                                                    genes_select=c("AGRP"), 
#                                                    es_metric="es_mu")
# p.mb_agrp <- plot_es.gene_centric.single_es_metric(df.es.gene, 
#                                                    annotations_highlight=c("DEINH6", "MSN3"),
#                                                    annotations_colormap=c("DEINH6"="black", "MSN3"="black"))
# p.mb_agrp
# 
# ### MSN
# df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj.mb, 
#                                                    genes_select=c("DRD2"), 
#                                                    es_metric="es_mu")
# p.mb_msn <- plot_es.gene_centric.single_es_metric(df.es.gene, 
#                                                   annotations_highlight=c("DEINH6", "MSN3"),
#                                                   annotations_colormap=c("DEINH6"="black", "MSN3"="black"))
# p.mb_msn
# 
# # ================================= TM PLOTS ============================= #
# 
# 
# ### Alpha
# df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj.tm, 
#                                                    genes_select=c("GCG"), 
#                                                    es_metric="es_mu")
# p.tm_alhpa <- plot_es.gene_centric.single_es_metric(df.es.gene, 
#                                                     annotations_highlight=c("Pancreas.pancreatic_A_cell","Liver.hepatocyte"),
#                                                     annotations_colormap=c("Pancreas.pancreatic_A_cell"="black","Liver.hepatocyte"="black"))
# p.tm_alhpa
# 
# 
# ### Liver
# df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj.tm, 
#                                                    genes_select=c("APOA2"), 
#                                                    es_metric="es_mu")
# p.tm_hepato <- plot_es.gene_centric.single_es_metric(df.es.gene, 
#                                                      annotations_highlight=c("Pancreas.pancreatic_A_cell","Liver.hepatocyte"),
#                                                      annotations_colormap=c("Pancreas.pancreatic_A_cell"="black","Liver.hepatocyte"="black"))
# p.tm_hepato
