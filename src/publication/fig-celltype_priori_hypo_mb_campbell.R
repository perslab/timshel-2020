############### SYNOPSIS ###################
### AIM: Make hypothalamus plot: Campbell+Mousebrain + BMI geneset cell-type enrichment

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

library(patchwork)

source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))



# ======================================================================= #
# ========================== CELL-TYPE BMI ENRICHMENT =================== #
# ======================================================================= #


### READ: all MB + Campbell cell-types
file.enrich <- here("results/es_enrichment--bmi_gene_lists.pvals.csv")
df.enrich <- read_csv(file.enrich)

### Make clean annotations (only needed for Campbell)
### Campbell annotations look like this because of Jon's code:
# n31.Fam19a2
# n32.Slc17a6.Trhr
# n33.unassigned.1.
# n34.unassigned.2.
tmp.split_str <- str_split(df.enrich$annotation, pattern="\\.") # list
df.enrich <- df.enrich %>% mutate(annotation_clean=tmp.split_str %>% purrr::map_chr(1))


# ======================================================================= #
# ================================ MOUSEBRAIN LDSC ================================ #
# ======================================================================= #

dataset_prefix <- "mousebrain_all"
filter.gwas <- "BMI_UKBB_Loh2018"

# ================== LOAD LDSC CTS RESULTS (multi GWAS) ================ #

### Read LDSC results
file.results <- here("results/prioritization_celltypes--mousebrain_all.multi_gwas.csv.gz")
df.ldsc_cts <- read_csv(file.results)

# =========================== FILTER GWAS =========================== #
df.ldsc_cts <- df.ldsc_cts %>% filter(gwas == filter.gwas)

# =========================== ADD METADATA =========================== #

df.metadata <- get_metadata(dataset_prefix)
df.ldsc_cts <- df.ldsc_cts %>% left_join(df.metadata, by="annotation")

### Add 'text' column
df.ldsc_cts <- df.ldsc_cts %>% mutate(text=str_split(Description, pattern=",") %>% purrr::map_chr(1))

### Add 'clean name' annotation (dummy for mb)
df.ldsc_cts <- df.ldsc_cts %>% mutate(annotation_clean=annotation)

# ======================== FILTER REGION =========================== #

filter.annotations <- get_annotations.mousebrain.hypothalamus()
df.ldsc_cts <- df.ldsc_cts %>% filter(annotation %in% filter.annotations)

# ======================== SELECT =========================== #
df.mb <- df.ldsc_cts

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
# df.ldsc_cts <- df.ldsc_cts %>% mutate(text=str_replace_all(text,"-", "\\\\")) # replace "-" with backslash ("\") [yes it is a little insane] | REF: https://github.com/STAT545-UBC/Discussion/issues/394



### Rename Neuron to Neurons to match with mousebrain
df.ldsc_cts <- df.ldsc_cts %>% mutate(taxonomy_lvl2=case_when(
  taxonomy_lvl2 == "Neuron" ~ "Neurons",
  TRUE ~ as.character(taxonomy_lvl2)
  )
)

# =========================== FILTER NEURONS =========================== #
# df.ldsc_cts <- df.ldsc_cts %>% filter(taxonomy_lvl2 == "Neurons")

# ======================== SELECT =========================== #
df.cb <- df.ldsc_cts

# ======================================================================= #
# ============================ MERGE AND PREP LDSC DATA ====================== #
# ======================================================================= #
list.bind <- list("Mousebrain"=df.mb %>% select(annotation, annotation_clean, p.value, category=Class, text, n_cells=NCells), 
     "Campbell"=df.cb %>% select(annotation, annotation_clean, p.value, category=taxonomy_lvl2, text, n_cells=n))
# TODO: make clean text column for each dataset before joining.

df.join <- bind_rows(list.bind, .id="dataset")

### Add pvalue
df.join <- df.join %>% mutate(p.value.mlog10 = -log10(p.value))

### Add enrichment data
df.join <- df.join %>% left_join(df.enrich %>% select(annotation_clean, combined_rare_mendelian_obesity), by="annotation_clean")

# ======================================================================= #
# ============================ GET DENDROGRAM =========================== #
# ======================================================================= #

source("fig-annotation_clustering_mb_campbell.R")
list.res <- get_mousebrain_campbell_integrated_dendrogram(save_fig=FALSE)
p.dendro.orig <- list.res[["plot"]]
dend <- list.res[["dendrogram"]]

# ======================================================================= #
# ============================ MODIFY DENDROGRAM =========================== #
# ======================================================================= #


# ============================ FUNCTION =========================== #
remove_geom <- function(ggplot2_object, geom_type) {
  # COPY FROM https://stackoverflow.com/a/46221054/6639640
  # ALTERNATIVE: https://rdrr.io/cran/ggedit/man/remove_geom.html
  # Delete layers that match the requested type.
  layers <- lapply(ggplot2_object$layers, function(x) {
    if (class(x$geom)[1] == geom_type) {
      NULL
    } else {
      x
    }
  })
  # Delete the unwanted layers.
  layers <- layers[!sapply(layers, is.null)]
  ggplot2_object$layers <- layers
  ggplot2_object
}

### TEST FUNCTION
# lapply(p.dendro.orig$layers, function(x) {print(class(x$geom)[1])})

# ============================ MODIFY =========================== #
n_leaf_nodes <- length(order.dendrogram(dend))

p.dendro <- p.dendro.orig # init
p.dendro <- p.dendro + theme_classic() # "reset" theme settings 
p.dendro <- p.dendro + coord_flip(clip="off",
                                  expand=F,
                                  xlim=c(-1,n_leaf_nodes-0.5) # OBS: not sure why the "-0.5" is needed, but it helps the alingment with LDSC plot.
                                  # ^ xlim is needed for perfect/correct alingment with ldsc plot
                                  # ^ xlim "zooms in" without removing data | REF https://stackoverflow.com/questions/25685185/limit-ggplot2-axes-without-removing-data-outside-limits-zoom
                                  # the initially plot is contructed with " p + xlim(c(-1,n_leaf_nodes+1))"
                                  )
# p.dendro <- p.dendro + xlim(c(0,n_leaf_nodes)) # this seems to help the alignment with LDSC if scale_y_continuous(position="bottom) is enabled
p.dendro <- p.dendro + scale_y_continuous(trans="reverse" # same as scale_y_reverse()
                                          # position="bottom" # put axis on top | works semi-ok
                                          )
p.dendro <- p.dendro + guides(color=F) # no guides for dataset color
p.dendro <- p.dendro + theme(axis.line.y=element_blank(), # remove all x-axis 
                             axis.title.y=element_blank(),
                             axis.ticks.y=element_blank(),
                             axis.text.y=element_blank())
p.dendro <- remove_geom(p.dendro, "GeomText") # remove old node text
p.dendro.w_labels <- p.dendro + geom_node_text(aes(filter=(leaf==TRUE), color=dataset, label=label), # add new node text [only needed temporarily to check that patchwork alings correctly]
                        size=rel(1.2), angle=0, hjust="left", nudge_y=0.01, 
                        show.legend=F)
p.dendro.w_labels <- p.dendro.w_labels + theme(plot.margin = unit(c(1,3,1,1), "cm")) # (t, r, b, l) widen margin
# p.dendro.w_labels

# ======================================================================= #
# ================================ PLOT LDSC ================================ #
# ======================================================================= #

### Order annotations by dendrogram
order.annotations <- labels(dend) # labels(dend) returns ordered labels | order.dendrogram(dend) returns ordered index.
order.annotations <- str_split(order.annotations, pattern="\\.") %>% purrr::map_chr(1) # clean names
df.join <- df.join %>% mutate(annotation_clean = factor(annotation_clean, levels=order.annotations))


### Colormap (*SHOULD BE THE SAME AS IN lib-annotation_clustering plot_es_dendrogram.mb_campbell()*)
dataset_names <- sort(unique(df.join$dataset))
colormap.dataset <- brewer.pal(name="Set1", n=9) # n=9 max for Set1 | "Set2"/"Dark2" ok for colorblind
colormap.dataset <- colormap.dataset[1:length(dataset_names)] # select the colors we need
names(colormap.dataset) <- dataset_names
# Campbell Mousebrain 
# "#E41A1C"  "#377EB8" 


fdr_threshold <- 0.05/nrow(df.join)
p.priori <- ggplot() +
  ### fdr line
  geom_hline(yintercept=-log10(fdr_threshold), linetype="dashed", color="darkgray") + 
  ### points
  geom_point(data=df.join, aes(x=annotation_clean, y=p.value.mlog10, fill=dataset, size=-log10(combined_rare_mendelian_obesity), color=dataset), shape=21, stroke=0.1) +
  # ^ shapes 21-24 have both stroke colour and a fill. The size of the filled part is controlled by size, the size of the stroke is controlled by stroke.
  geom_point(data=df.join %>% filter(combined_rare_mendelian_obesity<fdr_threshold), aes(annotation_clean, y=p.value.mlog10, size=-log10(combined_rare_mendelian_obesity), fill=dataset), stroke=0.8, color="black", shape=21, show.legend=F) + 
  ### text/description
  geom_text(data=df.join, aes(x=annotation_clean, y=p.value.mlog10, label=text), size=rel(2), hjust=0, nudge_y=0.08, show.legend=F) +
  # geom_text(data=df.join, aes(x=annotation_clean, y=p.value.mlog10, label=text), size=rel(2), angle=90, vjust=0.22, hjust=-0.08, show.legend=F) +
  # ^ this works when no using coord_flip()
  ### axes
  labs(x="", y=expression(-log[10](P[S-LDSC])), size=expression(-log[10](P[enrichment])), fill="") +
  # coord
  coord_flip(clip = 'off') + # This keeps the labels from disappearing 
  ### color
  scale_color_manual(values=colormap.dataset) +
  scale_fill_manual(values=colormap.dataset) +
  ### guides
  guides(stroke=F, color=F,
         fill=guide_legend(order = 1),
         size=guide_legend(order = 2)) + 
  ### theme
  theme_classic() + 
  theme(axis.text.y=element_text(size=rel(0.6)),
        axis.ticks.y=element_blank()) +
  # theme(axis.text.x=element_text(angle=65, hjust=1, size=rel(0.6))) + # no coord_flip
  theme(legend.position="top")
### Add margin to plot if displaying it directly (without pathwork)
p.priori <- p.priori + theme(plot.margin = unit(c(1,3,1,1), "cm")) # (t, r, b, l) widen margin
p.priori

file.out <- "figs/fig_celltypepriori.integrated_mb_campbell.priori.pdf"
ggsave(plot=p.priori, filename=file.out, width=8, height=8)



# ======================================================================= #
# ======================== PATCHWORK: PRIORI + DENDROGRAM ==================== #
# ======================================================================= #

### Check alignment ok
# p.dendro.w_labels + p.priori

### Make patchwork
p.patch <- p.dendro + p.priori + plot_spacer() + plot_layout(nrow=1, widths=c(0.3,1, 0.001))
p.patch

file.out <- "figs/fig_celltypepriori.integrated_mb_campbell.priori_with_dendrogram.pdf"
ggsave(plot=p.patch, filename=file.out, width=11, height=8)


### [OK DELETE] Debug alingment of plots:
# layer_scales(p.dendro)$x$range$range # get x-axis of barplot | REF: https://stackoverflow.com/a/35372274/6639640
# layer_scales(p.priori)$x$range$range


# =================================================================== #
# =========================== PLOT ENRICHMENT ======================= #
# =================================================================== #

### FILTER to match cell-types in df.join
df.enrich.filter <- df.enrich %>% filter(annotation_clean %in% df.join$annotation_clean)
df.enrich.filter

fdr_threshold <- 0.05/nrow(df.enrich.filter)
p.enrich <- ggplot() +
  geom_col(data=df.enrich.filter, aes(x=annotation_clean, y=-log10(combined_rare_mendelian_obesity), fill=dataset)) +
  ### extra
  geom_hline(yintercept=-log10(fdr_threshold), linetype="dashed", color="darkgray") + 
  ### axes
  labs(x="", y=expression(-log[10](P[enrichment]))) +
  # coord
  coord_flip(clip = 'off') + # This keeps the labels from disappearing 
  ### guides
  guides(fill=F) +
  ### color
  scale_fill_manual(values=colormap.dataset) +
  ### theme
  theme_classic() + 
  theme(axis.text.y=element_text(size=rel(0.6)))

p.enrich
file.out <- "figs/fig_celltype_geneset_enrichment.integrated_mb_campbell.pdf"
ggsave(plot=p.enrich, filename=file.out, width=8, height=8)


### Add margin to plot if displaying it directly (without pathwork)
p.enrich.mod <- p.enrich +
  scale_y_reverse() + 
  theme(axis.ticks.y = element_blank(), # wipe away all x-axis elements
        axis.text.y = element_blank(),
        axis.line.y = element_blank())
p.enrich.mod


# ======================================================================= #
# ======================== PATCHWORK: PRIORI + ENRIH ==================== #
# ======================================================================= #

p.patch <- p.enrich.mod + p.priori
p.patch

file.out <- "figs/fig_celltypepriori.integrated_mb_campbell.priori_with_enrichment.pdf"
ggsave(plot=p.patch, filename=file.out, width=8, height=8)



# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #


# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #


# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #









# ======================================================================= #
# ===================== LEFTOVER: no point size enrichment ===================== #
# ======================================================================= #
# WORKS!

# fdr_threshold <- 0.05/nrow(df.join)
# p.priori <- ggplot() +
#   geom_point(data=df.join, aes(x=annotation_clean, y=p.value.mlog10, color=dataset, size=-log10(combined_rare_mendelian_obesity))) +
#   ### text/description
#   # geom_text(data=df.join, aes(x=annotation_clean, y=p.value.mlog10, label=text), size=rel(2), angle=90, vjust=0.22, hjust=-0.08, show.legend=F) +
#   # ^ this works when no using coord_flip()
#   geom_text(data=df.join, aes(x=annotation_clean, y=p.value.mlog10, label=text), size=rel(2), hjust=0, nudge_y=0.08, show.legend=F) +
#   ### extra
#   geom_hline(yintercept=-log10(fdr_threshold), linetype="dashed", color="darkgray") + 
#   ### axes
#   labs(x="", y=expression(-log[10](P[S-LDSC]))) +
#   # coord
#   coord_flip(clip = 'off') + # This keeps the labels from disappearing 
#   ### guides
#   #...
#   ### color
#   scale_color_manual(values=colormap.dataset) +
#   ### theme
#   theme_classic() + 
#   theme(axis.text.y=element_text(size=rel(0.6))) +
#   # theme(axis.text.x=element_text(angle=65, hjust=1, size=rel(0.6))) + # no coord_flip
#   theme(legend.position="top")
# ### Add margin to plot if displaying it directly (without pathwork)
# p.priori
# 
# p.priori.margin <- p.priori + theme(plot.margin = unit(c(1,3,1,1), "cm")) # (t, r, b, l) widen margin
# p.priori.margin
