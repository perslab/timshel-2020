############### SYNOPSIS ###################
### AIM: Make hypothalamus plot
# With dendrogram

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
# ========================== BMI GENESET ENRICHMENT =================== #
# ======================================================================= #

### ---> NOT USED (but should work)
# ### READ: enrichment results for MB+hypothalamus cell-types
# file.enrich <- here("src/publication/tables/table-es_enrichment.combined.csv") # con
# df.enrich <- read_csv(file.enrich)
# df.enrich <- df.enrich %>% mutate(annotation_fmt=annotation) # copy column | only needed for legacy reasons because of _join() operation

# ======================================================================= #
# ================== LOAD LDSC CTS RESULTS (multi GWAS) ================ #
# ======================================================================= #

### Read LDSC results
file.results <- here("results/cellect_ldsc/prioritization.csv")
df.ldsc_cts <- read_csv(file.results)
df.ldsc_cts <- format_cellect_ldsc_results(df.ldsc_cts)
df.ldsc_cts <- df.ldsc_cts %>% filter(specificity_id %in% get_scrna_seq_dataset_prefixes("hypo"))
df.ldsc_cts %>% count(specificity_id)

# ======================================================================= #
# ============================ BMI point plot =============================== #
# ======================================================================= #

### Create 'plotting' data frame
df.ldsc_cts <- df.ldsc_cts %>% filter(gwas == "BMI_UKBB_Loh2018")  # filter BMI results

### Add fdr_significant flag (within GWAS)
df.ldsc_cts <- df.ldsc_cts %>% mutate(fdr_significant = if_else(p.value <= 0.05/n(), true=T, false=F))

### Add columns
df.ldsc_cts <- df.ldsc_cts %>% mutate(p.value.mlog10 = -log10(p.value))

# ======================================================================= #
# ============================ PROCESS: add meta-data ====================== #
# ======================================================================= #

### Init
df.join <- df.ldsc_cts

### Add clean names
df.join <- df.join %>% mutate(annotation_fmt = utils.rename_annotations.hypothalamus(annotation, specificity_id, check_all_matches=T))

### Add meta-data
df.metadata <- get_metadata("hypothalamus")
df.join <- df.join %>% left_join(df.metadata %>% select(-annotation, -specificity_id), by="annotation_fmt")

### Add extra columns
df.join <- df.join %>% mutate(
  category=taxonomy_lvl1, # *DETERMINES WHAT IS COLORED BY*
  text=annotation_fmt
  # text=annotation_marker
  
)


# ======================================================================= #
# ============================ GET DENDROGRAM =========================== #
# ======================================================================= #

source("fig-annotation_clustering_hypo.R")
list.res <- get_hypothalamus_integrated_dendrogram(save_fig=FALSE, clustering_method="wards")
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
                                          # position="bottom" # put axis on top | works esi-ok
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
p.dendro.w_labels


# ======================================================================= #
# ============================ PLOT =============================== #
# ======================================================================= #

### Init
df.plot <- df.join

### Order annotations by tax1
df.plot <- df.plot %>% arrange(taxonomy_lvl1, annotation_fmt) %>% mutate(annotation_fmt = factor(annotation_fmt, levels=annotation_fmt))

### Order annotations by dendrogram
order.annotations <- labels(dend) # labels(dend) returns ordered labels | order.dendrogram(dend) returns ordered index.
df.plot <- df.plot %>% mutate(annotation_fmt = factor(annotation_fmt, levels=order.annotations))

### Colormap
names.category <- sort(unique(df.plot$category))
# colormap.category <- brewer.pal(name="Set1", n=9) # n=9 max for Set1 | "Set2"/"Dark2" ok for colorblind
# colormap.category <- colormap.category[1:length(names.category)] # select the colors we need
colormap.category <- colorspace::qualitative_hcl(n=length(names.category), palette = "Dark 3")
names(colormap.category) <- names.category
colormap.category


fdr_threshold <- 0.05/nrow(df.plot)
p.priori <- ggplot() +
  ### fdr line
  geom_hline(yintercept=-log10(fdr_threshold), linetype="dashed", color="darkgray") + 
  ### points - no enrichment
  geom_point(data=df.plot, aes(x=annotation_fmt, y=p.value.mlog10, fill=category, color=category), shape=21, stroke=0.1) +
  ### points - w. enrichment
  # geom_point(data=df.plot, aes(x=annotation_fmt, y=p.value.mlog10, fill=category, size=-log10(combined_rare_mendelian_obesity), color=category), shape=21, stroke=0.1) +
  # ^ shapes 21-24 have both stroke colour and a fill. The size of the filled part is controlled by size, the size of the stroke is controlled by stroke.
  # geom_point(data=df.plot %>% filter(combined_rare_mendelian_obesity<fdr_threshold), aes(annotation_fmt, y=p.value.mlog10, size=-log10(combined_rare_mendelian_obesity), fill=category), stroke=0.8, color="black", shape=21, show.legend=F) + 
  ### text/description
  geom_text(data=df.plot %>% filter(fdr_significant), aes(x=annotation_fmt, y=p.value.mlog10, label=text), size=rel(2), hjust=0, nudge_y=0.08, show.legend=F) +
  # geom_text(data=df.plot, aes(x=annotation_fmt, y=p.value.mlog10, label=text), size=rel(2), angle=90, vjust=0.22, hjust=-0.08, show.legend=F) +
  # ^ this works when no using coord_flip()
  ### axes
  labs(x="", y=expression(-log[10](P[S-LDSC])), size=expression(-log[10](P[enrichment])), fill="") +
  # coord
  coord_flip(clip = 'off') + # This keeps the labels from disappearing 
  ### color
  scale_color_manual(values=colormap.category) +
  scale_fill_manual(values=colormap.category) +
  ### guides
  guides(stroke=F, color=F,
         fill=guide_legend(order = 1),
         size=guide_legend(order = 2)
  ) + 
  ### theme
  theme_classic() + 
  theme(axis.text.y=element_text(size=rel(0.1)),
        axis.ticks.y=element_blank()) +
  # theme(axis.text.x=element_text(angle=65, hjust=1, size=rel(0.6))) + # no coord_flip
  theme(legend.position="top")
### Add margin to plot if displaying it directly (without pathwork)
p.priori <- p.priori + theme(plot.margin = unit(c(1,3,1,1), "cm")) # (t, r, b, l) widen margin
p.priori

file.out <- "figs/fig_celltypepriori.hypothalamus.priori.pdf"
# ggsave(plot=p.priori, filename=file.out, width=8, height=8)

# ======================================================================= #
# ======================== PATCHWORK: PRIORI + DENDROGRAM ==================== #
# ======================================================================= #

### Check alignment ok
p.dendro.w_labels + p.priori

### Make patchwork
p.patch <- p.dendro + p.priori + plot_spacer() + plot_layout(nrow=1, widths=c(0.3,1, 0.001))
p.patch

file.out <- "figs/fig_celltypepriori.hypothalamus.priori.with_dendrogram.pdf"
ggsave(plot=p.patch, filename=file.out, width=11, height=8)


# ======================================================================= #
# ============================ XXXXXXXXX =============================== #
# ======================================================================= #


# ======================================================================= #
# ============================ XXXXXXXXX =============================== #
# ======================================================================= #


# ======================================================================= #
# ============================ LEFTOVERS =============================== #
# ======================================================================= #
