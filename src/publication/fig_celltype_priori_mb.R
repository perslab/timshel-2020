############### SYNOPSIS ###################
# Make mousebrain main fig


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

dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib/"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library

library(patchwork)

setwd(here("src/publication"))

# ======================================================================= #
# ============================ PARAMETERS =============================== #
# ======================================================================= #

# gwas <- "BMI_UPDATE_Yengo2018"
gwas <- "BMI_UKBB_Loh2018"
dataset_prefix <- "mousebrain_all"
genomic_annotation_prefix <- get_genomic_annotation_prefix(dataset_prefix)


# ======================================================================= #
# ============================ LOAD CLUSTER ORDER =============================== #
# ======================================================================= #

file.mb_dendrogram <- here("data/expression/mousebrain/mousebrain_table_s3-cluster_dendrogram.csv")
df.mb_dendrogram <- read_csv(file.mb_dendrogram)
df.mb_dendrogram <- df.mb_dendrogram %>% select(linnarson_order=`Dendrogram order`, annotation=`Cluster name`) # %>% mutate(linnarson_order = linnarson_order+1)

file.clustering_order <- here(sprintf("data/genes_cell_type_specific/%s.hclust_order.csv", dataset_prefix))
df.clustering_order <- read_csv(file.clustering_order)

# ======================================================================= #
# =============================== LOAD LDSC CTS RESULTS ================================= #
# ======================================================================= #

file.ldsc_cts <- sprintf("/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/%s__%s.cell_type_results.txt", genomic_annotation_prefix, gwas)
df.ldsc_cts <- load_ldsc_cts_results(file.ldsc_cts, dataset_prefix)
df.ldsc_cts <- df.ldsc_cts %>% filter(sem=="sem_mean")

# ======================================================================= #
# =============================== METADATA ================================= #
# ======================================================================= #

df.metadata <- get_metadata(dataset_prefix)
df.ldsc_cts <- df.ldsc_cts %>% left_join(df.metadata, by="annotation") # add meta data
df.ldsc_cts <- df.ldsc_cts %>% left_join(df.clustering_order, by="annotation") # add cluster order
df.ldsc_cts <- df.ldsc_cts %>% left_join(df.mb_dendrogram, by="annotation") # add cluster order

# ======================================================================= #
# ================================ PLOT: cell priori ================================= #
# ======================================================================= #

### Copy
df.plot <- df.ldsc_cts

### Rename
df.plot <- df.plot %>% mutate(Region = case_when(
  Region == "Midbrain dorsal" ~ "Midbrain",
  Region == "Midbrain dorsal,Midbrain ventral" ~ "Midbrain",
  Region == "Hippocampus,Cortex" ~ "Hippocampus/Cortex",
  TRUE ~ as.character(Region))
)

### Reorder
# df.plot <- df.plot %>% mutate(annotation=factor(annotation, levels=annotation[order(cluster_order)])) # hclust order USE ME!
# df.plot <- df.plot %>% mutate(annotation=factor(annotation, levels=annotation[order(linnarson_order)])) # linnarson order
df.plot <- df.plot %>% mutate(annotation=factor(annotation, levels=annotation[order(TaxonomyRank4)])) # linnarson order


fdr_threshold <- 0.05/nrow(df.plot)
p <- ggplot(df.plot, aes(x=annotation, y=-log10(p.value))) +
  geom_point(color="gray") + 
  geom_point(data=df.plot %>% filter(fdr_significant), aes(x=annotation, y=-log10(p.value), color=Region)) +
  # geom_segment(data=df.plot %>% filter(fdr_significant), aes(x=annotation, xend=annotation, y=0, yend=-log10(p.value)), color="gray", alpha=0.5) + # lollipop
  ggrepel::geom_text_repel(data=df.plot %>% filter(fdr_significant), aes(x=annotation, y=-log10(p.value), label=annotation, color=Region), hjust = 0, nudge_x = 1.5) + # OBS geom_text_repel
  geom_hline(yintercept=-log10(fdr_threshold), color="red") + 
  labs(x="Cell Type", y=expression(-log[10](P))) +
  ggtitle(sprintf("%s - %s", gwas, dataset_prefix)) +
  coord_flip() +
  theme_classic() + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  theme(legend.position="bottom")
p
# file.out <- sprintf("fig_1a.%s.%s.sem_mean.color_by_class.pdf", dataset_prefix, "BMI_UKBB_Loh2018")
# ggsave(p, filename=file.out, width=20, height=8)

# ======================================================================= #
# ============================ Annotations =============================== #
# ======================================================================= #

### Class
p_anno_class <- ggplot(df.plot, aes(x=annotation, y=0, fill=Class)) + 
  geom_tile() + 
  coord_flip() + 
  theme_minimal() + 
  theme(legend.position="bottom") +
  labs("Class annotation")
p_anno_class

### Tax 4
p_anno_tax4 <- ggplot(df.plot, aes(x=annotation, y=0, fill=TaxonomyRank4)) + 
  geom_tile() + 
  coord_flip() + 
  theme_minimal() + 
  theme(legend.position="none") + 
  labs("TaxonomyRank4")
p_anno_tax4

# ======================================================================= #
# ============================ patchwork =============================== #
# ======================================================================= #

p_anno_class + p + plot_layout(ncol = 2, widths = c(1,10)) 


p_anno_tax4 + p_anno_class + plot_layout(widths = c(5,5)) 

p_anno_tax4 + p_anno_class + p + plot_layout(ncol = 3, widths = c(1,1,10)) 

# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #

df.plot.tax_order <- df.ldsc_cts %>% mutate(annotation = factor(annotation, levels=annotation[order(TaxonomyRank4)]))

df.s <- df.plot.tax_order %>% count(TaxonomyRank4)
df.s <- df.s %>% 
  arrange(TaxonomyRank4) %>%
  mutate(pos_mid=cumsum(n)-n/2,
         pos_start=cumsum(n)-n,
         pos_end=cumsum(n),
         idx=1:n())
df.s <- df.s %>% mutate(flag_draw_rect=if_else(idx %% 2 == 0, TRUE, FALSE))
df.s

### Tax 4
p_anno_tax4 <- ggplot(df.plot.tax_order, aes(x=annotation, y=0, fill=TaxonomyRank4)) + 
  geom_tile() + 
  # annotate("text") + 
  geom_text(data=df.s, inherit.aes=F, aes(x=pos_mid, y=1, label=TaxonomyRank4)) +
  geom_rect(data=df.s %>% filter(flag_draw_rect), inherit.aes=F, aes(xmin=pos_start, xmax=pos_end, ymin=0, ymax=Inf), color="gray", alpha=0.3) +
  coord_flip() + 
  theme_minimal() + 
  theme(legend.position="none") + 
  labs("TaxonomyRank4")
p_anno_tax4

# inherit.aes=FALSE



# ======================================================================= #
# ============================ tmp =============================== #
# ======================================================================= #

df.tax4_cell_types <- df.ldsc_cts %>% 
  group_by(TaxonomyRank4) %>%
  summarise(celltypes=paste(annotation, collapse=","))
df.tax4_cell_types

df.tax4_cell_types.join <- df.tmp %>% left_join(df.tax4_cell_types)

df.tmp <- read_csv("TaxonomyRank4
Glutamatergic neuroblasts
Non-glutamatergic neuroblasts
Olfactory inhibitory neurons
Telencephalon inhibitory interneurons
Cholinergic and monoaminergic neurons
Peptidergic neurons
Spinal cord excitatory neurons
Spinal cord inhibitory neurons
Di- and mesencephalon excitatory neurons
Di- and mesencephalon inhibitory neurons
Hindbrain neurons
Cerebellum neurons
Telencephalon projecting excitatory neurons
Dentate gyrus granule neurons
Telencephalon projecting inhibitory neurons
Enteric neurons
Sympathetic noradrenergic neurons
Sympathetic cholinergic neurons
Peripheral sensory peptidergic neurons
Peripheral sensory neurofilament neurons
Peripheral sensory non-peptidergic neurons
Oligodendrocytes
Choroid epithelial cells
Subcommissural organ hypendymal cells
Ependymal cells
Dentate gyrus radial glia-like cells
Subventricular zone radial glia-like cells
Astrocytes
Olfactory ensheathing cells
Oligodendrocyte precursor cells
Schwann cells
Satellite glia
Enteric glia
Vascular and leptomeningeal cells
Vascular smooth muscle cells
Pericytes
Vascular endothelial cells
Perivascular macrophages
Microglia")
# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #


# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #
