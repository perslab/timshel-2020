############### SYNOPSIS ###################
### AIM: generate ES summary stats figures

### OUTPUT: 
# ....

### REMARKS:
# Most content is copied and modified from "src/expression_specificity/es_obj_summary_stats.R"

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
# =========================== FUNCTIONS ========================= #
# ======================================================================= #

# wrapper_get_n_es_genes_summary_stats <- function()


# ======================================================================= #
# ================================ PARAMS ============================== #
# ======================================================================= #

# dataset_prefix <- "tabula_muris"
# var_color_by <- sym("tissue")

# dataset_prefix <- "mousebrain"
# var_color_by <- sym("Class")

dataset_prefix <- "hypothalamus"
var_color_by <- sym("taxonomy_lvl1")


# ======================================================================= #
# ================================ LOAD DATA ============================ #
# ======================================================================= #

### Read ES summary - old way [does not work for hypothalamus because file does not exist]
# file.data <- sprintf("tables/table-annotation_summary.%s.csv", dataset_prefix)
# df.n_es <- read_csv(file.data)
# df.metadata <- get_metadata(dataset_prefix)
# df.violin <- df.n_es %>% left_join(df.metadata, by="annotation")

### Read ES summary
file.data <- sprintf("tables/table-metadata_combined.%s.csv", dataset_prefix)
df.n_es <- read_csv(file.data)

if (dataset_prefix == "hypothalamus") {
  df.n_es <- df.n_es %>% mutate(annotation=annotation_fmt)
}

### Make dfs for plotting
df.violin <- df.n_es
df.barplot <- df.n_es

# ======================================================================= #
# ============================ VIOLIN PLOT ============================== #
# ======================================================================= #

### order *var_color_by* by var_color_by (mutate it to make it an ordered factor)
df.violin <- df.violin %>% mutate(!!var_color_by:=factor(!!var_color_by, levels=sort(unique(!!var_color_by)))) # this is not really needed now, but we use it later

### Define color mapping [after converting var_color_by to a factor]
factor_levels <- df.violin %>% pull(!!var_color_by) %>% levels() # levels are ordered in the same order as annotations are ordered
colormap <- colorspace::qualitative_hcl(n=length(factor_levels), palette = "Dark 3")
names(colormap) <- reorder_factor_levels_to_first_last_pairs(factor_levels)

### Plot
p.violin <- ggplot(df.violin, aes(x=!!var_color_by, y=ESmu, fill=!!var_color_by)) +
  geom_violin() + 
  geom_jitter(height=0, width=0.1) +
  geom_hline(yintercept=mean(df.violin$ESmu), color="blue", linetype="dashed") +
  labs(x="", y="Number of ES genes") + 
  scale_fill_manual(values=colormap) + 
  guides(fill=F) + 
  theme_classic() + # 'base theme'
  theme(axis.text.x= element_text(angle=25, hjust=1, size=rel(1))) # smaller x-axis labels
# p.violin <- p.violin + theme(plot.margin = unit(c(1,1,3,1), "cm")) # (t, r, b, l) widen margin
p.violin

file.out <- sprintf("figs/fig_n_es_genes.violin.%s.pdf", dataset_prefix)
ggsave(plot=p.violin, filename=file.out, width=10, height=5)


# ======================================================================= #
# =============================== BAR PLOT ================================== #
# ======================================================================= #


### order annotations by var_color_by
df.barplot <- df.barplot %>% mutate(!!var_color_by:=factor(!!var_color_by, levels=sort(unique(!!var_color_by)))) # this is not really needed now, but we use it later
df.barplot <- df.barplot %>% 
  arrange(!!var_color_by) %>% # arrange() respect factor orders
  mutate(annotation=factor(annotation, levels=annotation))

### Add var_color_by_fct_reorder variable: var_color_by ordered to "first last pairs"
# ----> HERE THE LEGEND BECOMES IN THE SAME 'not nice' ORDER.
# factor_levels <- df.barplot %>% pull(!!var_color_by) %>% levels() # levels are ordered in the same order as annotations are ordered
# df.barplot <- df.barplot %>% mutate(var_color_by_fct_reorder=factor(!!var_color_by, levels=reorder_factor_levels_to_first_last_pairs(factor_levels)))
# df.barplot$var_color_by_fct_reorder %>% levels()


### Define color mapping [after converting var_color_by to a factor]
factor_levels <- df.barplot %>% pull(!!var_color_by) %>% levels() # levels are ordered in the same order as annotations are ordered
colormap <- colorspace::qualitative_hcl(n=length(factor_levels), palette = "Dark 3")
names(colormap) <- reorder_factor_levels_to_first_last_pairs(factor_levels)


### Plot
p.bar <- ggplot(df.barplot, aes(x=annotation, y=ESmu, fill=!!var_color_by)) +
  geom_col() + 
  geom_hline(yintercept=mean(df.barplot$ESmu), color="blue", linetype="dashed") +
  geom_text(aes(label=annotation), angle=90, hjust=0, size=rel(1.5)) +
  labs(x="Cell-type", y="Number of ES genes") + 
  scale_fill_manual(values=colormap) + 
  scale_y_continuous(expand=c(0, 0)) + 
  coord_cartesian(clip="off") + 
  theme_classic() + # 'base theme'
  theme(axis.ticks.x = element_blank(), # wipe away all x-axis elements
        axis.text.x = element_blank(),
        axis.line.x = element_blank())
  # theme(axis.text.x= element_text(angle=75, hjust=1, size=rel(0.15))) # smaller x-axis labels
p.bar <- p.bar + theme(plot.margin = unit(c(3,1,1,1), "cm")) # (t, r, b, l) widen margin
p.bar

file.out <- sprintf("figs/fig_n_es_genes.barplot.%s.pdf", dataset_prefix)
ggsave(plot=p.bar, filename=file.out, width=10, height=5)

# ======================================================================= #
# =============================== HIST PLOT ================================== #
# ======================================================================= #

n_mean <- mean(df.barplot$ESmu)
p.hist <- ggplot(df.barplot, aes(x=ESmu)) +
  geom_histogram(color="black", fill="darkgray", alpha=0.8) + 
  geom_vline(xintercept=n_mean, color="blue", linetype="dashed") + 
  labs(x="Number of ES genes", y="Count") +
  scale_y_continuous(expand=c(0, 0))

p.hist.standalone <- p.hist + annotate("text",  x=n_mean, y = Inf, label=sprintf("mean=%s",round(n_mean)), vjust=3, hjust=-0.1, color="blue")
p.hist.standalone
file.out <- sprintf("figs/fig_n_es_genes.histogram.%s.pdf", dataset_prefix)
ggsave(plot=p.hist.standalone, filename=file.out, width=10, height=5)


y_limits <- layer_scales(p.bar)$y$range$range # get y-axis of barplot | REF: https://stackoverflow.com/a/35372274/6639640
p.hist.mod <- p.hist +
  annotate("text",  x=n_mean, y=0, label=sprintf("mean=%s",round(n_mean)), vjust=-0.5, hjust=-0.2, size=3, color="blue", fontface="bold") +
  scale_x_continuous(expand=c(0, 0)) + 
  theme(axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.line.y=element_blank()) +
  coord_flip() + 
  xlim(y_limits)

p.hist.mod

# ======================================================================= #
# ======================= PATCHWORK: BAR+HIST =========================== #
# ======================================================================= #

### Make blank plot
# REASON: SEE my notes on "Patchwork 'outer' margin issue"
p.blank <- ggplot(tibble(x=seq(1,10), y=rep(1, 10))) + geom_blank()


### Patch
p.patch <- p.blank + 
  ((p.bar+theme(legend.position="left")) + p.hist.mod + plot_layout(nrow=1, widths=c(1, 0.2))) +
  p.blank + 
  plot_layout(ncol=1, heights=c(0.1, 1, 0.1))
p.patch
# & theme(plot.margin = unit(c(3,1,1,1), "cm")) # (t, r, b, l) widen margin --> does not work for pathwork

### Save
file.out <- sprintf("figs/fig_n_es_genes.combined.%s.pdf", dataset_prefix)
ggsave(plot=p.patch, filename=file.out, width=10, height=5)

### Here we order the factor levels of the "var_color_by" to be able to give them different colors.
# NB: we are not reordering annotation levels.



# ======================================================================= #
# =============== Empirical pval distribution [NOT IMPORTANT] =========== #
# ======================================================================= #

# ### Empirical pval distribution
# df.plot <- get_empirical_distribution(es_obj, es_name="ges", annotation=NULL) # ACBG (median=1)
# df.plot.gather <- df.plot %>% gather(key="distribution", value="es_value")
# df.plot.gather %>% ggplot(aes(x=es_value, fill=distribution)) + geom_density(alpha=0.2) + xlim(c(-10,10))
# df.plot.gather %>% ggplot(aes(x=es_value, fill=distribution)) + geom_density(alpha=0.2) + xlim(c(10,1000))
