############### SYNOPSIS ###################
# Correlate Yengo and Loh2018 S-LDSC results


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

library(ggpubr) 
library(ggrepel)

source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))


# ======================================================================= #
# ================== LOAD LDSC CTS RESULTS (multi GWAS) ================ #
# ======================================================================= #

### Read LDSC results
file.results <- here("results/prioritization_celltypes--mousebrain_all.multi_gwas.csv.gz")
df.ldsc_cts <- read_csv(file.results)


# ======================================================================= #
# =============================== FACET PLOT: Locke+Yengo ================================= #
# ======================================================================= #


### Filter: GWAS
filter.gwas <- c("BMI_UKBB_Loh2018", "BMI_Locke2015", "BMI_UPDATE_Yengo2018")

### Filter: annotations
filter.annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")
colormap.annotations <- get_color_mapping.prioritized_annotations_bmi(dataset="mousebrain")

### Prep
df.plot <- df.ldsc_cts %>% 
  filter(gwas %in% filter.gwas) %>%
  mutate(p.value=-log10(p.value)) %>%
  select(gwas, p.value, annotation) %>%
  spread(key="gwas", value="p.value") %>%
  mutate(flag_priori = if_else(annotation %in% filter.annotations, TRUE, FALSE))
df.plot

### Gather non-Loh GWAS
cols.gather <- filter.gwas[!filter.gwas %in% "BMI_UKBB_Loh2018"]
df.plot <- df.plot %>% gather(key="gwas", value="p.value", (!!cols.gather))

### Rename GWAS
tmp_gwas_vector <- df.plot$gwas # convenience selector. If you want a specific order of the result, this vector should contain (unique) ordered values
newnames <- utils.rename_gwas(tmp_gwas_vector, style="fullname_author_year") # "fullname_author_year","fullname","abrv_author_year","abrv_year","abrv" 
rename_vector <- newnames; names(rename_vector) <- tmp_gwas_vector
df.plot <- df.plot %>% mutate(gwas_fmt = recode_factor(gwas, !!!rename_vector)) # Use a named character vector to recode factors with unquote splicing. | REF: https://dplyr.tidyverse.org/reference/recode.html

### Add GIANT to name
df.plot <- df.plot %>% mutate(gwas_fmt = paste0("GIANT consortium ", gwas_fmt))

p <- ggplot(df.plot, aes(x=BMI_UKBB_Loh2018, y=p.value)) + 
  geom_point(color="grey") + 
  geom_point(data=df.plot %>% filter(flag_priori), aes(color=annotation)) +
  geom_text_repel(data=df.plot %>% filter(flag_priori), aes(label=annotation)) +
  geom_abline() + 
  scale_color_manual(values=colormap.annotations) +
  labs(
    x=expression(atop(-log[10](P[S-LDSC]),"BMI (Loh, 2018)")),
    y=expression(-log[10](P[S-LDSC]))
  ) + 
  guides(color=F) +
  coord_fixed()
p <- p + ggpubr::stat_cor(method = "pearson")
p <- p + facet_grid(~gwas_fmt) + theme(strip.background=element_blank(), strip.text=element_text(face="bold"))
p 

file.out <- "figs/fig_gwas.giant_vs_loh_bmi_results.mousebrain.pdf"
ggsave(plot=p, filename=file.out, width=10, height=6)



# ======================================================================= #
# =============================== PLOT: Yengo only ================================= #
# ======================================================================= #

### Filter: GWAS
filter.gwas <- c("BMI_UKBB_Loh2018", "BMI_UPDATE_Yengo2018")

### Filter: annotations
filter.annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")
colormap.annotations <- get_color_mapping.prioritized_annotations_bmi(dataset="mousebrain")

### Prep
df.plot <- df.ldsc_cts %>% 
  filter(gwas %in% filter.gwas) %>%
  mutate(p.value=-log10(p.value)) %>%
  select(gwas, p.value, annotation) %>%
  spread(key="gwas", value="p.value") %>%
  mutate(flag_priori = if_else(annotation %in% filter.annotations, TRUE, FALSE))
df.plot



p <- ggplot(df.plot, aes(x=BMI_UKBB_Loh2018, y=BMI_UPDATE_Yengo2018)) + 
  geom_point(color="grey") + 
  geom_point(data=df.plot %>% filter(flag_priori), aes(color=annotation)) +
  geom_text_repel(data=df.plot %>% filter(flag_priori), aes(label=annotation)) +
  geom_abline() + 
  scale_color_manual(values=colormap.annotations) +
  labs(
    x=expression(atop(-log[10](P[S-LDSC]),"BMI (Loh, 2018)")),
    y=expression(atop("BMI (Yengo, 2018)", -log[10](P[S-LDSC])))
  ) + 
  guides(color=F) +
  coord_fixed()
p <- p + ggpubr::stat_cor(method = "pearson")
p

file.out <- "figs/fig_gwas.yengo_vs_loh.mousebrain.pdf"
ggsave(plot=p, filename=file.out, width=8, height=6)


