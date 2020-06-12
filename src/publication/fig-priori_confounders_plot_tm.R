############### SYNOPSIS ###################
# Correlate Tabula Muris CELLECT celltype prioritization scores with confounders:
# n_cells in cluster, median UMI count, median gene count


### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# =============================== SETUP ================================= #
# ======================================================================= #

.libPaths(c(.libPaths(),"~/R/x86_64-pc-linux-gnu-library/"))
library(here)
library("vctrs", lib.loc = "~/R/x86_64-pc-linux-gnu-library/")
library("rlang", lib.loc = "~/R/x86_64-pc-linux-gnu-library/")
library("tidyverse")
library("pbkrtest", lib.loc = "/raid5/home/cbmr/wzx816/R/x86_64-redhat-linux-gnu-library/3.5")
library("ggpubr")
library("ggrepel")
library("patchwork")
library("data.table")

source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))


# ======================================================================= #
# ================== LOAD LDSC CTS RESULTS (multi GWAS) ================ #
# ======================================================================= #

# TODO: Move results to timshelcelltypes dir (concat)
file.prior_results <- "/scratch/rkm916/CELLECT/results/prioritization.csv"
#file.prior_results <- here("results","cellect_ldsc", "prioritization.csv")
# Read confounder results
file.confounder_stats <- here("out", "qc_checks", "tabula_muris_cell_type_confounderstats.csv")
# Read esmu 
file.esmu.tabula_muris <- here("out", "es", "tabula_muris.mu.csv.gz")

df_ldsc_results = read_csv(file.prior_results)
dt_esmu = fread(file.esmu.tabula_muris)
df_confounder_stats = read_csv(file.confounder_stats)

#df.metadata <- read_csv(file.metadata)

df_ldsc_results <- df_ldsc_results %>% filter(specificity_id == "tabula_muris")

# ======================================================================= #
# =============================== FACET PLOT: Locke+Yengo ================================= #
# ======================================================================= #


### Filter: GWAS
filter.gwas <- paste0("1KG_phase3_EUR_null_gwas.P",1:1000)

### Filter: annotations
#filter.annotations <- get_prioritized_annotations_bmi(dataset="tabula_muris")
#colormap.annotations <- get_color_mapping.prioritized_annotations_bmi(dataset="tabula_muris")

### Prep
df.plot <- df_ldsc_results %>% 
  filter(gwas %in% filter.gwas) %>%
  mutate(p.value=-log10(pvalue)) %>%
  select(gwas, p.value, annotation) %>%
  spread(key="gwas", value="p.value") #%>%
#mutate(flag_priori = if_else(annotation %in% filter.annotations, TRUE, FALSE)) 

df.plot
# 
# prefix invalid first character (number)
colnames(df.plot) = gsub("^1KG", "thousandG", colnames(df.plot))

# verify matching order
all.equal(df.plot$annotation, df_confounder_stats$tissue_celltype)
# TRUE 
all.equal(df.plot$annotation, colnames(dt_esmu)[2:ncol(dt_esmu)])
# [1] "54 string mismatches"
# correct mismatch
vec_idx = match(df.plot$annotation, colnames(dt_esmu))
dt_esmu = dt_esmu[,c(1,vec_idx), with=F]

all.equal(df.plot$annotation, colnames(dt_esmu)[2:ncol(dt_esmu)])
# [1] TRUE

vec_n.esmu = colSums(dt_esmu[,2:ncol(dt_esmu)]>0)

# compute correlation between no. cells and no. esmu genes
cor(vec_n.esmu, df_confounder_stats$NCells)
# [1] 0.3631044

# regress out no. esmu genes from no. cells
# purpose: show that the bias from no. cells on prioritization p-values arises through no. esmu genes
df_confounder_stats <- add_column(df_confounder_stats, 
                                  "n_esmu_genes"=vec_n.esmu)

regr = lm(formula = "NCells ~ n_esmu_genes", data = df_confounder_stats)
vec_residuals = summary(regr)$residuals#
df_confounder_stats <- add_column(df_confounder_stats, "n_cells_residuals" = vec_residuals)

# compute correlation between P-LDSC and residual n cells
mat_cor_gwas_ncells_residuals <- t(sapply(gsub("^1KG", "thousandG", filter.gwas), function(gwas) {
  testout = cor.test(df.plot[[gwas]], df_confounder_stats[["n_cells_residuals"]])
  return(c("estimate"=testout$estimate, "p-value"=testout$p.value))
}))

t.test(mat_cor_gwas_ncells_residuals[,"estimate.cor"])
# p-value = 0.08159

# check correlation between P-LDSC and non=zero esmu genes
mat_cor_gwas_ncells_n.esmu <- t(sapply(gsub("^1KG", "thousandG", filter.gwas), function(gwas) {
  testout = cor.test(df.plot[[gwas]], df_confounder_stats[["n_esmu_genes"]])
  return(c("estimate"=testout$estimate, "p-value"=testout$p.value))
}))

t.test(mat_cor_gwas_ncells_n.esmu)
# p-value = 0.002145

df.plot = add_column(.data = df.plot, 
                     "NCells"=df_confounder_stats$NCells, 
                     "nCount_RNA_median"=df_confounder_stats$nCount_RNA_median, 
                     "nCount_gene_median"=df_confounder_stats$nCount_gene_median,
                     "n_esmu_genes" = vec_n.esmu,
                     "n_cells_residuals" = df_confounder_stats$n_cells_residuals)

# compute correlation matrix
list_cor_gwas_confounder <- lapply(gsub("^1KG", "thousandG", filter.gwas), function(gwas) {
  mat_out = t(sapply(c("NCells", "nCount_RNA_median","nCount_gene_median"), function(confounder) {
    testout = cor.test(df.plot[[gwas]], df.plot[[confounder]])
    c("gwas" = gwas, "confounder"=confounder, "estimate"=testout$estimate, "p-value"=testout$p.value)
  }))
})

mat_cor_gwas_confounder = Reduce(x=list_cor_gwas_confounder, f=rbind) 
dt_cor_gwas_confounder <- data.table(mat_cor_gwas_confounder)

# file.out <- ("/projects/jonatan/pub-perslab/20-BMI-brain/output/null_gwas_confounder_analysis_tabula_muris.csv")
# fwrite(dt_cor_gwas_confounder, file.out)
# file.out.xlsx <- ("/projects/jonatan/pub-perslab/20-BMI-brain/output/null_gwas_confounder_analysis_tabula_muris.xlsx")
# openxlsx::write.xlsx(dt_cor_gwas_confounder, file.out.xlsx)


### PLOTS ###
# Description: histogram of CELLECT celltype LDSC p-values across gwas, (no -log10 transform)
# Purpose: show that p-values for null gwas are uniformly distributed

df_ldsc_results %>% 
  filter(gwas %in% filter.gwas) 

phisto_pvalue <- ggplot(df_ldsc_results, aes(x=pvalue)) + geom_histogram(breaks=seq(0, 1, by = 0.01), 
                                                            fill="grey") + 
  labs(x=expression(atop(P[S-LDSC]),"1KG phase3 EUR null GWAS"),
       y="count") +
  theme_classic() +
  theme(strip.background=element_blank(), strip.text=element_text(face="bold"))

# Description: P-LDSC vs confounder (NCells, nCount_RNA_medianfa, nCount_gene_median), correlation Rho
# Purpose: show that correlations between -log10(LDSC p-values) and confounders tend towards 0

dt_cor_gwas_confounder$confounder <- factor(dt_cor_gwas_confounder$confounder)
dt_cor_gwas_confounder$estimate.cor <- as.numeric(dt_cor_gwas_confounder$estimate.cor)

anno1 = paste0("*** P=",format(t.test(dt_cor_gwas_confounder[confounder=="NCells", estimate.cor])$p.value, digits=1))
anno2 = paste0("P=",format(t.test(dt_cor_gwas_confounder[confounder=="nCount_gene_median", estimate.cor])$p.value, digits=1))
anno3 = paste0("P=",format(t.test(dt_cor_gwas_confounder[confounder=="nCount_RNA_median", estimate.cor])$p.value, digits=1))

pbox_rho <- ggplot(data=dt_cor_gwas_confounder, aes(x=confounder,y=estimate.cor)) +
  geom_boxplot() + 
  scale_x_discrete(labels=c("no. cells", "median no. genes", "median no. UMI")) +
  labs(
    y=expression(atop("Pearson correlation with",-log[10](P[S-LDSC])))) +
  theme_classic() +
  theme(strip.background=element_blank(), 
        strip.text=element_text(face="bold"),
        axis.title.x=element_blank()) + 
  geom_signif(stat="signif", 
              test="t.test",
              #map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05))#, 
              annotations = c(anno1, anno2, anno3),
              y_position=c(0.6,0.6,0.6), 
              xmin=c(1,2,3), xmax = c(1,2,3))
  

# now plot the esmu correlation with P-LDSC
dt_cor_gwas_confounder2 = data.table("confounder"="n_esmu_genes", mat_cor_gwas_ncells_n.esmu)
dt_cor_gwas_confounder2 = rbind(dt_cor_gwas_confounder2, data.table("confounder"="n_cells_residuals", mat_cor_gwas_ncells_residuals))

anno4 = paste0("** P=",format(t.test(dt_cor_gwas_confounder2[confounder=="n_esmu_genes", estimate.cor])$p.value, digits=1))
anno5 = paste0("P=",format(t.test(dt_cor_gwas_confounder2[confounder=="n_cells_residuals", estimate.cor])$p.value, digits=1))

pbox_rho_n.esmu <- ggplot(data=dt_cor_gwas_confounder2, aes(x=confounder,y=estimate.cor)) +
  geom_boxplot() + 
  scale_x_discrete(labels=c("no. esmu genes", "no. cells (regression residuals)")) +
  labs(
    y=expression(atop("Pearson correlation with",-log[10](P[S-LDSC])))) +
  theme_classic() +
  theme(strip.background=element_blank(), 
        strip.text=element_text(face="bold"),
        axis.title.x=element_blank()) +
  geom_signif(stat="signif", 
              test="t.test",
              #map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05))#, 
              annotations = c(anno4, anno5),
              y_position=c(0.6,0.6), 
              xmin=c(1,2), xmax = c(1,2)) 
  

pout = phisto_pvalue + pbox_rho + pbox_rho_n.esmu

file.out <- here("src", "publication","figs","tabula_muris_nullgwas_confounders_panel.pdf")
ggsave(plot=pout, filename=file.out, width=14, height=6)




##### 
# LEFTOVERS
#####

if (F) {
  p1 <- ggplot(df.plot, aes(x=NCells,y=v1KG_phase3_EUR_null_gwas_P1)) +
    scale_x_continuous(trans = "log10") + 
    geom_point(color="grey") + 
    #geom_point(data=df.plot, aes(color=annotation)) +
    #geom_text_repel(data=df.plot %>% filter(flag_priori), aes(label=annotation), size=2) +
    #geom_abline() + 
    #scale_color_manual(values=colormap.annotations) +
    labs(
      x="n cells in cluster",
      y=expression(atop(-log[10](P[S-LDSC]),"1KG phase3 EUR null GWAS P1"))
      )
     #+ guides(color=F) +
     #coord_fixed()
  #p1 <- p1 + coord_fixed()
  p1 <- p1 + ggpubr::stat_cor(method = "pearson") # this is the correlation of the log x axis
  p1 <- p1 + theme_classic()
  p1 <- p1 + theme(strip.background=element_blank(), strip.text=element_text(face="bold"))
  
  #p <- p + facet_grid(~gwas_fmt) + theme(strip.background=element_blank(), strip.text=element_text(face="bold"))
  p1
  
  ### median RNA count ####
  
  p2 <- ggplot(df.plot, aes(x=nCount_RNA_median,y=v1KG_phase3_EUR_null_gwas_P1)) +
    #scale_x_continuous(trans = "log10") + 
    geom_point(color="grey") + 
    #geom_point(data=df.plot, aes(color=annotation)) +
    #geom_text_repel(data=df.plot %>% filter(flag_priori), aes(label=annotation), size=2) +
    #geom_abline() + 
    #scale_color_manual(values=colormap.annotations) +
    labs(
      x="median RNA count",
      y=expression(atop(-log[10](P[S-LDSC]),"1KG phase3 EUR null GWAS P1"))
    )
  #+ guides(color=F) +
  #coord_fixed()
  #p2 <- p2 + coord_fixed()
  p2 <- p2 + ggpubr::stat_cor(method = "pearson") # this is the correlation of the log x axis
  p2 <- p2 + theme_classic()
  p2 <- p2 + theme(strip.background=element_blank(), strip.text=element_text(face="bold"))
  #p <- p + facet_grid(~gwas_fmt) + theme(strip.background=element_blank(), strip.text=element_text(face="bold"))
  p2
  
  
  ### median gene count ### 
  
  p3 <- ggplot(df.plot, aes(x=nCount_gene_median,y=v1KG_phase3_EUR_null_gwas_P1)) +
    #scale_x_continuous(trans = "log10") + 
    geom_point(color="grey") + 
    #geom_point(color="grey") + 
    #geom_point(data=df.plot, aes(color=annotation)) +
    #geom_text_repel(data=df.plot %>% filter(flag_priori), aes(label=annotation), size=2) +
    #geom_abline() + 
    #scale_color_manual(values=colormap.annotations) +
    labs(
      x="median number of genes detected",
      y=expression(atop(-log[10](P[S-LDSC]),"1KG phase3 EUR null GWAS P1"))
    )
  #+ guides(color=F) +
  #coord_fixed()
  #p3 <- p3 + coord_fixed()
  p3 <- p3 + ggpubr::stat_cor(method = "pearson") # this is the correlation of the log x axis
  p3 <- p3 + theme_classic()
  p3 <- p3 + theme(strip.background=element_blank(), strip.text=element_text(face="bold"))
  #p <- p + facet_grid(~gwas_fmt) + theme(strip.background=element_blank(), strip.text=element_text(face="bold"))
  p3
  
  
  p <- p1+p2+p3
  file.out <- here("src", "publication","figs","tabula_muris_tissue_celltype_confounderstats.pdf")
  ggsave(plot=p, filename=file.out, width=18, height=6)
  
  
  #### check correlation after removing lowest counts #### 
  
  df.plot <- arrange(df.plot,nCount_RNA_median)
  
  sapply(0:200, function(i) {
    cor.test_out <- cor.test(x = log10(df.plot$nCount_RNA_median[i:nrow(df.plot)]), 
             y = df.plot$filter.gwas[i:nrow(df.plot)]) 
    c("n_dropped" = i, "estimate"=cor.test_out$estimate, "p.value"=cor.test_out$p.value)
    
  }) %>% t %>% as.data.frame -> df_cortest_nCount_RNA
  
  
  p4 <- ggplot(df_cortest_nCount_RNA, aes(x=as.numeric(n_dropped),y=estimate.cor)) +
    geom_point(color="grey") 
  p4
  
  ### remove cell clusters, starting from high ###
  
  df.plot <- arrange(df.plot,desc(nCount_RNA_median))
  
  sapply(0:200, function(i) {
    cor.test_out <- cor.test(x = log10(df.plot$nCount_RNA_median[i:nrow(df.plot)]), 
                             y = df.plot$filter.gwas[i:nrow(df.plot)]) 
    c("n_dropped" = i, "estimate"=cor.test_out$estimate, "p.value"=cor.test_out$p.value)
    
  }) %>% t %>% as.data.frame -> df_cortest_nCount_RNA_2
  
  
  p5 <- ggplot(df_cortest_nCount_RNA_2, aes(x=as.numeric(n_dropped),y=estimate.cor)) +
    geom_point(color="grey") 
  p5
  
  
  model.matrix_TaxonomyRank2 = model.matrix(object = formula("~ TaxonomyRank2"), data = df.metadata)
  model.matrix_TaxonomyRank2 = model.matrix_TaxonomyRank2[match(df.plot$annotation,df.metadata$ClusterName),]
  
  lm1 = lm(df.plot$BMI_UKBB_Loh2018 ~ model.matrix_TaxonomyRank2 + log10(df.plot$nCount_RNA_median))
  
  summary(lm1)
  
  cor.test(model.matrix_TaxonomyRank2[,"TaxonomyRank2CNS neurons"], log10(df.plot$nCount_RNA_median))
  
  df.taxonomy = data.frame(model.matrix_TaxonomyRank2, select(.data = df.plot, nCount_RNA_median,BMI_UKBB_Loh2018))
  
  df.taxonomy %>%
    select(BMI_UKBB_Loh2018,NCells, nCount_RNA_median, TaxonomyRank2CNS.neurons, TaxonomyRank2Immune.cells, TaxonomyRank2Neural.crest.like.glia,TaxonomyRank2PNS.neurons, TaxonomyRank2Vascular.cells) %>%
    mutate(NnCount_RNA_median = log10(nCount_RNA_median), nCount_gene_median=log10(nCount_gene_median)) %>%
    plot
  # ### Compare all confoudners 
  # 
  # df.plot %>% 
  #   select(BMI_UKBB_Loh2018,NCells, nCount_RNA_median, nCount_gene_median) %>%
  #   mutate(NCells= log10(NCells), nCount_RNA_median = log10(nCount_RNA_median), nCount_gene_median=log10(nCount_gene_median)) %>%
  #   plot
}