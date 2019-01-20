

library(tidyverse)


wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/wgcna_modules/"
setwd(wd)



# ======================= MACA  ======================= #


### Read LDSC results
file.ldsc_cts <- "/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/wgcna.maca.BMI_Yengo2018.cell_type_results.txt"
df.ldsc_cts <- read_tsv(file.ldsc_cts) %>% rename(BETA=Coefficient, SE=Coefficient_std_error, P=Coefficient_P_value)
df.ldsc_cts <- df.ldsc_cts %>% mutate(module_id=stringr::str_split_fixed(Name,pattern="\\.",n=2)[,2]) # add module ID from 'Name' string, e.g. maca_tissue_cell_type.slateblue4
n.top_modules <- 100
top_modules <- df.ldsc_cts %>% top_n(n.top_modules, -P) %>% pull(module_id)


### Read cell embedding (takes some minutes)
file.cell_activity <- "/projects/jonatan/tmp-maca/tables/maca_kIM_embed_1_cellmodEmbed.csv"
df.cem <- data.table::fread(file=file.cell_activity, nThread=24, data.table=F, showProgress=T)
## V1                    V2 aliceblue_2 antiquewhite_3 aquamarine_4
## 1 Bladder_mesenchymal cell A12.D041914.3_8_M.1.1  -3.4248447      -9.890296    1.3593829
## 2     Bladder_bladder cell B16.D041914.3_8_M.1.1   7.2769709     -21.595570   -2.2994921
## 3     Bladder_bladder cell C18.D041914.3_8_M.1.1  19.9109020     -15.226803   -1.8486564
df.cem <- df.cem %>% rename(cell_type=V1, cell_id=V2)


# ======================= MOUSE BRAIN ======================= #


### Read LDSC results
file.ldsc_cts <- "/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/wgcna.mousebrain.BMI_Yengo2018.cell_type_results.txt"
df.ldsc_cts <- read_tsv(file.ldsc_cts) %>% rename(BETA=Coefficient, SE=Coefficient_std_error, P=Coefficient_P_value)
df.ldsc_cts <- df.ldsc_cts %>% mutate(module_id=stringr::str_split_fixed(Name,pattern="\\.",n=2)[,2]) # add module ID from 'Name' string, e.g. maca_tissue_cell_type.slateblue4
n.top_modules <- 100
top_modules <- df.ldsc_cts %>% top_n(n.top_modules, -P) %>% pull(module_id)

### Read cell embedding (takes some minutes)
file.cell_activity <- "/projects/jonatan/tmp-mousebrain/tables/mousebrain_Neurons_ClusterName_1b_kIM_cellModEmbed.csv"
df.cem <- data.table::fread(file=file.cell_activity, nThread=24, data.table=F, showProgress=T)
## ClusterName              orig.ident deepskyblue2 indianred3 mediumorchid3
## 1        ENT9 10X82_2_TCTCTCACCAGTTA-    -34.98637 -2.1115722     0.3407592
## 2        ENT9 10X82_2_TATTATCTACCAGA-    -42.80877 -1.9136726     7.4475931
df.cem <- df.cem %>% rename(cell_type=ClusterName, cell_id=orig.ident)

# ======================= CALC ACTIVITY ======================= #

### Extract specific columns
cols_top_modules <- colnames(df.cem) %in% top_modules
df.cem.top <- df.cem %>% select(cell_type, which(cols_top_modules)) # which(...) is needed to convert logical to numeric. REF: https://stackoverflow.com/a/26877820/6639640

#### MEAN SUMMARISE
df.summary <- df.cem.top %>% group_by(cell_type) %>% summarise_all(mean, na.rm=T)


# Pre-calculation [not finished/thought through]
# list.cell_type <- for each cell_type: # try %dopar% on outer loop [but note that this will copy the expression data matrix n_par times in memory!]
#   # n_c = n(cell_type)
#   # n_o = n(other)
#   df.gene = data.frame(n_genes, 4)
# for each gene: # inner loop should be the one with most variables
#   df.gene[gene, var_c] = var(x[cell_type,gene])
# df.gene[gene, mean_c] = mean(x[cell_type,gene])
# df.gene[gene, var_o] = var(x[other,gene])
# df.gene[gene, mean_o] = mean(x[other,gene])
# return(df.gene)
# df = bind_rows(list_cell_type) # will give ~2M rows (20k genes x 100 cell-types


library(doParallel)
library(foreach) # lib is also loaded as part of doParallel, but we just make sure to load it. Needed for "%dopar%"
registerDoParallel(30)
time.start <- proc.time()
iter_control <- df.cem.top %>% select(-cell_type) %>% colnames()
# iter_control <- iter_control[1:10]
list.res <- foreach (idx.foreach=1:length(iter_control), # columns
                     .inorder=TRUE, # .inorder=FALSE can increase performance, but we want things to be in order so we can set the names
                     .packages=c("tidyverse")) %dopar% {
  #list.par_analysis <- foreach (i=foreach_min:foreach_max, .packages=c("plyr")) %dopar% {
  time.loop.start <- proc.time()
  name_element.foreach <- iter_control[idx.foreach] # an element in iter_control, e.g. "slateblue4"
  print(paste("processing #", idx.foreach, "/#", length(iter_control), sep="" ))
  print(paste("name_element.foreach is:", name_element.foreach))
  df.res.t <- data.frame(cell_type=unique(df.cem.top$cell_type), t_test=NA, pval=NA, df=NA)
  for (cell_type_loop in unique(df.cem.top$cell_type)) {
    print(cell_type_loop)
    x <- df.cem.top %>% filter(cell_type==cell_type_loop) %>% select(UQ(rlang::sym(name_element.foreach)))
    y <- df.cem.top %>% filter(cell_type!=cell_type_loop) %>% select(UQ(rlang::sym(name_element.foreach)))
    t <- t.test(x,y, alternative="greater",var.equal=T)
    df.res.t[df.res.t$cell_type==cell_type_loop, "t_test"] <- t$statistic
    df.res.t[df.res.t$cell_type==cell_type_loop, "pval"] <- t$p.value
    df.res.t[df.res.t$cell_type==cell_type_loop, "df"] <- t$parameter
  }
  df.res.t
  
  # df.x <- df.cem.top %>% 
  #   group_by(cell_type) %>% 
  #   select(UQ(rlang::sym(name_element.foreach))) %>% # select module
  #   summarise(var_c=var(UQ(rlang::sym(name_element.foreach))),
  #             mean_c=mean(UQ(rlang::sym(name_element.foreach)))
  #             )
  # ### "return value"
  
}
### Setting names of result list
# foreach(..., .inorder=TRUE, ...): logical flag indicating whether the .combine function requires the task results to be combined in the same order that they were submitted. If the order is not important, then it setting .inorder to FALSE can give improved performance. The default value is TRUE.
# ---> since the ".inorder" is true per default, the list is combined in the same order that the elements were submitted.
names(list.res) <- iter_control
### Binding
df.res.foreach <- bind_rows(list.res, .id="module_id")
df.res.foreach



def.res.foreach.mat <- df.res.foreach %>% 
  select(module_id, cell_type, t_test) %>% 
  spread(key=module_id, value=t_test) %>% 
  column_to_rownames(var="cell_type")

library(plotly)
Sys.setenv("plotly_username"="pascal.timshel")
Sys.setenv("plotly_api_key"="9IldOzzfWvMyzLmy7t2Y")
## Row- and column-wise clustering 
hr <- hclust(as.dist(1-cor(t(def.res.foreach.mat), method="pearson")), method="complete") # rows
hc <- hclust(as.dist(1-cor(def.res.foreach.mat, method="pearson")), method="complete")  # columns
df.mat <- def.res.foreach.mat[hr$order,hc$order] # Re-arrange based on order [DOES NOT MATER IF YOU MELT AFTERWARDS]
df.gather <- df.mat %>% rownames_to_column(var="cell_type") %>% gather(key="module_id", value="value", -cell_type)
df.gather <- df.gather %>% mutate( # ordering factors
  cell_type = factor(cell_type, levels=rownames(def.res.foreach.mat)[hr$order]),
  module_id = factor(module_id, levels=colnames(def.res.foreach.mat)[hc$order])
)

### 'log' transform
df.gather <- df.gather %>% mutate(value_log = sign(value)*log10(abs(value)+1)) 
# ^log10(..+1): to avoid t-stats very close to zero of 'exploding' when doing log-transformation.
# ^sign(): returns 1, 0, or -1 if the number is positive, zero, or negative
# > x <- 22.000001; sign(x)*log10(abs(x)) ---> 1.342423 | with + 1, no real effect
# > x <- 22.000001; sign(x)*log10(abs(x)+1) ---> 1.361728 | no +1, no real effect
# > x <- 0.000001; sign(x)*log10(abs(x)) ---> -6 | no +1, large result number, which it should not be
# > x <- 0.000001; sign(x)*log10(abs(x)+1) ---> 4.342943e-07 | with +1, small result number as it should be

# Plot
p <- ggplot(df.gather, aes(cell_type, module_id, fill = value_log)) + geom_tile() + scale_fill_gradient2() # limits=c(-100,100) does not work with ggplotly
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
pp <- ggplotly(p)
pp


### mousebrain
# pp.mousebrain.top100 <- pp
# df.res.foreach.mousebrain.top100 <- df.res.foreach
### maca
# pp.maca.top100 <- pp
# df.res.foreach.maca.top100 <- df.res.foreach


# api_create(pp, filename = "bmi-brain/module_activity.maca") # (geom_tile?) does not work for some reason

library(superheat) # https://rlbarter.github.io/superheat/basic-usage.html
superheat(def.res.foreach.mat, 
          # add clustering
          pretty.order.rows = T, pretty.order.cols = T, 
          row.dendrogram = T, col.dendrogram = T,
          bottom.label.text.angle = 90,
          # change the size of the labels
          left.label.size = 0.2,
          bottom.label.size = 1,
          force.left.label = TRUE, force.bottom.label = TRUE,
          scale = F) # a logical specifying whether or not to center and scale the columns of X



#### T-TEST SUMMARISE
# df.cem.join.top %>% select(1:4, -CellID) %>% group_by(ClusterName)





# ======================= OLD PLOTS ======================= #


df.summary.plot <- df.summary %>% filter(complete.cases(.)) %>% column_to_rownames(var="ClusterName") %>% as.data.frame() # ALT: drop_na()
head(df.summary.plot)
rownames(df.summary.plot)


library(superheat) # https://rlbarter.github.io/superheat/basic-usage.html
superheat(df.summary.plot, 
          # add clustering
          pretty.order.rows = T, pretty.order.cols = T, 
          row.dendrogram = T, col.dendrogram = T,
          bottom.label.text.angle = 90,
          # change the size of the labels
          left.label.size = 0.2,
          bottom.label.size = 1,
          force.left.label = TRUE, force.bottom.label = TRUE,
          scale = T) # a logical specifying whether or not to center and scale the columns of X



library(plotly)
# Do hierarchical clustering of cells (rows)
clust <- df.summary.plot %>% 
  dist() %>% 
  hclust()
names(clust)
ord <- clust$order # Get order
df <- df.summary.plot[ord,] # Re-arrange based on order
df <- df %>% rownames_to_column(var="ClusterName") %>% gather(key="module", value="value", -ClusterName)
# Plot
p <- ggplot(df, aes(ClusterName, module, fill = value)) + geom_tile() + scale_fill_gradient2()
p <- ggplotly(p)
p



# df.summary.top <- df.summary %>% gather(key="module", value="activity", -cell_type) %>% group_by(module) %>% arrange(activity) %>% slice(c(1:5,(n()-5):n()))
# ALT1: top_n(n=5, wt=activity) # --> can only give you top
# ALT2: filter(row_number()==1 | row_number()==n()) # --> should be able to work.



# ======================= ENRICHMENT TESTS ======================= #

activity_enrichment <- function(category.now, df, variable_to_test_for_enrichment) {
  # print(category.now)
  df.now <- df %>% mutate(category.tmp = if_else(UQ(rlang::sym(variable_to_test_for_enrichment))==category.now, category.now, "Other"))
  x <- df.now %>% filter(category.tmp==category.now) %>% pull(bt_value)
  y <- df.now %>% filter(category.tmp=="Other") %>% pull(bt_value)
  ttest <- t.test(x,y , alternative="greater") # alternative hypothesis: x > y.
  # wtest <- wilcox.test(bt_value ~ category.tmp, data = df.now, alternative="greater") # *OBS*: using formula is DANGEROUS (when using a one-sided test) because the levels in 'category.tmp' are ordered alphabetically
  return(ttest$statistic)
}


df.wtest_res <- data.frame(enrichment=sapply(df.dataset %>% 
                                               distinct(UQ(rlang::sym(variable_to_test_for_enrichment))) %>% 
                                               pull(), enrichment_category, df=df.dataset, variable_to_test_for_enrichment=variable_to_test_for_enrichment)) # SWITCH = Class
df.wtest_res <- df.wtest_res %>% rownames_to_column(var="category")
df.wtest_res
