
library(tidyverse)
library(loomR)

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/wgcna_modules/"
setwd(wd)



# ======================= READ WGCNA cell embedding (ME per cell) ======================= #
file.in <- "/scratch/tmp-wgcna/maca_mods_scaled_kME_1_cellModEmbed.loom" 
lfile <- connect(file.in)

### Explore
names(lfile)
lfile[["col_attrs"]]
lfile[["col_attrs/cell_names"]][] # cell ids
lfile[["col_attrs/ident"]][] # tissue cell-type
lfile[["row_attrs"]]

### Get all cell-types
cell_type.names <- lfile[["col_attrs/ident"]][] # tissue cell-type
### Access all module names
module.names <- lfile[["row_attrs/module"]][] # same as lfile[["row_attrs/gene_names"]][]

### Make data frame
df <- as.data.frame(lfile[["matrix"]][,]) # OBS: cells x modules
colnames(df) <- module.names


# ======================= Calculate correlations amoung MEs ======================= #
# ensure the results are repeatable
set.seed(1)
library(caret)


# load(file="tmp.correlationMatrix.RData") # correlationMatrix 
# correlationMatrix <- cor(df) # calculate correlation matrix
# save(correlationMatrix, file="tmp.correlationMatrix.RData")
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.9, verbose=T) # 0.9 cutoff gives 276 modules that will be removed
print(highlyCorrelated) # print indexes of highly correlated attributes

colnames(df)[highlyCorrelated]

### Some output:
# Combination row 2088 and column 3020 is above the cut-off, value = 0.921 
# Flagging column 3020 
##  colnames(df)[c(2088,3020)] ---> 
# [1] "maca_mods_scaled_kME_1.tissue_cell_type_kME.Mammary_luminal epithelial cell of mammary gland.black_2072"
# [2] "maca_mods_scaled_kME_1.tissue_cell_type_kME.Trachea_endothelial cell.mediumseagreen_594

### find overlap
file.clusters <- "/raid5/projects/jonatan/tmp-maca/tables/maca_tissue_cell_type_kME_cell_cluster_module_genes.csv"
df.clusters <- read_csv(file.clusters)
df.clusters.sub <- df.clusters %>% filter(module %in% c("black_2072", "mediumseagreen_594"))
res <- df.clusters.sub %>% group_by(module) %>% do(var=.$ensembl)
intersect(res$var[[1]], res$var[[2]])

### Hypergeometric test
n_totalpop <- 15000 # total number of genes
n_overlap <- length(intersect(res$var[[1]], res$var[[2]]))
n_genesA <- length(res$var[[1]])
n_genesB <- length(res$var[[2]])
phyper(q=n_overlap, m=n_genesA, n=n_totalpop-n_genesA, k=n_genesB, lower.tail=FALSE) # https://stats.stackexchange.com/a/16259
print(sprintf("%s, %s, %s", n_overlap, n_genesA, n_genesB))

### We define genesA as the set of white balls. Then we 'draw a sample' genesB and look at the number of 'white balls' in the 'sample'.
# x/q = no. white balls in sample --> overlap
# m = no. white balls in urn ---> n_genesA (white balls)
# n = no. black balls in urn ---> n_totalpop-n_genesA (non-white balls ==> black balls)
# k = no. balls sampled (without replacement) ---> n_genesB




