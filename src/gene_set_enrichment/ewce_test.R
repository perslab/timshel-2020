### TMP TESTING EWCE FOR . 
### CONCLUSION (90% sure): EWCE uses get_summed_proportions() to calculate background and hence is NOT VALID for other ES input values.

### DOCS
# bootstrap.enrichment.test: ---> main function to calculate enrichment. Very simple code.
# get_summed_proportions: Given the target geneset, randomly sample genelists of equal length, obtain the specificity of these and then obtain the mean specificity in each sampled list (and the target list)
# 

# library(devtools)
# install_github("nathanskene/ewce")
library(EWCE)
library(tidyverse)

### Load
load("/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain.sem_obj.RData")

### Gene list
df.genes <- read_tsv("/raid5/projects/timshel/sc-genetics/sc-genetics/data/genes_obesity/yazdi2015_table1_mouse_obesity_genes.mapped.txt")
df.genes$ensembl_gene_id

genes.hits <- df.genes$ensembl_gene_id[!is.na(df.genes$ensembl_gene_id)]
genes.hits <- genes.hits[genes.hits %in% sem_obj$genes]

### Make CTD
ctd = list("dummy_annot_lvl1" = list(
  "mean_exp"=as.matrix(sem_obj$data$mean), # matrix, genes x cells. with rownames and columnnames
  "annot"=colnames(sem_obj$sem_meta$mean),
  "specificity"=as.matrix(sem_obj$sem_meta$mean)
))


### Params
reps <- 1000 # <- Use 1000 bootstrap lists so it runs quickly, for publishable analysis use >10000
level <- 1 # <- Use level 1 annotations (i.e. Interneurons)


# Bootstrap significance testing, without controlling for transcript length and GC content
full_results = bootstrap.enrichment.test(sct_data=ctd,
                                         hits=genes.hits,
                                         bg=sem_obj$genes,
                                         reps=reps,
                                         annotLevel=level)

print(full_results$results[order(full_results$results$p),3:5][1:6,])
print(ewce.plot(full_results$results,mtc_method="BH"))

f <- bootstrap.enrichment.test
f1 <- check.ewce.genelist.inputs

### Plotting results from multiple gene lists
# o achieve this the results data frames are just appended onto each other, with an additional list column added detailing which analysis they relate to.
full_res2 = data.frame(full_results$results,list="Alzheimers")
scnd_res2 = data.frame(second_results$results,list="Second")
merged_results = rbind(full_res2,scnd_res2)
ewce.plot(total_res=merged_results,mtc_method="BH")