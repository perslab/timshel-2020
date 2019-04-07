
library(tidyverse)

### Load the celltype data
load(file="rsession_celltyping.mousebrain.RData")

### Meta data
file.metadata <- "/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain-agg_L5.metadata.csv"
df.metadata <- read_csv(file.metadata)
cols.keep.metadata <- c("Class",
                        "ClusterName",
                        "Description",
                        "Developmental_compartment",
                        "Neurotransmitter",
                        "TaxonomyRank1",
                        "TaxonomyRank2",
                        "TaxonomyRank3",
                        "TaxonomyRank4",
                        "TaxonomySymbol",
                        "NCells")
df.metadata <- df.metadata %>% select(cols.keep.metadata)
df.metadata <- df.metadata %>% mutate(cell_type_id=paste0(Class, "_",ClusterName)) %>% select(cell_type_id, everything())
df.metadata

### Join results + metadata
df <- df.results %>% left_join(df.metadata, by=c("COVAR"="cell_type_id"))
df <- df %>% arrange(desc(P)) %>% mutate(COVAR = factor(COVAR, unique(COVAR))) # Prepating for barplot: ORDERING name_clean by the 'category'
# The key to (re)ordering factor levels, is that it MUST be a unique set of values supplied to 'levels'

### Write CSV
df %>% arrange(P) %>% write_csv("out.magma_celltyping.results_with_metadata.mousebrain.csv")

### Plot
for (col in cols.keep.metadata) {
  print(col)
  p = ggplot(data=df)+ # ctAssocs[[annotLevel]]$results
    geom_bar(aes(y=P,x=COVAR, fill=get(col) ),stat="identity")+
    scale_y_log10()+
    coord_flip()+
    ylab(expression(-log[10](pvalue)))+
    xlab("Cell type") + 
    ggtitle(col)
  p <- p + geom_hline(yintercept=as.numeric(0.05/ctAssocs$total_baseline_tests_performed),colour="black")
  p
  ggsave(sprintf("out.tmp.mousbrain_color-%s.pdf", col), w=20, h=40)
}

p


df.metadata %>% filter(ClusterName=="HBGLU3") %>% transpose()

