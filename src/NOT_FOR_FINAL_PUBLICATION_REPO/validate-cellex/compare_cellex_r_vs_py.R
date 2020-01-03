############### SYNOPSIS ###################
# Compare CELLEX-py vs R implementation


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

source(here("src/lib/load_functions.R")) # load sc-genetics library

# ======================================================================= #
# =============================== LOAD DATA ================================= #
# ======================================================================= #

### PT CELLEX results (CELLEX v1.0, first Pip release)
### Tobias: /nfsdata/projects/tobias/cellex/tutorials/tutorial_mousebrain_out/mousebrain.esmu.csv.gz
file.py.mouse <- here("tmp-data/cellex/mousebrain_cellex.mouse.esmu.csv.gz")
df.py.mouse <- read_csv(file.py.mouse)
df.py.mouse <- df.py.mouse %>% rename(gene=X1)

### R implementation
file.r.human <- here("out/es/mousebrain.mu.csv.gz")
df.r.human <- read_csv(file.r.human)

# ======================================================================= #
# ======================== Map CELLEX to orthologs ====================== #
# ======================================================================= #

### map mouse to human orthologs
df.py.mouse.tmp <- df.py.mouse %>% as.data.frame() %>% column_to_rownames("gene")
df.py.human <- mouse_to_human_ortholog_gene_expression_mapping(df.py.mouse.tmp, type_mouse_gene_ids="ensembl")
# [1] "Reading gene mapping files..."
# [1] "Subsetting on orthologs..."
# [1] "Number of genes not mapped to human ortholog: 5260 out of 20331 genes (25.87 pct)"
df.py.human <- df.py.human %>% rownames_to_column(var="gene") %>% as_tibble()



# ======================================================================= #
# ======================== Compute correlation =========================== #
# ======================================================================= #

stopifnot(all(colnames(df.py.human)==colnames(df.r.human))) # --> TRUE, cols match perfectly
mat.cor <- cor(df.py.human %>% select(-gene), df.r.human %>% select(-gene))
summary(diag(mat.cor))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9841  0.9944  0.9964  0.9958  0.9982  0.9999 

# ======================================================================= #
# =============================== Compare ================================= #
# ======================================================================= #

get_df_for_cmp <- function(df.py, df.r, annotation) {
  df <- full_join(df.py %>% select(gene, py=!!sym(annotation)), 
                  df.r %>% select(gene, r=!!sym(annotation)), 
                  by="gene")
  c <- with(df.plot, cor(py, r))
  p <- ggplot(df, aes(x=py, y=r)) + geom_point() + labs(title=sprintf("%s | Pearson r=%s", annotation, round(c,3)), x="CELLEX-Py", y="CELLEX-R")
  print(p)
  return(df)
}

df.plot <- get_df_for_cmp(df.py.human, df.r.human, "ABC")
df.plot <- get_df_for_cmp(df.py.human, df.r.human, "MEINH2")
df.plot <- get_df_for_cmp(df.py.human, df.r.human, "MEINH3")




# ======================================================================= #
# =============================== XXXXX ================================= #
# ======================================================================= #


# ======================================================================= #
# =============================== XXXXX ================================= #
# ======================================================================= #

