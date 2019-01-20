

### Install
# source("https://bioconductor.org/biocLite.R")
# biocLite("ComplexHeatmap")
### REFs:
# http://www.bioconductor.org/packages/3.6/bioc/html/ComplexHeatmap.html (release 3.6 for R 3.4)
# https://github.com/jokergoo/ComplexHeatmap
# https://jokergoo.github.io/ComplexHeatmap-reference/book/

library(tidyverse)
library(ComplexHeatmap)

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/"
setwd(wd)

load(file="mousebrain.sem_obj.RData") # human

# ======================= meta SEM - per anno ======================= #

### REF 'Error: node stack overflow': https://github.com/renozao/NMF/issues/10#issuecomment-42807994
# the base heatmap function also breaks on this dataset, actually due to plot.dendrogram

for (anno_name in names(sem_obj[["group_by_annotation.sem_bin"]])) {
  print(anno_name)
  df.plot <- sem_obj[["group_by_annotation.sem_bin"]][[anno_name]] %>% 
    mutate(gene=sem_obj[["genes"]]) %>% 
    # slice(1:10) %>%
    column_to_rownames(var="gene") %>% as.data.frame() # this step must come last to keep the rownames.
  col_fun <- circlize::colorRamp2(c(0, 100), c("white", "red"))
  
  dummy <- tryCatch({ 
    message(sprintf("running ID=%s", anno_name))
    h1 <- Heatmap(df.plot, 
                  cluster_rows = T, 
                  cluster_columns = FALSE, 
                  show_row_names = F,
                  col = col_fun,
                  show_row_dend = F, # do not show row dendrogram | this might solve 'Error: node stack overflow'
                  use_raster = TRUE, 
                  raster_device = "png")
    file.out <- sprintf("heatmap_complex.sem_bin.%s.pdf", anno_name)
    pdf(file.out, width = 8, height = 20) # {ncols=5}{w=8}
    time_elapsed <- system.time(draw(h1))
    print(time_elapsed)
    dev.off()
  }, warning = function(w) {
    message(sprintf("warning: %s", w))
  }, error = function(e) {
    message(sprintf("************* error *************: %s", e))
  }, finally = {
    message(sprintf("done with ID=%s", anno_name))
  }) # end tryCatch
}

