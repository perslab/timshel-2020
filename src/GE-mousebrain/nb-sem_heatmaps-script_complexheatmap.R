

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

# ======================= meta SEM ======================= #

for (sem_meta_name in names(sem_obj[["sem_meta"]])) {
  # sem_meta_name <- "mean"
  print(sem_meta_name)
  df.plot <- sem_obj[["sem_meta"]][[sem_meta_name]] %>% 
    mutate(gene=sem_obj[["genes"]]) %>% 
    # slice(1:10) %>%
    column_to_rownames(var="gene") %>% as.data.frame() # this step must come last to keep the rownames.
  col_fun <- circlize::colorRamp2(c(0, 100), c("white", "red"))
  h1 <- Heatmap(df.plot, 
          cluster_rows = T, 
          cluster_columns = FALSE, 
          show_row_names = F,
          col = col_fun,
          use_raster = TRUE, 
          raster_device = "png")
  file.out <- sprintf("heatmap_complex.sem_meta.%s.pdf", sem_meta_name)
  pdf(file.out, width = 40, height = 20)
  time_elapsed <- system.time(draw(h1))
  print(time_elapsed)
  dev.off()
}

# # ======================= Jon code ======================= #
# ht1 <- Heatmap(cellModEmbed_mat[,], 
#                cluster_rows = F,
#                #row_order = ,
#                cluster_columns = T, 
#                #show_row_dend = F, 
#                show_column_dend = F, 
#                show_heatmap_legend = T, 
#                show_row_names = F, 
#                show_column_names = T,
#                heatmap_legend_param = list(title = "Expression"), use_raster = T)
# pdf(sprintf("%s%s_%s_%s_cellModEmbed_cellOrder.pdf", dir_plots, data_prefix, run_prefix, fuzzyModMembership),h=max(10, nrow(cellModEmbed_mat) %/% 1000), w=max(8, ncol(cellModEmbed_mat) %/% 2))
