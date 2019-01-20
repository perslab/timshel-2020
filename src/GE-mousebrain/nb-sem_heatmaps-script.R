
library(tidyverse)
# library(superheat) # https://rlbarter.github.io/superheat/basic-usage.html
# library(pheatmap) # https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
# library(d3heatmap) # interactive, https://cran.r-project.org/web/packages/d3heatmap/vignettes/Introduction.html
library(heatmaply) # https://github.com/talgalili/heatmaply | Install this lib before you can save figs --> webshot::install_phantomjs()

# wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/"
# setwd(wd)

load(file="mousebrain.sem_obj.RData") # human

for (sem_meta_name in names(sem_obj[["sem_meta"]])) {
  print(sem_meta_name)
  df.plot <- sem_obj[["sem_meta"]][[sem_meta_name]] %>% 
    mutate(gene=sem_obj[["genes"]]) %>% 
    # slice(1:10) %>%
    column_to_rownames(var="gene") %>% as.data.frame() # this step must come last to keep the rownames.
  file.out <- sprintf("heatmap.sem_meta.%s.html", sem_meta_name)
  heatmaply(df.plot,
            main=sem_meta_name,
            showticklabels = c(FALSE, FALSE),
            file=file.out)
}

print("Script done")