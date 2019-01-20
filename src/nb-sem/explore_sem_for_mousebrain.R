### INFO: dirty script to make plot of some key stats for the mousebrain SEM matrix files


library(tidyverse)
library(skimr) # REF: https://cran.r-project.org/web/packages/skimr/vignettes/Using_skimr.html


DIR.X <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/mousebrain"
filenames <- list.files(path=DIR.X,  pattern="*.hsapiens_orthologs.csv.gz")
list.x <- lapply(file.path(DIR.X, filenames), read_csv)
list.x <- lapply(list.x, function(df) df %>% select(-gene))
names(list.x) <- filenames
list.sum <- lapply(list.x, skim) # summarise
list.sum
list.sum.sum <- lapply(list.sum, function(df) {
  df %>% select(-formatted) %>% 
    tidyr::spread(key="stat", value="value") %>%
    select(variable,mean,sd,p100) %>%
    mutate(cv=sd/mean)
})
list.sum.sum[[1]]

df.sum.sum <- bind_rows(list.sum.sum, .id = "file") # combine
df.sum.sum

### Inspect
df.x <- df.sum.sum %>% group_by(file) %>%
  top_n(5, p100)
df.x

### Plot
ggplot(df.sum.sum, aes(x=cv, fill=file)) + geom_density(alpha=0.4) + xlim(0,15) # CV
ggplot(df.sum.sum, aes(x=file, y=log10(cv), fill=file)) + geom_boxplot() # CV
ggplot(df.sum.sum, aes(x=file, y=log(p100), fill=file)) + geom_boxplot() # p100




############## TMP ##################
df <- list.sum[[1]]
head(df)
df %>% 
  select(-formatted) %>% 
  tidyr::spread(key="stat", value="value") %>%
  select(variable,mean,sd,p100) %>%
  mutate(cv=sd/mean)
  # group_by(variable) %>%
  # summarise(mean_mean=mean())
?spread


skim(list.x[[1]])
skim(list.x[[2]])
skim(list.x[[3]])

