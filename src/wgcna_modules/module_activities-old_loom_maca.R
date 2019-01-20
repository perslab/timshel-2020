

library(tidyverse)
library(loomR)

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/wgcna_modules/"
setwd(wd)


# ======================= MACA  ======================= #
file.in <- "/scratch/tmp-wgcna/maca_mods_scaled_kME_1_cellModEmbed.loom"  # ME per cell
# file.in <- "/projects/jonatan/tmp-maca/RObjects/maca_kIM_embed_1_cellModEmbed.loom" # kIM ---> MISSING MODULE ID
lfile <- connect(file.in)

### Make data frame
df.me <- as.tibble(lfile[["matrix"]][,]) # OBS: cells x modules. Takes some time ...

### Explore
names(lfile)
lfile[["col_attrs"]]
lfile[["col_attrs/cell_names"]][] # cell ids
lfile[["col_attrs/ident"]][] # tissue cell-type
lfile[["row_attrs"]]

### Access all module names
module.names <- lfile[["row_attrs/module"]][] # same as lfile[["row_attrs/gene_names"]][]
# tmp.match_matrix <- stringr::str_match(module.names, "^maca_kIM_embed_1.tissue_cell_type_kIM.(.*)_(.*)\\.(.*)$") # for kIM EMBEDDING --> MISSING MODULE ID # E.g. "maca_kIM_embed_1.tissue_cell_type_kIM.Lung_leukocyte32"
tmp.match_matrix <- stringr::str_match(module.names, "^maca_mods_scaled_kME_1.tissue_cell_type_kME.(.*)_(.*)\\.(.*)$") # for ME EMBEDDING # E.g. "maca_mods_scaled_kME_1.tissue_cell_type_kME.Fat_myeloid cell.mediumspringgreen_596"
tissue <- tmp.match_matrix[,2]
cell_type <- tmp.match_matrix[,3]
module_id <- tmp.match_matrix[,4]
module_name <- paste(tissue, cell_type, module_id, sep="-")
df.metadata.modules <- tibble(tissue=tissue, cell_type=cell_type, module_id=module_id, module_name=module_name)


### Get all cell-types
cell_type.names <- lfile[["col_attrs/ident"]][] # tissue cell-type

### Set metadata of df.me
colnames(df.me) <- df.metadata.modules$module_id
df.me <- df.me %>% mutate(cell_type = cell_type.names) %>% select(cell_type, everything())


### Summarise
modules_to_summarise <- c("slateblue4", "bisque3", "salmon", "blueviolet_269", "mediumorchid_2398") # MACA
df.summary <- df.me %>% group_by(cell_type) %>% summarise_at(.vars=modules_to_summarise, mean)

### RENAME modules: add origin
names(modules_to_summarise) <- df.metadata.modules %>% slice(match(modules_to_summarise, module_id)) %>% pull(module_name) # named vector for renaming
modules_to_summarise
df.summary <- df.summary %>% rename(UQS(modules_to_summarise)) # rename

### CODE TO EXTRACT TOP AND BOTTOM ACTIVE CELL-TYPES FOR EACH MODULE
df.summary.top <- df.summary %>% gather(key="module", value="activity", -cell_type) %>% group_by(module) %>% arrange(activity) %>% slice(c(1:5,(n()-5):n()))
# ALT1: top_n(n=5, wt=activity) # --> can only give you top
# ALT2: filter(row_number()==1 | row_number()==n()) # --> should be able to work.




# ======================= MOUSE BRAIN OLD (JOINING CELL_TYPE INFORMATION) ---> DELETE ======================= #

### Read LDSC results
file.ldsc_cts <- "/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/wgcna.mousebrain.BMI_Yengo2018.cell_type_results.txt"
df.ldsc_cts <- read_tsv(file.ldsc_cts) %>% rename(BETA=Coefficient, SE=Coefficient_std_error, P=Coefficient_P_value)
df.ldsc_cts <- df.ldsc_cts %>% mutate(module_id=stringr::str_split_fixed(Name,pattern="\\.",n=2)[,2]) # add module ID from 'Name' string, e.g. maca_tissue_cell_type.slateblue4
n.top_modules <- 100
top_modules <- df.ldsc_cts %>% top_n(n.top_modules, -P) %>% pull(module_id)

### Read cell metadata
file.mousebrain_cell_metadata <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain-L5.metadata_cluster_name.csv.gz"
df.mousebrain_cell_metadata <- read_csv(file.mousebrain_cell_metadata)
# sum(duplicated(df.mousebrain_cell_metadata$CellID)) # ---> 118
# sum(duplicated(df.mousebrain_cell_metadata)) # --> 16
df.mousebrain_cell_metadata <- df.mousebrain_cell_metadata %>% filter(!duplicated(CellID)) # remove any duplicates

### Read cell embedding (takes some minutes)
file.cell_activity <- "/projects/jonatan/tmp-mousebrain/tables/mousebrain_Neurons_ClusterName_1b_kIM_cellModEmbed.csv"
# df.cem <- read_csv(file.cell_activity)
# df.cem <- df.cem %>% rename(CellID = X1)
sum(!df.cem$CellID %in% df.mousebrain_cell_metadata$CellID) # ---> 41 cells missing in cell_metadata, that's ok.

### Left join: keep only observations with cell-type assignments from Mousebrain metadata
df.cem.join <- left_join(df.mousebrain_cell_metadata, df.cem, by="CellID")
head(df.cem.join)


### Extract specific columns
cols_top_modules <- (str_split_fixed(colnames(df.cem.join), "__", n=2)[,2] %in% top_modules)
df.cem.join.top <- df.cem.join %>% select(CellID, ClusterName, which(cols_top_modules)) # which(...) is needed to convert logical to numeric. REF: https://stackoverflow.com/a/26877820/6639640

### save data
# saveRDS(df.cem.join.top, file="module_embedding.mousebrain_neurons.top_modules.Rds")
# saveRDS(df.cem.join, file="module_embedding.mousebrain_neurons.all_modules.Rds")

### Write tables
# data.table::fwrite(df.cem.join.top, file="module_embedding.mousebrain_neurons.top_modules.csv")
# data.table::fwrite(df.cem.join, file="module_embedding.mousebrain_neurons.all_modules.csv", nThread=24)

### Read tables
system.time(df.cem.join.top <- data.table::fread(file="module_embedding.mousebrain_neurons.top_modules.csv", nThread=24, data.table=F, showProgress=T))
df.cem.join.top <- as.tibble(df.cem.join.top)

#### MEAN SUMMARISE
df.summary <- df.cem.join.top %>% group_by(ClusterName) %>% summarise_if(is.numeric, mean, na.rm=T)


