# TMP DEBUGGING ROLYPLOY.
# OK DELETE


library(rolypoly)

if (!exists('is_rolypoly_pascaltimshel_fork', mode='function')) {
  # 'is_rolypoly_pascaltimshel_fork' function only exists in pascaltimshel fork.
  # exists() returns TRUE only if function is found.
  stop("You have not loaded the pascaltimshel forked version of rolypoly. Univariate functionality is only supported in pascaltimshel forked version. Will exit...")
}


### load GWAS data
load("/scratch/rolypoly_out_scratch/out.rolypoly_objs-v3.univariate/rolypoly_objs.BMI_UPDATE_Yengo2018_no_mhc.tss.100kb.none.nboot500.gwas_linked.RData")

### load expressio data
df.expr <- read_csv("/nfsdata/projects/timshel/sc-genetics/sc-genetics/data/genes_cell_type_specific/mousebrain_all.mean.csv.gz")
df.expr <- df.expr %>% 
  as.data.frame() %>%
  column_to_rownames(var = "gene")
name.annotation <- "ACNT2"
df.annotation <- df.expr %>% select(name.annotation) %>% as.data.frame()
head(df.annotation)

rp.precomputed <- list.gwas_linked[["rp"]] # --> object
block_data <- df.annotation

tmp.time <- system.time(rolypoly <- rolypoly:::rolypoly_load_block_data(rolypoly = rp.precomputed, # RolyPoly object precomputed (linked GWAS and gene annotation data)
                                                                        block_data = block_data)) 


# ======================= TMPTMP ======================= #

rp.tmp <- tryCatch(
  {
    tmp.time <- system.time(rolypoly <- rolypoly:::rolypoly_load_block_data(rolypoly = rp.precomputed, # RolyPoly object precomputed (linked GWAS and gene annotation data)
                                                                            block_data = block_data)) # load expression data
  },
  error=function(cond) {
    message("An error happened")
    message("Error: ", cond)
    message() # gives extra space.
  } # end error handling
) # end tryCatch
