############### SYNOPSIS ###################
# From WGCNA module to biology
# VERSION 2: update to '__' convention.

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# =============================== SETUP ================================ #
# ======================================================================= #

library(here)
library(tidyverse)
library(gProfileR)
# library(GOplot)

setwd(here("src/wgcna_modules"))

# ======================================================================= #
# =============================== PARAMS ================================ #
# ======================================================================= #

gwas <- "BMI_UKBB_Loh2018"
# gwas <- "BMI_UPDATE_Yengo2018"

name.dataset <- "mousebrain"
# prefix_genomic_annot <- "wgcna.mousebrain-190218.fdr_sign_celltypes.continuous" # deepsplit1
prefix_genomic_annot <- "wgcna.mousebrain-190213.fdr_sign_celltypes.continuous" # lavenderblush
### prefix_genomic_annot <- "wgcna.mousebrain-190111.fdr_sign_celltypes.continuous"
### prefix_genomic_annot <- "wgcna.mousebrain-181214.fdr_sign_celltypes.continuous"


# name.dataset <- "tabula_muris"
# prefix_genomic_annot <- "wgcna.tabula_muris-190111.fdr_sign_celltypes.continuous"
### prefix_genomic_annot <- "wgcna.tabula_muris-181214.fdr_sign_celltypes.continuous"

# ======================================================================= #
# =============================== FUNCTIONS ================================ #
# ======================================================================= #

### Map:  ENTREZ --> ENSEMBL
map_entrez2ensembl <- function(df.magma) {
  file.mapping <- here("data/gene_annotations/gene_id_mapping.hsapiens.ensembl_entrez.txt.gz")
  df.mapping <- read_tsv(file.mapping)

  genes_mapped <- df.mapping$ensembl_gene_id[match(df.magma$GENE, df.mapping$entrezgene)]
  print(sprintf("Number of genes mapped: %s",sum(!is.na(genes_mapped))))
  print(sprintf("Number of genes not mapped: %s",sum(is.na(genes_mapped)))) # number of not mapped genes
  df.magma.clean <- df.magma %>% 
    mutate(ensembl_gene_id=genes_mapped) %>%
    filter(!is.na(ensembl_gene_id))
  return(df.magma.clean)
}


do_gprofiler <- function(df, ordered_query){
  ### INPUT: data frame with columns: 
  # pkME
  # ensembl_gene_id
  # module_id
  ### NOTES
  # df.ortholog_genes_background$ensembl_gene_id is also used (but not parsed to the function)
  
  df <- df %>% arrange(pkME) # make sure genes are sorted by kME
  module_id <- unique(df$module_id) 
  stopifnot(length(module_id) == 1) # there should only be one module_id
  
  df.res = tryCatch({ # to avoid gprofiler error "transfer closed with outstanding read data remaining"
    message(sprintf("running gprofiler for module_id=%s", module_id))
    df.res <- gProfileR::gprofiler(df$ensembl_gene_id, organism="hsapiens", ordered_query=ordered_query, significant=T, custom_bg=df.ortholog_genes_background$ensembl_gene_id)
    ### Custom background (not ordered)
    # df.res <- gProfileR::gprofiler(df$ensembl_gene_id, organism="hsapiens", ordered_query=F, significant=T, custom_bg=df.ortholog_genes_background$ensembl_gene_id)
    ### Ordered
    # df.res <- gProfileR::gprofiler(df$ensembl_gene_id, organism="hsapiens", ordered_query=T, significant=T)
    df.res <- df.res %>% mutate(call_status="Successful")
    return(df.res)
  }, warning = function(w) {
    message(sprintf("warning: %s", w))
  }, error = function(e) {
    message(sprintf("error: %s", e))
    return(data.frame(call_status="Failed")) # return single-row data frame
  }, finally = {
    message(sprintf("done with module_id=%s", module_id))
  })
  return(df.res)
}


# ======================================================================= #
# =============================== FILES INPUT ================================ #
# ======================================================================= #

# file.magma_gwas <- "/scratch/tmp-magma_gwas/BMI_Yengo2018.txt.10UP.1.5DOWN.genes.out"
# file.magma_gwas <- "/nfsdata/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_magma/BMI_Yengo2018.txt.10UP.1.5DOWN.genes.out"
file.magma_gwas <- here("out/magma/gene_based_scores/BMI_UKBB_Loh2018_no_mhc.resid_correct_all.gsa.genes.mapped.out") # 100 KB window
# file.gwas_loci <- "/projects/timshel/DEPICT/BMI_Yengo2018/results/BMI_Yengo2018.1e-5.depict_tissues_loci.txt" # 1e-5
# file.gwas_loci <- "/projects/timshel/DEPICT/BMI_Yengo2018/results/BMI_Yengo2018.5e-8.depict_tissues_loci.txt" # 5e-8
file.gwas_loci <- here("out/depict/results/BMI_UKBB_Loh2018_no_mhc.5e-8.depict_tissues_loci.txt")

file.genes_mendelian <- here("data/genes_obesity/turcot2018_s21.mendelian_obesity_genes.mapped.txt") # human_genes = ensembl_gene_id
file.genes_rare_variant <- here("data/genes_obesity/turcot2018_table1.rare_variants.mapped.txt") # human_genes = ensembl_gene_id
file.genes_mouse_obesity <- here("data/genes_obesity/yazdi2015_table1.mouse_obesity_genes.mapped.txt") # human_genes = ensembl_gene_id

file.ortholog_genes_background <- here("data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz") # human_genes = ensembl_gene_id

# ======================================================================= #
# ============================ *SWITCH* ================================ #
# ======================================================================= #

### Set LDSC specific files
file.ldsc_cts <- sprintf("/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/%s__%s.cell_type_results.txt", prefix_genomic_annot, gwas)
file.module_geneset <- sprintf("/scratch/sc-ldsc/%s/log.%s.multi_geneset.txt", prefix_genomic_annot, prefix_genomic_annot)


if (name.dataset == "tabula_muris") {
  ### MACA
  file.module_origin_metadata <- "/projects/timshel/sc-genetics/sc-genetics/data/expression/tabula_muris/tabula_muris_facs.tissue_celltype.celltype_metadata.csv"
  df.module_origin_metadata <- read_csv(file.module_origin_metadata) %>% select(module_origin=tissue_celltype,
                                                                                module_origin_ncells=n_cells, 
                                                                                module_origin_desc=tissue)
  df.module_origin_metadata <- df.module_origin_metadata %>% mutate(module_origin = stringr::str_replace_all(module_origin, pattern="\\s+", replacement="_")) 
} else if (name.dataset == "mousebrain") {
  ### Mousebrain
  file.module_origin_metadata <- "/projects/timshel/sc-genetics/sc-genetics/data/expression/mousebrain/mousebrain-agg_L5.metadata.csv"
  df.module_origin_metadata <- read_csv(file.module_origin_metadata) %>% select(module_origin=annotation,
                                                                                module_origin_ncells=NCells, 
                                                                                module_origin_desc=Description)
}


# ======================================================================= #
# ============================ READ DATA ================================ #
# ======================================================================= #

### Orthologs
df.ortholog_genes_background <- read_tsv(file.ortholog_genes_background) # contains 1-1 unique orthologs


### MAGMA
df.magma_gwas.raw <- map_entrez2ensembl(read_tsv(file.magma_gwas)) # magma annotation window: 10kb up, 1.5kb down
# [1] "Number of genes mapped: 17467"
# [1] "Number of genes not mapped: 158"
df.magma_gwas <- df.magma_gwas.raw %>% transmute(ensembl_gene_id=ensembl_gene_id,
                                             magma_zstat = ZSTAT,
                    magma_zstat_rank_all_genes = rank(ZSTAT)) %>%
  filter(ensembl_gene_id %in% df.ortholog_genes_background$ensembl_gene_id) %>% # filter to ortholog genes
  mutate(magma_zstat_rank_orthologs = rank(magma_zstat)) 
# df.magma_gwas %>% arrange(P)

### GWAS LOCI
df.gwas_loci <- read_tsv(file.gwas_loci)
genes_in_loci.rows <- df.gwas_loci %>% pull(genes_in_locus)
genes_in_loci.ext <- stringr::str_split(paste(genes_in_loci.rows,collapse=";"), pattern=";")[[1]]
df.genes_in_gwas_loci <- tibble(ensembl_gene_id=unique(genes_in_loci.ext)) # 5e-8=1695 genes, 1e-5=3396 genes

### OBESITY GENES: MENDELIAN / RARE / MOUSE
df.genes_mendelian <- read_tsv(file.genes_mendelian)
df.genes_rare_variant <- read_tsv(file.genes_rare_variant)
df.genes_mouse_obesity <- read_tsv(file.genes_mouse_obesity)

### GENESET
df.module_geneset <- read_tsv(file.module_geneset)
df.module_geneset <- df.module_geneset %>% rename(ensembl_gene_id=gene, module_id=annotation, module_origin=cell_cluster)
df.module_geneset <- df.module_geneset %>% mutate(module_origin = stringr::str_replace_all(module_origin, pattern="\\s+", replacement="_"))  # *HACK FOR TABULA MURIS*
df.module_geneset <- df.module_geneset %>% rename(pkME=annotation_value) # *NEW*
### GETSET METADATA
df.module_metadata <- df.module_geneset %>% select(module_origin, module_id) %>% distinct()

### LDSC module results
df.ldsc_cts <- read_tsv(file.ldsc_cts)
mat_tmp_split_str <- stringr::str_split_fixed(df.ldsc_cts$Name,pattern="\\__",n=Inf)
df.ldsc_cts <- df.ldsc_cts %>% mutate(module_id=mat_tmp_split_str[,ncol(mat_tmp_split_str)]) # add module ID from 'Name' string, e.g. maca_tissue_cell_type.slateblue4
df.ldsc_cts <- df.ldsc_cts %>% rename(BETA=Coefficient, SE=Coefficient_std_error, P=Coefficient_P_value)

# ======================================================================= #
# =============================== TRAIT-specificity ================================ #
# ======================================================================= #



# ======================================================================= #
# =============================== GO analysis ================================ #
# ======================================================================= #

### Filter top10 LDSC modules
n.top_modules <- 15
top_modules <- df.ldsc_cts %>% top_n(n.top_modules, -P) %>% pull(module_id)


### split data frame
list_of_dfs_per_module <- df.module_geneset %>% 
  filter(module_id %in% top_modules) %>% # filter top modules
  base::split(.$module_id) # ALT1: plyr::dlply(df, "V1", identity)
names(list_of_dfs_per_module) # --> names of modules

### Run GO analysis
# TODO: consider making this a %dopar% loop. http://www.vikparuchuri.com/blog/monitoring-progress-inside-foreach-loop/
library(parallel)
cl <- makeCluster(n.top_modules, type = "FORK") # outfile="" (no redirection) prints to console [only works for socket clusters?]
clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(gProfileR))
list_of_dfs_per_module.gprofiler.ordered <- parLapply(cl, list_of_dfs_per_module, do_gprofiler, ordered_query=T) # PARALLEL VERSION. returns a list of data frames (not named - only 'indexes')
df.gprofiler.ordered <- bind_rows(list_of_dfs_per_module.gprofiler.ordered, .id="module_id")
list_of_dfs_per_module.gprofiler.unordered <- parLapply(cl, list_of_dfs_per_module, do_gprofiler, ordered_query=F) # PARALLEL VERSION. returns a list of data frames (not named - only 'indexes')
df.gprofiler.unordered <- bind_rows(list_of_dfs_per_module.gprofiler.unordered, .id="module_id")
stopCluster(cl)
# names(list_of_dfs_per_module.gprofiler) <- names(list_of_dfs_per_module) #: not needed because do_gprofiler() never returns NULL

### gprofiler WIKI
# term.size: the number of genes in a term, 
# query.size: number of recognized genes in the input query 
# ^----> I cannot make sense of the query.size numbers. I observe that for the same the module, the "query.size" changes for different GO terms analyzed. 
# ^----> The query.size should be the same for every GO term analyzed for a module, since it is the same set of query genes. (Custom background option does not change this)
# ^----> SOLUTION/EXPLANATION: the query.size changes within a module when running with "ordered_query=T".
# overlap.size: the overlap of the query and term. 
# If a custom statistical background has been passed, these sets have been intersected with the background, resulting in a value smaller or equal than without a background.
# precision: or predictive rate, i.e. the proportion of genes in the input list that have the function.
# recall: i.e. the proportion of functionally annotated genes that the query recovers.


### Filter TF domain
# df.gprofiler <- df.gprofiler %>% filter(!domain=="tf")

df.gprofiler.ordered.meta <- df.gprofiler.ordered %>% 
  filter(!domain=="tf") %>%
  left_join(df.module_metadata, by=c("module_id")) %>% 
  select(module_origin, -query.number, everything())

df.gprofiler.unordered.meta <- df.gprofiler.unordered %>% 
  filter(!domain=="tf") %>%
  left_join(df.module_metadata, by=c("module_id")) %>% 
  select(module_origin, -query.number, everything())


# ======================================================================= #
# =============================== Gene based biology  ================================ #
# ======================================================================= #

df.bio_genes <- df.module_geneset # copy
### add kme rank
df.bio_genes <- df.bio_genes %>% 
  group_by(module_id) %>% 
  mutate(kme_rank_within_module=rank(-pkME)) %>%
  ungroup()
### add ldsc module p-value
df.bio_genes <- df.bio_genes %>% 
  left_join(df.ldsc_cts %>% select(module_id, module_ldsc_pval=P), by="module_id")

### add gene flags
df.bio_genes <- df.bio_genes %>% mutate( 
  flag_gene_gwas_loci = ensembl_gene_id %in% df.genes_in_gwas_loci$ensembl_gene_id,
  flag_gene_mendelian = ensembl_gene_id %in% df.genes_mendelian$ensembl_gene_id,
  flag_gene_rare_variant = ensembl_gene_id %in% df.genes_rare_variant$ensembl_gene_id,
  flag_gene_mouse_obesity = ensembl_gene_id %in% df.genes_mouse_obesity$ensembl_gene_id
  )
### add magma gene pvalues
df.bio_genes <- df.bio_genes %>% left_join(df.magma_gwas, by="ensembl_gene_id")
### Order
df.bio_genes <- df.bio_genes %>% arrange(module_ldsc_pval, kme_rank_within_module)
### copy and filter (we copy to make be able to list all modules in df.bio_modules)
df.bio_genes.all <- df.bio_genes
df.bio_genes <- df.bio_genes %>% filter(module_ldsc_pval <= 0.05) # only write out module gene-based info for esi-significant modules (otherwise the list gets too long)

# ======================================================================= #
# =============================== Module biology summary ================================ #
# ======================================================================= #

names(df.bio_genes)

df.bio_modules <- df.bio_genes.all %>% 
  group_by(module_id) %>%
  summarise(
    module_origin=unique(module_origin),
    #module_ldsc_pval=unique(module_ldsc_pval), # we join ldsc results later, so no need for this
    n_genes_module=n(),
    n_genes_gwas_loci=sum(flag_gene_gwas_loci),
    n_genes_mendelian=sum(flag_gene_mendelian),
    n_genes_rare_variant=sum(flag_gene_rare_variant),
    n_genes_mouse_obesity=sum(flag_gene_mouse_obesity),
    mean_magma_zstat=mean(magma_zstat, na.rm=T)
  )

### add ALL modules from ldsc (full join) [this makes us independent on the module_ldsc_pval filtering we did for the df.bio_genes]
df.bio_modules <- df.bio_modules %>% full_join(df.ldsc_cts %>% select(module_id, module_ldsc_pval=P), by="module_id")
### add meta data (e.g. n_cells) | full join
df.bio_modules <- df.bio_modules %>% full_join(df.module_origin_metadata, by="module_origin")
### Order
df.bio_modules <- df.bio_modules %>% arrange(module_ldsc_pval)

# ======================================================================= #
# =============================== Export Excel ================================ #
# ======================================================================= #

library(openxlsx) # better than xlsx library because it does not depend on Java (which causes "java.lang.OutOfMemoryError" errors). Also openxlsx is being actively developed

do.excel_export <- function(df, sheet_name, xlsx.workbook) {
  # write df to excel
  # remove_existing_excel_sheet(xlsx.workbook, sheet_name)
  addWorksheet(wb=xlsx.workbook, sheetName=sheet_name)
  writeData(wb=xlsx.workbook, sheet=sheet_name, as.data.frame(df), colNames=T, rowNames=F)
}


# save.image(file = sprintf("out.module_to_biology.%s.RData", name.dataset))

### Write to excel file
xlsx.workbook <- createWorkbook(creator="PTimshel") # start excel
do.excel_export(df.gprofiler.ordered.meta, sheet_name="GO analysis - GSEA", xlsx.workbook)
do.excel_export(df.gprofiler.unordered.meta, sheet_name="GO analysis", xlsx.workbook)
do.excel_export(df.bio_genes, sheet_name="Module bio. gene-based", xlsx.workbook)
do.excel_export(df.bio_modules, sheet_name="Module bio. summary", xlsx.workbook)
file.out <- here("results/", sprintf("modules--biology.%s-%s.xlsx", gwas, prefix_genomic_annot))
saveWorkbook(xlsx.workbook, file=file.out, overwrite=TRUE) # write file




# ======================================================================= #
# =============================== XXXXX ================================ #
# ======================================================================= #



# ======================================================================= #
# =============================== LEFTOVER ================================ #
# ======================================================================= #



