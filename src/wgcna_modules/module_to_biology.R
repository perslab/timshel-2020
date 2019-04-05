############### SYNOPSIS ###################
# From WGCNA module to biology

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# =============================== SETUP ================================ #
# ======================================================================= #

library(tidyverse)
library(gProfileR)
# library(GOplot)

wd <- "/projects/timshel/sc-genetics/sc-genetics/src/wgcna_modules/"
setwd(wd)

# ======================================================================= #
# =============================== PARAMS ================================ #
# ======================================================================= #

# name.dataset <- "maca"
# name.dataset <- "mousebrain_Neurons"
name.dataset <- "mousebrain"
# name.dataset <- "hypothalamus_mette_thesis"

# ======================================================================= #
# =============================== FUNCTIONS ================================ #
# ======================================================================= #

### Map:  ENTREZ --> ENSEMBL
map_entrez2ensembl <- function(df.magma) {
  file.mapping <- "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_id_mapping.hsapiens.ensembl_entrez.txt.gz"
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

file.magma_gwas <- "/scratch/tmp-magma_gwas/BMI_Yengo2018.txt.10UP.1.5DOWN.genes.out"
# file.gwas_loci <- "/projects/timshel/DEPICT/BMI_Yengo2018/results/BMI_Yengo2018.1e-5.depict_tissues_loci.txt" # 1e-5
file.gwas_loci <- "/projects/timshel/DEPICT/BMI_Yengo2018/results/BMI_Yengo2018.5e-8.depict_tissues_loci.txt" # 5e-8
file.genes_mendelian <- "/projects/timshel/sc-genetics/sc-genetics/data/genes_obesity/turcot2018_s21_mendelian_obesity_genes.mapped.txt" # human_genes = ensembl_gene_id
file.genes_rare_variant <- "/projects/timshel/sc-genetics/sc-genetics/data/genes_obesity/turcot2018_table1_rare_variants.mapped.txt" # human_genes = ensembl_gene_id
file.genes_mouse_obesity <- "/projects/timshel/sc-genetics/sc-genetics/data/genes_obesity/yazdi2015_table1_mouse_obesity_genes.mapped.txt" # human_genes = ensembl_gene_id

file.ortholog_genes_background <- "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz" # human_genes = ensembl_gene_id

# ======================================================================= #
# ============================ *SWITCH* ================================ #
# ======================================================================= #

if (name.dataset == "maca") {
  ### MACA
  file.ldsc_cts <- "/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/wgcna.maca.BMI_Yengo2018.cell_type_results.txt"
  file.module_geneset <- "/scratch/sc-ldsc/maca/log.maca_tissue_cell_type.multi_geneset.txt"
  
  file.module_origin_metadata <- "/projects/timshel/sc-genetics/sc-genetics/src/GE-maca/maca_facs.tissue_cell_type.metadata.csv"
  df.module_origin_metadata <- read_csv(file.module_origin_metadata) %>% select(module_origin=tissue_cell_type,
                                                                                module_origin_ncells=n, 
                                                                                module_origin_desc=tissue)
} else if (name.dataset == "mousebrain_Neurons") {
  ### Mousebrain Neurons
  file.ldsc_cts <- "/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/wgcna.mousebrain_Neurons.BMI_Yengo2018.cell_type_results.txt"
  file.module_geneset <- "/scratch/sc-ldsc/mousebrain_Neurons/log.Neurons_ClusterName.multi_geneset.txt"
  file.module_origin_metadata <- "/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain-agg_L5.metadata.csv"
  df.module_origin_metadata <- read_csv(file.module_origin_metadata) %>% select(module_origin=ClusterName,
                                                                                module_origin_ncells=NCells, 
                                                                                module_origin_desc=Description)
} else if (name.dataset == "mousebrain") {
  ### Mousebrain
  file.module_geneset <- "/scratch/sc-ldsc/log.mousebrain_ALL_CLASSES.multi_geneset.txt"
  # wgcna_run       cell_cluster    module_color_id      gene_input      hgnc    pkME    gene    annotation
  # Astrocytes      ACBG    blue1   ENSMUSG00000038292      Ccdc155 0.962783147428202       ENSG00000161609 Astrocytes.blue1
  # Astrocytes      ACBG    blue1   ENSMUSG00000027368      Dusp2   0.962783147428202       ENSG00000158050 Astrocytes.blue1

  
  file.module_origin_metadata <- "/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain-agg_L5.metadata.csv"
  df.module_origin_metadata <- read_csv(file.module_origin_metadata) %>% select(module_origin=ClusterName,
                                                                                module_origin_ncells=NCells, 
                                                                                module_origin_desc=Description)
  
} else if (name.dataset == "hypothalamus_mette_thesis") {
  ### hypothalamus_mette_thesis
  file.ldsc_cts <- "/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/wgcna.hypothalamus_mette_thesis.BMI_Yengo2018.cell_type_results.txt"
  file.module_geneset <- "/scratch/sc-ldsc/hypothalamus_mette_thesis/log.hypothalamus_mette_thesis.multi_geneset.txt"
  
  df.module_origin_metadata <- tibble(module_origin="DUMMY",
                                      module_origin_ncells="DUMMY", 
                                      module_origin_desc="DUMMY")
}


# ======================================================================= #
# ============================ READ DATA ================================ #
# ======================================================================= #

### Orthologs
df.ortholog_genes_background <- read_tsv(file.ortholog_genes_background) # contains 1-1 unique orthologs


### MAGMA
df.magma_gwas.raw <- map_entrez2ensembl(read_table(file.magma_gwas)) # magma annotation window: 10kb up, 1.5kb down
# [1] "Number of genes mapped: 17467"
# [1] "Number of genes not mapped: 158"
df.magma_gwas <- df.magma_gwas.raw %>% transmute(ensembl_gene_id=ensembl_gene_id,
                                             magma_pval = P,
                    magma_pval_rank_all_genes = rank(P)) %>%
  filter(ensembl_gene_id %in% df.ortholog_genes_background$ensembl_gene_id) %>% # filter to ortholog genes
  mutate(magma_pval_rank_orthologs = rank(magma_pval)) 
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
df.module_geneset

df.module_metadata <- df.module_geneset %>% select(module_origin, module_id) %>% distinct()

### LDSC module results
if (name.dataset == "mousebrain") {
  files <- Sys.glob("/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/wgcna.mousebrain_*.BMI_Yengo2018.cell_type_results.txt")
  df.ldsc_cts <- lapply(files, read_tsv) %>% bind_rows()
  df.ldsc_cts <- df.ldsc_cts %>% mutate(module_id=stringr::str_replace_all(Name, pattern="mousebrain_", "")) # SPECIFIC FOR 'MOUSEBRAIN ALL': e.g mousebrain_Astrocytes.purple2 --> Astrocytes.purple2
} else {
  df.ldsc_cts <- read_tsv(file.ldsc_cts)
  df.ldsc_cts <- df.ldsc_cts %>% mutate(module_id=stringr::str_split_fixed(Name,pattern="\\.",n=2)[,2]) # add module ID from 'Name' string, e.g. maca_tissue_cell_type.slateblue4
}
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



# ### OLD non-parallel (works, but may give time-out error)
# df.gprofiler <- df.module_geneset %>% 
#   filter(module_id %in% top_modules) %>% # filter top modules
#   group_by(module_id) %>%
#   arrange(pkME) %>% # order by kME
#   do(gProfileR::gprofiler(.$ensembl_gene_id, organism="hsapiens", ordered_query=T, significant=T, custom_bg=df.ortholog_genes_background$ensembl_gene_id))

### Filter TF domain
# df.gprofiler <- df.gprofiler %>% filter(!domain=="tf")


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
### add ldsc module p-value and filter
df.bio_genes <- df.bio_genes %>% 
  left_join(df.ldsc_cts %>% select(module_id, module_ldsc_pval=P), by="module_id") %>%
  filter(module_ldsc_pval <= 0.05) # only write out module gene-based info for semi-significant modules (otherwise the list gets too long)
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

# ======================================================================= #
# =============================== Module biology summary ================================ #
# ======================================================================= #

names(df.bio_genes)

df.bio_modules <- df.bio_genes %>% 
  group_by(module_id) %>%
  summarise(
    module_origin=unique(module_origin),
    #module_ldsc_pval=unique(module_ldsc_pval), # we join ldsc results later, so no need for this
    n_genes_module=n(),
    n_genes_gwas_loci=sum(flag_gene_gwas_loci),
    n_genes_mendelian=sum(flag_gene_mendelian),
    n_genes_rare_variant=sum(flag_gene_rare_variant),
    n_genes_mouse_obesity=sum(flag_gene_mouse_obesity),
    mean_minus_log10_magma_pval=mean(-log10(magma_pval), na.rm=T)
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

xlsx.workbook <- createWorkbook(creator="PTimshel") # start excel
do.excel_export(df.gprofiler.ordered.meta, sheet_name="GO analysis - GSEA", xlsx.workbook)
do.excel_export(df.gprofiler.unordered.meta, sheet_name="GO analysis", xlsx.workbook)
do.excel_export(df.bio_genes, sheet_name="Module bio. gene-based", xlsx.workbook)
do.excel_export(df.bio_modules, sheet_name="Module bio. summary", xlsx.workbook)
saveWorkbook(xlsx.workbook, file = sprintf("out.module_to_biology.%s.xlsx", name.dataset), overwrite = TRUE) # write file



# Write to excel file

# ======================================================================= #
# =============================== XXXXX ================================ #
# ======================================================================= #



# ======================================================================= #
# =============================== LEFTOVER ================================ #
# ======================================================================= #



