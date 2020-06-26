############### SYNOPSIS ###################
# Download and prepare GTEx v8 bulk gene expression data

### OUTPUT: 
# normalized gene * sample matrix, one version with all samples, another only with brain
# metadata

### REMARKS:
# ....

### REFERENCE:

# * sample annotation overview: https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx
# * Gene-level expression
# *Based on collapsed GENCODE annotation
# * Quantified with RNA-SeQC
# * Read counts
# * *No covariate correction*
  
# # Filtering
#   
# * The GTEx Consortium. Genetic effects on gene expression across human tissues. https://doi.org/10.1038/nature24277 (2017)
# 
# # Normalization
# 
# * DESeq2 vignette: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts
# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #

library("here")
library("tidyverse")
library("data.table")
library("Biobase", lib.loc = "~/R/x86_64-pc-linux-gnu-library/")
library("SummarizedExperiment", lib.loc = "~/R/x86_64-pc-linux-gnu-library/")
library("DESeq2", lib.loc = "~/R/x86_64-pc-linux-gnu-library/")
#library("openxlsx")
library("WGCNA")
# ======================================================================= #
# ================================ CONSTANTS ================================ #
# ======================================================================= #

path_gene_counts_gct = here("tmp-data", "expression", "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
path_sample_attrs = here("tmp-data", "expression", "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

# ======================================================================= #
# ========================== Download data ================================ #
# ======================================================================= #
### You only need to run this step once

if (!file.exists(path_gene_counts_gct)){ 
  download.file(url = "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz", destfile=path_gene_counts_gct)
} 

if (!file.exists(path_sample_attrs)) {
  download.file(url = "https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", destfile=path_sample_attrs)
}

# ======================================================================= #
# ========================== Read in data ================================ #
# ======================================================================= #

dt_gene_counts = fread(input = path_gene_counts_gct, sep = "\t")
dt_sample_attrs = fread(path_sample_attrs) 

dt_sample_attrs_RNASEQ = dt_sample_attrs[SMAFRZE=="RNASEQ"]

# ======================================================================= #
# ========================== modify strings ================================ #
# ======================================================================= #

# replace non-safe characters with safe ones 

dt_sample_attrs_RNASEQ$SMTS <- tolower(gsub("-| |\\(|\\)","_", dt_sample_attrs_RNASEQ$SMTS))
dt_sample_attrs_RNASEQ$SMTS <- gsub("___","_", dt_sample_attrs_RNASEQ$SMTS)
dt_sample_attrs_RNASEQ$SMTS <- gsub("^_|_$","", dt_sample_attrs_RNASEQ$SMTS)
dt_sample_attrs_RNASEQ$SMTS <- gsub("__","_", dt_sample_attrs_RNASEQ$SMTS)

dt_sample_attrs_RNASEQ$SMTSD <- tolower(gsub("-| |\\(|\\)","_", dt_sample_attrs_RNASEQ$SMTSD))
dt_sample_attrs_RNASEQ$SMTSD <- gsub("___","_", dt_sample_attrs_RNASEQ$SMTSD)
dt_sample_attrs_RNASEQ$SMTSD <- gsub("^_|_$","", dt_sample_attrs_RNASEQ$SMTSD)
dt_sample_attrs_RNASEQ$SMTSD <- gsub("__","_", dt_sample_attrs_RNASEQ$SMTSD)

colnames(dt_sample_attrs_RNASEQ) = gsub("__", "_", colnames(dt_sample_attrs_RNASEQ))

# ======================================================================= #
# ========================== Filter data ================================ #
# ======================================================================= #

# convert reads table to matrix

mat_gene_counts = as.matrix(dt_gene_counts[,3:ncol(dt_gene_counts)])
rownames(mat_gene_counts) = dt_gene_counts$Name 

# check metadata and expression matrix match

all.equal(dt_sample_attrs_RNASEQ$SAMPID, colnames(mat_gene_counts))
# [1] TRUE

### filter samples 

# filter on variable value cutoffs as in 
# The GTEx Consortium. Genetic effects on gene expression across human tissues. https://doi.org/10.1038/nature24277 (2017)

vec_logical_sampleOK = with(data = dt_sample_attrs_RNASEQ, expr= {
  SMMPPD >= 10e6 &# mapped reads (Mapped: Total number of reads aligned/mapped)
    SMMAPRT >= 0.2 &# read mapping rate (Mapping Rate: Ratio of total mapped reads to total reads)
    SMNTERRT <= 0.3 &# intergenic mapping rate (Intergenic Rate: The fraction of reads that map to the genomic space between genes)
    SMBSMMRT <= 0.008 & # base mismatch rate (Base Mismatch Rate: Number of bases not matching the reference divided by the total number of bases aligned.)
    SMRRNART <= 0.3})# rRNA read rate (rRNA Rate: Ratio of all reads aligned to rRNA regions to total reads)

table(vec_logical_sampleOK)
# FALSE  TRUE 
#   173 17209 

mat_gene_counts = mat_gene_counts[,vec_logical_sampleOK]
dt_sample_attrs_RNASEQ = dt_sample_attrs_RNASEQ[vec_logical_sampleOK,]

# Additionally, outlier samples were identiﬁed based on expression proﬁle using a correlation-based statistic and sex incompatibility checks, following methods described in [14].

#Wright et al, Nature (2014) Heritability and genomics of gene expression in peripheral blood  
#https://www-nature-com/articles/ng.2951

#"Second, we examined the pairwise correlation matrix of expression profiles. Using rij as the correlation between arrays i and j, we computed
#the average correlation of array i with all others of the total n arrays. Lower values correspond to lower quality and were expressed in terms of median absolute deviations to provide a sense of distance from the grand correlation mean."

#Wright et al, Nature (2014) Heritability and genomics of gene expression in peripheral blood  https://www.nature.com/articles/ng.2951

enableWGCNAThreads(nThreads = 20)

IQ_range_multiplier_cutoff = 1.5

list_mat_SMTSD_filter <- lapply(unique(dt_sample_attrs_RNASEQ$SMTSD), function(SMTSD_lvl){
  mat_gene_counts_SMTSD_sub = mat_gene_counts[,match(dt_sample_attrs_RNASEQ$SAMPID[dt_sample_attrs_RNASEQ$SMTSD==SMTSD_lvl], colnames(mat_gene_counts))]
  mat_cor_tmp = WGCNA::cor(x = mat_gene_counts_SMTSD_sub, use = "pairwise.complete.obs", nThreads = 20)
  mat_cor_tmp = mat_cor_tmp - Matrix::diag(x=1, nrow=nrow(mat_cor_tmp))
  vec_cor_mean = colMeans2(x = mat_cor_tmp) 
  vec_dist = vec_cor_mean - mean(vec_cor_mean)
  vec_dist_scaled = vec_dist/median(abs(vec_dist))
  names(vec_dist_scaled) = colnames(mat_gene_counts_SMTSD_sub)
  
  Q1 = as.numeric(quantile(x=vec_dist_scaled, probs=0.25))
  Q3 = as.numeric(quantile(x=vec_dist_scaled, probs=0.75))
  outlier_threshold = Q1-IQ_range_multiplier_cutoff*(Q3-Q1)
  vec_outliers = names(vec_dist_scaled)[vec_dist_scaled<outlier_threshold]
  
  return(mat_gene_counts_SMTSD_sub[,!colnames(mat_gene_counts_SMTSD_sub) %in% vec_outliers])
})

mat_gene_counts_filtOutliers = Reduce(x=list_mat_SMTSD_filter, f=cbind)

dt_sample_attrs_RNASEQ_filtOutliers = dt_sample_attrs_RNASEQ[dt_sample_attrs_RNASEQ$SAMPID %in% colnames(mat_gene_counts_filtOutliers)] 

vec_order = match(colnames(mat_gene_counts_filtOutliers), dt_sample_attrs_RNASEQ_filtOutliers$SAMPID)
dt_sample_attrs_RNASEQ_filtOutliers = dt_sample_attrs_RNASEQ_filtOutliers[vec_order]

### filter genes

#use read cutoff used in 
#The GTEx Consortium. Genetic effects on gene expression across human tissues. https://doi.org/10.1038/nature24277 (2017)

vec_logical_keep <- rowSums(mat_gene_counts_filtOutliers >= 6)>=10 # at least six reads in at least ten samples
table(vec_logical_keep)
# FALSE  TRUE 
#  9807 46393 

mat_gene_counts_filtOutliers <- mat_gene_counts_filtOutliers[vec_logical_keep,]

## Normalization

#DESeq2: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts

#Do normalized counts correct for variables in the design?
#  No. The design variables are not used when estimating the size factors, and counts(dds, normalized=TRUE) is providing counts scaled by size or normalization factors. The design is only used when estimating dispersion and log2 fold changes.

#Correct within broader tissue (SMTS level) to avoid inducing normalization-related effects (trade-off..)

#"To perform the median of ratios method of normalization, DESeq2 has a single estimateSizeFactors() function that will generate size factors for us. We will use the function in the example below, but in a typical RNA-seq analysis this step is automatically performed by the DESeq() function, which we will see later"

design = as.formula(paste0("~ 1")) # just a dummy, not using it in DESEq

list_mat_SMTS_normLog = lapply(unique(dt_sample_attrs_RNASEQ_filtOutliers$SMTS), 
                               function(SMTS_lvl){
                                 vec_logical = dt_sample_attrs_RNASEQ_filtOutliers$SMTS==SMTS_lvl
                                 dds <- DESeqDataSetFromMatrix(
                                   countData  = mat_gene_counts_filtOutliers[,vec_logical],
                                   colData = dt_sample_attrs_RNASEQ_filtOutliers[vec_logical,],
                                   design = design)
                                 dds <- estimateSizeFactors(dds)
                                 mat_gene_counts_filtOutliers_norm <- counts(dds, normalized=TRUE)
                                 mat_gene_counts_filtOutliers_normLog <- log2(mat_gene_counts_filtOutliers_norm + 1)
                                 rownames(mat_gene_counts_filtOutliers_normLog) = rownames(mat_gene_counts_filtOutliers_normLog)
                                 return(mat_gene_counts_filtOutliers_normLog)
                               })

names(list_mat_SMTS_normLog) = unique(dt_sample_attrs_RNASEQ_filtOutliers$SMTS)

mat_gene_counts_filtOutliers_normLog = Reduce(x=list_mat_SMTS_normLog, f=cbind)

vec_order = match(colnames(mat_gene_counts_filtOutliers_normLog), dt_sample_attrs_RNASEQ_filtOutliers$SAMPID)
dt_sample_attrs_RNASEQ_filtOutliers = dt_sample_attrs_RNASEQ_filtOutliers[vec_order]

#### clear up gene names

#remove version number appended to gene names (GENCODE)
# https://www.gencodegenes.org/pages/data_format.html

# would there be any duplicates? 

rownames(mat_gene_counts_filtOutliers_normLog)[duplicated(substr(x=rownames(mat_gene_counts_filtOutliers_normLog),start =1,stop=15))]
# [1] character(0)
rownames(mat_gene_counts_filtOutliers_normLog) = substr(x=rownames(mat_gene_counts_filtOutliers_normLog),1,15)

# ======================================================================= #
# =============== EXPORT CELL-TYPE/ANNOTATION METADATA TO CSV =========== #
# ======================================================================= #


df.metadata <- dt_sample_attrs[,c("sample_id", "tissue", "tissue_sub", "RIN"):=.(SAMPID,SMTS,SMTSD,SMRIN)]
df.metadata <- dt_sample_attrs_RNASEQ_filtOutliers %>% as_tibble %>% distinct(tissue, tissue_sub)
df.metadata <- df.metadata %>% mutate(annotation=str_replace_all(tolower(make.names(tissue_sub)), "(\\.)+", "_"))
df.metadata <- df.metadata %>% mutate(annotation=str_replace_all(annotation, "_$", ""))
df.metadata <- df.metadata %>% select(annotation, everything())

### Export
df.metadata %>% write_csv(here("data/expression/gtex/gtex-metadata.csv"))

# ======================================================================= #
# ===================== EXPORT DATA AND ANNOTATIONS TO TMP ===================== #
# ======================================================================= #

dt_gene_counts_filtOutliers_normLog = data.table("gene"=rownames(mat_gene_counts_filtOutliers_normLog),mat_gene_counts_filtOutliers_normLog)

file.out.data = here("tmp-data","expression","gtex_v8_filtered_normLog_all.csv")
fwrite(x =  dt_gene_counts_filtOutliers_normLog, file = file.out.data)
R.utils::gzip(file.out.data, overwrite=TRUE) # gzip

# write cell metadata to tmp
file.out.meta <- here("tmp-data","expression","gtex_v8_filtered_normLog.metadata_all.csv")
data.table::fwrite(dt_sample_attrs_RNASEQ_filtOutliers, file=file.out.meta,  # fwrite cannot write gziped files
                   nThread=24, verbose=T) # write file ---> write to scratch

# file.out.meta = here("tmp-data", "expression", paste0(dirOut,prefixData,"_annotations_",params$date,".csv"))
# fwrite(dt_attrs_combined_RNASEQ_sub_filtOutliers, file = file.out.meta)