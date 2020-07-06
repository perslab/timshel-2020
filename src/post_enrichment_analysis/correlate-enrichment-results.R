############### SYNOPSIS ###################
# Run correlation analysis in CELLECT and rare / mendelian variant results

### OUTPUT:
#' TODO
#' 
### REMARKS:
# ....

### REFERENCE:
# * 

if (F) {
  opt= list(
    path_CELLECT_results = "/projects/jonatan/pub-perslab/timshel-bmicelltypes2019/results/cellect_ldsc/prioritization.csv",
    dir_geneset_results= "/projects/jonatan/pub-perslab/timshel-bmicelltypes2019/out/es_enrichment_test/per_geneset/",
    specificity_id=  "campbell2017_lvl2",
    geneset_results_regex = "BMI_23_genes",
    col_geneset_results = "p.value",
    col_CELLECT = "pvalue",
    vec_geneset_name = 'c("protein_mono_early_extr_obesity")',
    vec_gwas = 'c("BMI_UKBB_Loh2018")',
    method= "pearson",
    minLog10transform = T,
    alternative = "two-sided",
    filename_out = "cor_CELLECT_geneset_MNS.csv",
    append_results=  T
  )
}
######################################################################
########################### OptParse #################################
######################################################################

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--path_CELLECT_results", type="character",
              help = "character, absolute path to prioritization.csv file with CELLECT outs, [default %default]"), 
  make_option("--dir_geneset_results", type="character", 
              help = "character, absolute path to dir containing ES geneset enrichment results, [default %default]"), 
  make_option("--specificity_id", type="character",
              help = "e.g. 'mousebrain'"), 
  make_option("--geneset_results_regex", type="character",
              help = "string to use alongside specificity_id to match correct file in dir_geneset_results"), 
  make_option("--col_geneset_results", type="character", default= "p.value",
              help = "which column of the geneset results should be correlated? "), 
  make_option("--col_CELLECT", type="character", default= "pvalue",
              help = "which column of the CELLECT results should be correlated? "), 
  make_option("--vec_geneset_name", type="character", 
              help = "quoted vector of character, select this in geneset_name column of the geneset results csv. If geneset_name column doesn't exist it uses all rows [default $default]"), 
  make_option("--vec_gwas", type="character",
              help = "quoted vector of gwas among CELLECT results to include [default $default]"), 
  make_option("--method", type="character", default= "pearson",
              help = "correlation coefficient, passed to stats::cor, one of pearson, kendall, spearman, [default $default]"), 
  make_option("--minLog10transform", type="logical", default= T,
              help = "transform values by -log10 before correlating? [default $default]"), 
  make_option("--alternative", type="character", default= "two.sided",
              help = "alternative hypothesis for p.value [default $default]"), 
  make_option("--filename_out", type="character",
              help = "character, name for output files (including the filetype suffix)"),
  make_option("--append_results", type="logical", default=TRUE,
              help = "If file already exists, append results? [default %default] ")
)


opt <- parse_args(OptionParser(option_list=option_list))

path_CELLECT_results  <- opt$path_CELLECT_results
dir_geneset_results <- opt$dir_geneset_results
specificity_id <- opt$specificity_id
geneset_results_regex <- opt$geneset_results_regex
vec_geneset_name <- eval(parse(text=opt$vec_geneset_name))
col_geneset_results <- opt$col_geneset_results
col_CELLECT <- opt$col_CELLECT
vec_gwas <- eval(parse(text=opt$vec_gwas))
method <- opt$method
minLog10transform <- opt$minLog10transform
alternative = opt$alternative
filename_out <- opt$filename_out
append_results <- opt$append_results


######################################################################
########################### PACKAGES #################################
######################################################################

message("Loading packages")

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("broom"))
suppressPackageStartupMessages(library("openxlsx"))

######################################################################
########################### SET OPTIONS ##############################
######################################################################

options(stringsAsFactors = F, 
        use="pairwise.complete.obs", 
        warn=1, 
        verbose=F,
        mc.cores=40 # for parallel computation
) 

######################################################################
############################ CONSTANTS ###############################
######################################################################

flag_date = substr(gsub("-","",as.character(Sys.Date())),3,1000)

randomSeed = 12345
set.seed(randomSeed)

######################################################################
########################## LOAD FILES ################################
######################################################################

dt_CELLECT <- fread(path_CELLECT_results)

path_geneset_results <- dir(path=dir_geneset_results, 
                            pattern=geneset_results_regex, 
                            full.names=T)
path_geneset_results <- grep(pattern = specificity_id, 
                             x= path_geneset_results,
                             value=T)
path_geneset_results <- grep(pattern = "\\.csv\\.gz", 
                             x= path_geneset_results,
                             value=T)
dt_geneset_results <- fread(path_geneset_results)

######################################################################
########################## FILTER AND REORDER  #######################
######################################################################

list_dt_cor.test <- lapply(vec_gwas, function(gwas) {
  
  dt_cor.test = data.table(t(sapply(vec_geneset_name, function(geneset_name){
    
    # cell_type,statistic,parameter,p.value,p.value_emp,alternative
    vec_geneset_value <- as.numeric(dt_geneset_results[[col_geneset_results]])

    names(vec_geneset_value) <- dt_geneset_results[["annotation"]]
    
    if (!is.null(dt_geneset_results$geneset_name)) vec_geneset_value <- vec_geneset_value[dt_geneset_results$geneset_name == geneset_name]
    
    df_CELLECT = setDF(dt_CELLECT)
    # gwas,specificity_id,annotation,tau,se,pvalue
    vec_CELLECT_value <- df_CELLECT[df_CELLECT$gwas == gwas & df_CELLECT$specificity_id == specificity_id,][[col_CELLECT]] %>% as.numeric
      
      #[eval(condition),][[col_CELLECT]] %>% as.numeric
    names(vec_CELLECT_value) <- df_CELLECT[df_CELLECT$gwas == gwas & df_CELLECT$specificity_id == specificity_id,][["annotation"]]
    
    
    # check the two vectors
    stopifnot(all(names(vec_CELLECT_value) %in%  names(vec_geneset_value)))
    stopifnot(length(vec_CELLECT_value) == length(vec_geneset_value))
    
    # put vectors in same order
    vec_CELLECT_value <- vec_CELLECT_value[match(names(vec_geneset_value),names(vec_CELLECT_value))]
    
    # logtransform p values
    if (minLog10transform) {
      vec_geneset_value <- -log10(vec_geneset_value)
      vec_CELLECT_value <- -log10(vec_CELLECT_value)
    }
    
    # compute correlation
    # cor(x=vec_geneset_value,
    #     y=vec_CELLECT_value,
    #     method=method)
    # 
    cor.test.out <- cor.test(x=vec_geneset_value,
                                     y=vec_CELLECT_value,
                                     method=method,
                                     alternative=alternative)
    # tbbl_out <-broom::tidy(cor.test(x=vec_geneset_value,
    #         y=vec_CELLECT_value,
    #         method=method,
    #         alternative=alternative))
    # vec_out<-as.character(tbbl_out)
    
    vec_out = c("estimate"=cor.test.out$estimate,
                "statistic"=cor.test.out$statistic,
                "p.value"=cor.test.out$p.value,
                "parameter"=cor.test.out$parameter,
                "conf.low"=cor.test.out$conf.low,
                "conf.high"=cor.test.out$conf.high,
                "method"=cor.test.out$method,
                "alternative"=cor.test.out$alternative)
    
    
    #names(vec_out)<-colnames(tbbl_out)    
    vec_out
    
  })))
  dt_cor.test = data.table("gwas"=gwas, "specificity_id" = specificity_id, "geneset_name"=vec_geneset_name,dt_cor.test)

})

#print(list_dt_cor.test)
#print(list_cor)
# mat_cor <- matrix(data=mat_cor, nrow=length(vec_geneset_name), ncol = length(vec_gwas))
# colnames(mat_cor) <- vec_gwas
# 
# dt_cor <- data.table(
#   geneset_name = vec_geneset_name,
#   mat_cor 
# )
# 
# dt_cor_long <- melt.data.table(dt_cor, 
#                     id.vars= "geneset_name",
#                     variable.name = "gwas",
#                     value.name = paste0(method, "_cor"))

dt_cor_long = Reduce(f = rbind, x=list_dt_cor.test)

dt_cor_long$minLog10_vals <- minLog10transform

######################################################################
#############################  WRAP UP  ##############################
######################################################################

file.out.results <- here("out", "post_enrichment_analysis", filename_out)
fwrite(x = dt_cor_long, file = file.out.results,  append = if (file.exists(file.out.results)) append_results else F, nThread=24, verbose=T)
message("script done!")

