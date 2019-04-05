############### SYNOPSIS ###################
# Extract list of RolyPoly GWAS SNPs

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

############################################

wd <- "/projects/timshel/sc-genetics/sc-genetics/src/misc/"
setwd(wd)

# ======================================================================================================= #
# ==========================================       SETUP      =========================================== #
# ======================================================================================================= #


# suppressMessages(library(rolypoly))
suppressMessages(library(tidyverse))
library(stringr)

dir.ldfiles <- "/projects/timshel/sc-genetics/sc-genetics/data/rolypoly/EUR_LD_FILTERED_NONAN_R"

# ======================================================================================================= #
#=============================================== READ ROLYPOLY LD data ================================================== #
# ======================================================================================================= #


## RolyPoly code
# ld_data <- readRDS(chrom_ld_file_path)[(R**2 > r2_threshold), .(SNP_A, SNP_B, R)]
# setkey(ld_data, SNP_A) # just in case

### INFO
# Each .Rds file contrains a data.table / data.frame object.
# CHR_A      BP_A       SNP_A     MAF_A CHR_B      BP_B       SNP_B     MAF_B           R
# 1:     1     10177 rs367896724 0.4055670     1     10352 rs555500075 0.4264410  0.21699800
# 2:     1     10352 rs555500075 0.4264410     1     11008 rs575272151 0.0884692 -0.11132800

path.ldfiles <- file.path(dir.ldfiles, list.files(dir.ldfiles, pattern=".Rds")) # full path
# names(path.ldfiles) <- sapply(str_split(basename(path.ldfiles), pattern="\\."), '[[', 1)  # extracting chromosome info. alternative --> tools::file_path_sans_ext()
path.ldfiles
list.df_ld_data <- lapply(path.ldfiles, readRDS) # read files. 1.4 GB RAM.
df.ld_data <- bind_rows(list.df_ld_data)
df.ld_data.fmt <- bind_rows(df.ld_data %>% select(CHR_A, SNP_A), 
                            df.ld_data %>% select(CHR_B, SNP_B) %>% rename(CHR_A=CHR_B, SNP_A=SNP_B)) # renaming
df.ld_data.fmt.unique <- df.ld_data.fmt %>% distinct() # this takes < 30 sec


head(df.ld_data.fmt.unique)
# CHR_A       SNP_A
# 1:     1 rs367896724
# 2:     1 rs555500075
# 3:     1 rs575272151
nrow(df.ld_data.fmt.unique) # TOTAL NUMBER OF SNPs---> 8,737,808
# n_distinct(df.ld_data.fmt.unique$SNP_A) # --> 8,737,807 (should be same as nrow, but we do not care...)




# ======================================================================================================= #
#=============================================== EXPORT ================================================== #
# ======================================================================================================= #

# write_csv(df.ld_data.fmt.unique, path="rolypoly_ld_snps.csv.gz")

# ======================================================================================================= #
#============================================ SEMI-colon SNPs =========================================== #
# ======================================================================================================= #

### *OBS* *NOT UNDERSTOOD*:  some entries contain multiple SNPs (semi-colon sepereret)

df.multi_line_snps <- df.ld_data %>% filter(grepl(";", SNP_A) | grepl(";", SNP_B))
# nrow(df.multi_line_snps) # --> 175,319
# CHR_A    BP_A                   SNP_A     MAF_A CHR_B    BP_B       SNP_B     MAF_B       R
# 1     1 1232206 rs544125251;rs557094073 0.0149105     1 1232233 rs115237750 0.0238569 0.81915
# 2     1 1232206 rs544125251;rs557094073 0.0149105     1 1232229 rs116263241 0.0238569 0.81915

names(df.ld_data)
