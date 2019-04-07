############### SYNOPSIS ###################
# Analyze LDSC WGCNA results on TB+MB modules with kme cut-offs

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# =============================== SETUP ================================ #
# ======================================================================= #

library(tidyverse)

wd <- "/projects/timshel/sc-genetics/sc-genetics/src/wgcna_modules/"
setwd(wd)

# ======================================================================= #
# =============================== PARAMS ================================ #
# ======================================================================= #

### TB
# dir.data <- "/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/"
# filenames <- list.files(path=dir.data,  pattern="wgcna.tabula_muris-190111.fdr_sign_celltypes_min_kme_(.*).continuous__BMI_Yengo2018.cell_type_results.txt")
# filenames <- c("wgcna.tabula_muris-190111.fdr_sign_celltypes.continuous__BMI_Yengo2018.cell_type_results.txt")

### MB
dir.data <- "/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/"
filenames <- list.files(path=dir.data,  pattern="wgcna.mousebrain-190111.fdr_sign_celltypes_min_(.*).continuous__BMI_Yengo2018.cell_type_results.txt")
filenames <- c(filenames, "wgcna.mousebrain-190111.fdr_sign_celltypes.continuous__BMI_Yengo2018.cell_type_results.txt")

### Load
list.dfs <- lapply(file.path(dir.data, filenames), read_tsv)
tmp_match <- stringr::str_match(filenames, pattern="min_(.*).continuous__BMI_Yengo2018.cell_type_results.txt")[,2] # ALT: filenames
tmp_match[is.na(tmp_match)] <- "kme_0"
names(list.dfs) <- tmp_match
names(list.dfs)
df.ldsc_cts <- list.dfs %>% bind_rows(.id="cutoff")

### Process
df.ldsc_cts <- df.ldsc_cts %>% mutate(module_id = str_split_fixed(Name, "__", n=Inf)[,2]) %>% select(-Name, -Coefficient, -Coefficient_std_error)
df.ldsc_cts.spread <- df.ldsc_cts %>% spread(key = cutoff, value = Coefficient_P_value)

### Plot violon
ggplot(df.ldsc_cts, aes(x=cutoff, y=-log10(Coefficient_P_value))) + geom_violin()

### Plot top modules
df.ldsc_cts.spread.top10 <- df.ldsc_cts.spread %>% arrange(kme_0) %>% slice(1:10)
df.ldsc_cts.spread.top10 <- df.ldsc_cts.spread.top10 %>% gather(key="kme", value="p", -module_id)
p <- ggplot(df.ldsc_cts.spread.top10, aes(x=kme, y=-log10(p), group=module_id, color=module_id, label=module_id)) + geom_line()
p + labs(x="kme cut-off", y="-log10(LDSC P-val)", title="LDSC significance as a function of module kME cut-off")
# p + geom_label_repel(nudge_x = 1,na.rm = TRUE)

library(ggrepel)
# ======================================================================= #
# =============================== FUNCTIONS ================================ #
# ======================================================================= #
