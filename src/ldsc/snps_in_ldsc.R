

# ### LDSCORE FILES
# Chr 20
# 36301 Hapmap SNPs
# 31360 NN individual ANNOT (no pandas concat)
# 31360 MACA COMBINED ANNOT
# 31102 BASELINE
# 31102 1000G_Phase3_cell_type_groups
# 31102 Cahoy CONTROL
# 31102 Cahoy 3
# 31102 baselineLD_v2.0
# 31102 baseline_v1.1
# 221626 1000G BIM
# Chr 21
# 19331 Hapmap SNPs
# 16994 NN individual ANNOT (no pandas concat)
# 16994 MACA COMBINED ANNOT
# 17041 BASELINE
# 17041 1000G_Phase3_cell_type_groups
# 17041 Cahoy CONTROL
# 17041 Cahoy 3
# 17041 baselineLD_v2.0
# 17041 baseline_v1.1
# 138712 1000G BIM
# Chr 22
# 20106 Hapmap SNPs
# 17262 NN individual ANNOT (no pandas concat)
# 17262 MACA COMBINED ANNOT
# 17490 BASELINE
# 17490 1000G_Phase3_cell_type_groups
# 17490 Cahoy CONTROL
# 17490 Cahoy 3
# 17490 baselineLD_v2.0
# 17490 baseline_v1.1
# 141123 1000G BIM

chromosome <- 22
file.hm <- sprintf("/raid5/projects/timshel/sc-genetics/ldsc/data/hapmap3_snps/hm.%s.snp", chromosome)
file.baseline <- sprintf("/raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_baseline/baseline.%s.l2.ldscore.gz", chromosome)
file.1kg <- sprintf("/raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC.%s.bim", chromosome)
file.hm_print_snps <- "/raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_baseline/print_snps.txt"
file.hm_merge_alleles <- "/raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist"


df.hm <- read_table(file.hm,col_names="SNP")
df.hm_print_snps <- read_table(file.hm_print_snps,col_names="SNP")
df.hm_merge_alleles <- read_tsv(file.hm_merge_alleles)
df.baseline <- read_tsv(file.baseline) %>% select(CHR, SNP, BP)
df.1kg <- read_tsv(file.1kg, col_names=c("CHR", "SNP", "CM", "BP", "A1", "A2"))


length(df.baseline$SNP) # 17489
length(df.1kg$SNP) # 141123
length(df.hm$SNP) # 20106
length(df.hm_print_snps$SNP) # 1217311
length(df.hm_merge_alleles$SNP) # 1217311 --> we use these to filter the GWAS sumstats
sum(df.hm$SNP %in% df.baseline$SNP) # 15958
sum(df.baseline$SNP %in% df.hm$SNP) # 15958 | Not all baseline SNPs are in hapmap per chromosome file. This means that the baseline model cannot have used these SNPs
sum(df.baseline$SNP %in% df.1kg$SNP) # 17489 | Ok, all baseline SNPs are in 1kg.
sum(df.hm$SNP %in% df.1kg$SNP) # 17261 | Hmmm, ~3000 of the hapmap per chromosome SNPs are missing from the 1kg. 
sum(df.baseline$SNP %in% df.hm_print_snps$SNP) # 17489 | HURRA, all baseline SNPs are in the hapmap print SNPs!!!! 



# ======================= COMMANDS ======================= #
# ### LDSCORE files
# for i in {20..22}; do
# echo Chr $i
# echo `cat /raid5/projects/timshel/sc-genetics/ldsc/data/hapmap3_snps/hm.$i.snp | wc -l` "Hapmap SNPs"
# echo `zcat /scratch/sc-ldsc/nn_lira_sema/nn_lira_sema.yellow.$i.l2.ldscore.gz | wc -l` "NN individual ANNOT (no pandas concat)"
# #echo `zcat /scratch/sc-ldsc/nn_lira_sema_1000kb/nn_lira_sema.COMBINED_ANNOT.$i.l2.ldscore.gz | wc -l` "NN COMBINED ANNOT 1000kb"
# #echo `zcat /scratch/sc-ldsc/nn_lira_sema-combined/nn_lira_sema.COMBINED_ANNOT.$i.l2.ldscore.gz | wc -l` "NN COMBINED ANNOT 100kb"
# echo `zcat /scratch/sc-ldsc/maca/maca_tissue_cell_type.COMBINED_ANNOT.$i.l2.ldscore.gz | wc -l` "MACA COMBINED ANNOT"
# echo `zcat /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_baseline/baseline.$i.l2.ldscore.gz | wc -l` "BASELINE"
# echo `zcat /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_Phase3_cell_type_groups/cell_type_group.9.$i.l2.ldscore.gz | wc -l` "1000G_Phase3_cell_type_groups"
# echo `zcat /raid5/projects/timshel/sc-genetics/ldsc/data/LDSC_SEG_ldscores/Cahoy_1000Gv3_ldscores/Cahoy.control.$i.l2.ldscore.gz | wc -l` "Cahoy CONTROL"
# echo `zcat /raid5/projects/timshel/sc-genetics/ldsc/data/LDSC_SEG_ldscores/Cahoy_1000Gv3_ldscores/Cahoy.3.$i.l2.ldscore.gz | wc -l` "Cahoy 3"
# echo `zcat /raid5/projects/timshel/sc-genetics/ldsc/data/baselineLD_v2.0/baselineLD.$i.l2.ldscore.gz | wc -l` "baselineLD_v2.0"
# echo `zcat /raid5/projects/timshel/sc-genetics/ldsc/data/baseline_v1.1/baseline.$i.l2.ldscore.gz | wc -l` "baseline_v1.1"
# echo `cat /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i.bim | wc -l` "1000G BIM"
# done
# 
# 
# ### ANNOT
# for i in {20..22}; do
# echo Chr $i
# echo `zcat /scratch/sc-ldsc/nn_lira_sema/nn_lira_sema.yellow.$i.annot.gz | wc -l` "NN individual ANNOT (no pandas concat)"
# echo `zcat /scratch/sc-ldsc/nn_lira_sema_1000kb/nn_lira_sema.COMBINED_ANNOT.$i.annot.gz | wc -l` "NN COMBINED ANNOT 1000kb"
# echo `zcat /scratch/sc-ldsc/nn_lira_sema-combined/nn_lira_sema.COMBINED_ANNOT.$i.annot.gz | wc -l` "NN COMBINED ANNOT 100kb"
# echo `zcat /scratch/sc-ldsc/maca/maca_tissue_cell_type.COMBINED_ANNOT.$i.annot.gz | wc -l` "MACA COMBINED ANNOT"
# echo `zcat /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_baseline/baseline.$i.annot.gz | wc -l` "BASELINE"
# echo `zcat /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_Phase3_cell_type_groups/cell_type_group.9.$i.annot.gz | wc -l` "1000G_Phase3_cell_type_groups"
# echo `zcat /raid5/projects/timshel/sc-genetics/ldsc/data/LDSC_SEG_ldscores/Cahoy_1000Gv3_ldscores/Cahoy.control.$i.annot.gz | wc -l` "Cahoy CONTROL"
# echo `zcat /raid5/projects/timshel/sc-genetics/ldsc/data/LDSC_SEG_ldscores/Cahoy_1000Gv3_ldscores/Cahoy.3.$i.annot.gz | wc -l` "Cahoy 3"
# echo `zcat /raid5/projects/timshel/sc-genetics/ldsc/data/baselineLD_v2.0/baselineLD.$i.annot.gz | wc -l` "baselineLD_v2.0"
# echo `zcat /raid5/projects/timshel/sc-genetics/ldsc/data/baseline_v1.1/baseline.$i.annot.gz | wc -l` "baseline_v1.1"
# echo `cat /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i.bim | wc -l` "1000G BIM"
# done
