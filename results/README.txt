

###################################### External file source ######################################

### ES gene enrichments
cp /projects/jonatan/applied/18-mousebrain_7/tables/mb_GSA_SEMall_WGCNAtop_2_analyt_mat_wilcoxonPval_genesetTests.txt es_enrichment--wgcna_modules.pvals.txt

### Module kME table
cp /projects/jonatan/applied/18-mousebrain_7/tables/Neurons_sub_ClusterName_7.2_run1_kMs_full_join.csv.gz modules--kme_table.csv.gz

### Module genes and metadata
cp /scratch/sc-ldsc/wgcna.mousebrain-190213.fdr_sign_celltypes.continuous/log.wgcna.mousebrain-190213.fdr_sign_celltypes.continuous.multi_geneset.txt modules--metadata.txt

### Module genetic prioritization: BMI, all cell-types
cp /projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/wgcna.mousebrain-190213.fdr_sign_celltypes.continuous__BMI_UKBB_Loh2018.cell_type_results.txt prioritization_modules--mousebrain.BMI_UKBB_Loh2018.csv.gz

