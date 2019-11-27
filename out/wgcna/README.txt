
# ======================================================================= #
# ==================== FILE ORIGIN INFO 19.11.27 ======================= #
# ======================================================================= #

Date: 19.11.27
WGCNA on n=22 MB BMI cell-types
cp /projects/jonatan/pub-perslab/18-mousebrain/18-mousebrain_7/tables/ClusterName_prior_191127a_cell_cluster_module_genes.csv.gz modules.mousebrain_bmicelltypes.rwgcna_table.csv.gz
M1 is back as floralwhite_427




# ======================================================================= #
# ==================== FILE ORIGIN INFO [OUTDATED] ====================== #
# ======================================================================= #

### Modules from FDR significant cell-types [v3, 190213] | NEW: WGCNA run on n=11 BMI_UKBB_Loh2018 FDR cell-types 
# Feb 13 [6:15 PM]
# Hej Pascal,
# Nu er den nye mb WGCNA analyse  endelig klar. Jeg gik tilbage til den tidligere udgave af pipelinen fra Januar.
# Koerslen hedder Neurons_sub_ClusterName_7.2_run1 og den relevante output fil er
# /projects/jonatan/mousebrain_7/tables/Neurons_sub_ClusterName_7.2_run1_cell_cluster_module_genes.csv.gz
# Dodgerblue modulet gaar igen uden aendring, nu under navnet “lavenderblush”.
# (har i senere aendringer i scriptet soerget for at alle random seeds, inklusive navne, er reproducible) 
# dict_genomic_annot = {"wgcna.mousebrain-190213.fdr_sign_celltypes.continuous": 
# 					  	{"dataset":"mousebrain", 
# 					  	"file_multi_gene_set":"/projects/jonatan/mousebrain_7/tables/Neurons_sub_ClusterName_7.2_run1_cell_cluster_module_genes.csv.gz"}
# 					 }
# New file path [August 2019] ----> /projects/jonatan/pub-perslab/18-mousebrain/18-mousebrain_7/tables/Neurons_sub_ClusterName_7.2_run1_cell_cluster_module_genes.csv.gz

### Shell commands to generate out/wgcna/modules.mousebrain_bmicelltypes.rwgcna_table.csv
### NB: we filter out TEGLU4 (not significant after fix)
# FILE_M=/projects/jonatan/pub-perslab/18-mousebrain/18-mousebrain_7/tables/Neurons_sub_ClusterName_7.2_run1_cell_cluster_module_genes.csv.gz
# zcat $FILE_M | grep -v TEGLU4 | gzip > /projects/timshel/sc-genetics/timshel-bmicelltypes2019/out/wgcna/modules.mousebrain_bmicelltypes.rwgcna_table.csv.gz
