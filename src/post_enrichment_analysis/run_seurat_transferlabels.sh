#!/usr/bin/bash

################# SET PARAMETERS #################

flag_date="191203"

# test data
dir_data="../../tmp-data/expression/"
declare -a id_datExpr_tests=('moffitt2018' 'campbell2017' 'mikkelsen2019') 
#declare -a path_seuratObj_tests=('tmp-data/expression/moffitt2018_seuratObj_subset.RDS.gz' 'tmp-data/expression/campbell2017_seuratObj_subset.RDS.gz' 'tmp-data/expression/mikkelsen2019_seuratObj_subset.RDS.gz' 'tmp-data/expression/mousebrain_seuratObj_subset.RDS.gz')
#declare -a id_metadata_tests=('moffitt2018.metadata' 'campbell2017.metadata' 'mikkelsen2019.metadata' 'mousebrain.metadata')
#metadata_test_cell_id_col='cell_id'
declare -a metadata_test_annot_cols=('cell_type' 'cell_type_all_lvl2' 'cell_type')
declare -a vec_metadata_test_annots_subsets=('c("i9_Gaba","e8_Cck_Ebf3")' 'c("n29.Nr5a1-Adcyap1")' 'c("Glut_5")')
# 'c("DEINH3","MEGLU11","MEINH2","MEGLU10","DEGLU5","MEGLU1","MEGLU3","TEGLU23","DEGLU4","TEGLU19","HBSER5","HBGLU2","MEGLU2","TEGLU17","MEGLU7","MEINH3","HBGLU5","TEINH12","TEINH1","TEGLU4","TEGLU21","HBINH8")') 
# ref data

###
#declare -a id_datExpr_refs=('moffitt2018' 'campbell2017' 'mikkelsen2019', 'mousebrain')
###

#declare -a path_seuratObj_refs=('tmp-data/expression/moffitt2018_seuratObj.RDS.gz' 'tmp-data/expression/campbell2017_seuratObj.RDS.gz' 'tmp-data/expression/mikkelsen2019_seuratObj.RDS.gz' 'tmp-data/expression/mousebrain_seuratObj.RDS.gz')

###
#declare -a id_metadata_refs=('moffitt2018.metadata' 'campbell2017.metadata' 'mikkelsen2019.metadata' 'mousebrain.metadata')
###

#metadata_ref_cell_id_col='cell_id'
#declare -a metadata_ref_annot_cols=('cell_type' 'cell_type_all_lvl2' 'cell_type' 'ClusterName') 
#'ClusterName')
#declare -a vec_metadata_ref_annots_subsets=('c("i9_Gaba","e8_Cck_Ebf3")' 'c("n29.Nr5a1-Adcyap1")' 'c("Glut_5")' 'c("DEINH3","MEGLU11","MEINH2","MEGLU10","DEGLU5","MEGLU1","MEGLU3","TEGLU23","DEGLU4","TEGLU19","HBSER5","HBGLU2","MEGLU2","TEGLU17","MEGLU7","MEINH3","HBGLU5","TEINH12","TEINH1","TEGLU4","TEGLU21","HBINH8")')

# other parameters
nComp=50
#regressOut_percent_mt=TRUE
output_label="CELLECTprior_vs_CELLECTprior"
append_results=TRUE

################# RUN R SCRIPT #################

echo "date: " ${flag_date}

# we need to loop over each test dataset 3 times 
declare -a vec_k=(0 1 2)
#declare -a vec_k=(0 0 0 0 1 1 1 1 2 2 2 2)
# we need to loop 4 times over the vector of ref datasets

#declare -a vec_j=(0 1 2 3 0 1 2 3 0 1 2 3)

for (( i=0; i<${#vec_k[@]}; i++ )); do
    
    k=${vec_k[$i]}
    j=${vec_j[$i]}
    
    echo 
    echo "i =" ${i} "k =" ${k} "j =" ${j}
    
    ### TEST DATASET ###
    
  	id_datExpr_test=${id_datExpr_tests[$k]}
	  echo "id_datExpr_test: " ${id_datExpr_test}
	  
	  #path_seuratObj_test=${path_seuratObj_tests[$k]}
	  
	  #id_metadata_test=${id_metadata_tests[$k]}
	  
	  metadata_test_annot_col=${metadata_test_annot_cols[$k]}
	  
	  vec_metadata_test_annots_subset=${vec_metadata_test_annots_subsets[$k]}
	  echo "vec_metadata_test_annots_subset: " ${vec_metadata_test_annots_subset}
	  
	  ### REF DATASET ###
	  id_datExpr_ref=mikkelsen2019
	  #id_datExpr_ref=${id_datExpr_refs[$j]}
	  echo "id_datExpr_ref :" ${id_datExpr_ref}
	  
	  #path_seuratObj_ref=${path_seuratObj_refs[$j]}
	  
	  #id_metadata_ref=${id_metadata_refs[$j]}
	  metadata_ref_annot_col=cell_type
	  #metadata_ref_annot_col=${metadata_ref_annot_cols[$j]}
	  
	  #vec_metadata_ref_annots_subset=${vec_metadata_ref_annots_subsets[$j]}
	  #echo "vec_metadata_ref_annots_subset: " ${vec_metadata_ref_annots_subset}	  
	  
	  time Rscript ./seurat_transferlabels.R --dir_data ${dir_data} \
	                                         --id_datExpr_test ${id_datExpr_test} \
	                                         --metadata_test_annot_col ${metadata_test_annot_col} \
	                                         --vec_metadata_test_annots_subset ${vec_metadata_test_annots_subset} \
	                                         --id_datExpr_ref ${id_datExpr_ref} \
	                                         --metadata_ref_annot_col ${metadata_ref_annot_col}  \
	                                         --nComp ${nComp} \
	                                         --output_label ${output_label} \
	                                         --append_results ${append_results}

done


echo "bash script done!"



