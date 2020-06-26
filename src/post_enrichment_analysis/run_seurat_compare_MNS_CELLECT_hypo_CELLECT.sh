#!/usr/bin/bash

################# SET PARAMETERS #################

flag_date="191203"

# test data
dir_data="../../tmp-data/expression/"
declare -a id_datExpr_tests=('mousebrain') 
declare -a metadata_test_annot_cols=('ClusterName')
declare -a vec_metadata_test_annots_subsets=('c("DEINH3","MEGLU11","MEINH2","MEGLU10","DEGLU5","MEGLU1","MEGLU3","TEGLU23","DEGLU4","TEGLU19","HBSER5","HBGLU2","MEGLU2","TEGLU17","MEGLU7","MEINH3","HBGLU5","TEINH12","TEINH1","TEGLU4","TEGLU21","HBINH8")')

# ref data
declare -a id_datExpr_refs=('moffitt2018' 'campbell2017' 'mikkelsen2019')
declare -a metadata_ref_annot_cols=('cell_type' 'cell_type_all_lvl2' 'cell_type') 


# other parameters
nComp=50
#regressOut_percent_mt=TRUE
output_label="MNSCELLECT_vs_hypoCELLECT"
append_results=TRUE

################# RUN R SCRIPT #################

echo "date: " ${flag_date}


declare -a vec_k=(0 0 0)
declare -a vec_j=(0 1 2)

for (( i=0; i<${#vec_k[@]}; i++ )); do
    
    k=${vec_k[$i]}
    j=${vec_j[$i]}
    
    echo 
    echo "i =" ${i} "k =" ${k} "j =" ${j}
    
    ### TEST DATASET ###
    
  	id_datExpr_test=${id_datExpr_tests[$k]}
	  echo "id_datExpr_test: " ${id_datExpr_test}
	  metadata_test_annot_col=${metadata_test_annot_cols[$k]}
	  
	  vec_metadata_test_annots_subset=${vec_metadata_test_annots_subsets[$k]}
	  echo "vec_metadata_test_annots_subset: " ${vec_metadata_test_annots_subset}
	  
	  ### REF DATASET ###
 
	  id_datExpr_ref=${id_datExpr_refs[$j]}
	  echo "id_datExpr_ref :" ${id_datExpr_ref}
	  metadata_ref_annot_col=${metadata_ref_annot_cols[$j]}
	  
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



