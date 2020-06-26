#!/usr/bin/bash

# calls the correlate-enrichment-results.R script 

#### GWAS (Fig 3A) ####

# BMI (BMI_UKBB_Loh2018)
# Educational Attainment (EA3_Lee2018)
# Intelligence (INTELLIGENCE_Savage2018)
# Schizophrenia (SCZ_Pardinas2018)
# Insomnia (INSOMNIA_Jansen2018)
# Multiple Schlerosis (MS_Patsopoulos2011) #TODO check 
# Rheumatoid Arthritis (RA_Okada2014) #TODO check 
# Low Density Lipoprotein (LIPIDS_LDL_Teslovich2010) #TODO check 
# Waist-hip Ratio (adj. BMI) (WHRadjBMI_UKBB_Loh2018) 
# Height (HEIGHT_UKBB_Loh2018) 


flag_date="200617"
declare -a specificity_ids=("mousebrain" "campbell2017_lvl2" "chen2017" "romanov2017" "mikkelsen2019" "moffitt2018" "kimVMH2019_smartseq")

for (( i=0; i<${#specificity_ids[@]}; i++ )); do
  
  specificity_id=${specificity_ids[$i]}
  echo $specificity_id
  
  Rscript ./correlate-enrichment-results.R --path_CELLECT_results ../../results/cellect_ldsc/prioritization.csv --dir_geneset_results ../../out/es_enrichment_test/per_geneset/  --specificity_id $specificity_id --geneset_results_regex BMI_23_genes  --col_geneset_results p.value   --col_CELLECT pvalue  --vec_geneset_name 'c("protein_mono_early_extr_obesity")'  --vec_gwas 'c("BMI_UKBB_Loh2018")'  --method pearson  --alternative two.sided --filename_out cor_CELLECT_geneset_results_${flag_date}.csv --append_results T
  

done
# vec_GWAS 'c("BMI_UKBB_Loh2018", "EA3_Lee2018", "INTELLIGENCE_Savage2018", "SCZ_Pardinas2018",  "INSOMNIA_Jansen2018", "MS_Patsopoulos2011","RA_Okada2014", "LIPIDS_LDL_Teslovich2010","WHRadjBMI_UKBB_Loh2018","HEIGHT_UKBB_Loh2018")' 
echo "bash script done"
