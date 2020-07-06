#!/usr/bin/bash

# Download CELLECT and run the mtag_munge.py script from the ldsc subdir

# https://github.com/perslab/CELLECT/wiki/CELLECT-LDSC-Tutorial
# setup up and activate munge_ldsc environment before running this command

DIR_CELLECT=/projects/jonatan/tools/sc-genetics/CELLECT
N=50 # number of files to process in parallel

if [ -d "./GWAS_munged" ]; then
  mkdir ./GWAS_munged
fi


#filename=./NULL_GWAS/plink_out.1KG_phase3_EUR_null_gwas.P97.qassoc_with_alleles
for filename in ./plink_out/*qassoc_with_alleles*; do
	((i=i%N)); ((i++==0)) && wait
	gwas_id=$(echo $filename | egrep -o 'P[[:digit:]]+' | head -n1)
	echo $filename
	echo $gwas_id
	python ${DIR_CELLECT}/ldsc/mtag_munge.py \
   --sumstats $filename \
 	 --merge-alleles $DIR_CELLECT/data/ldsc/w_hm3.snplist \
   --snp SNP \
	 --n-value 503 \
   --keep-pval \
	 --a1 A1 \
	 --a2 A2 \
   --p P \
	 --stdout-off \
   --out ./GWAS_munged/1KG_phase3_EUR_null_gwas.${gwas_id}  & #> ./GWAS_munged.1KG_phase3_EUR_null_gwas.${gwas_id}_log.txt
done
