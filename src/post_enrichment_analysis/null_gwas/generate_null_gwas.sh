#!usr/bin/bash

# generate qassoc files 
# requires ./genotypes/ALL.chr_merged.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bim, too large to include in repo
# need to adjust plink path
################## RUN PLINK 1.9 ##################

# point towards plink installation
export PLINK_EXEC="/raid5/projects/timshel/bin/plink1.9_2015-11-26/plink"


export BFILE="./phenotypes/ALL.chr_merged.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.fam"
export PHENO="./phenotypes/null_gwas_guassian_phenotypes.n1000.phe"
export PLINK_OUT="./plink_out/plink_out.1KG_phase3_EUR_null_gwas"

if [ -d "./plink_out" ]; then
  mkdir plink_out
fi

$PLINK_EXEC --bfile $BFILE --threads 50 --assoc --pheno $PHENO --allow-no-sex --all-pheno --out $PLINK_OUT

################## Adding allele information ##################
### REASON: needed for LDSC munge summary stats

for i in {1..1000}; do
	echo $i
        FILE_OUT=./plink_out/plink_out.1KG_phase3_EUR_null_gwas.P${i}.qassoc_with_alleles
        echo -e "A1\tA2" > $FILE_OUT.tmp
        cut -f 5-6 ./genotypes/ALL.chr_merged.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bim >> $FILE_OUT.tmp
        paste ./plink_out/plink_out.1KG_phase3_EUR_null_gwas.P${i}.qassoc $FILE_OUT.tmp > $FILE_OUT
        rm $FILE_OUT.tmp 
done 
