
### Command to get number of genes in loci
FILE_IN=/projects/timshel/sc-genetics/sc-genetics/src/depict/results/BMI_UKBB_Loh2018_no_mhc.5e-8.mousebrain_sem_mean_loci.txt
cat $FILE_IN | perl -lane '@T=split(";", $F[6]); print(join("\n",@T))' | sort -u | wc -l
1693 ---> subtract 1 ("genes_in_locus") --> 1692 genes