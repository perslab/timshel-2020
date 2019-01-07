
###################################### munge_sumstats.py must be run with pandas 0.20.1? [at least for INSOMNIA_Jansen2018] ######################################
----> SOLUTION [I had pandas 0.23.4 before]: conda install pandas=0.20.1

REF: https://groups.google.com/forum/#!msg/ldsc_users/9Mxw4vAB2l8/Muc_vtmYAwAJ
Hi Max,
Yes, pandas versioning is almost certainly the issue here. My understanding is this error occurs with pandas 0.21, but works correctly if you revert to e.g. 0.20.1.
Cheers,
Raymond

@@@@@@@@@@@@@@@@@@ Error when using pandas 0.23.4 @@@@@@@@@@@@@@@@@@
Reading list of SNPs for allele merge from /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist
Read 1217311 SNPs for allele merge.
Reading sumstats from /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/INSOMNIA_Jansen2018/Insomnia_sumstats_Jansenetal.txt.gz into memory 5000000 SNPs at a time.
... done
Read 10862567 SNPs from --sumstats file.
Removed 9698237 SNPs not in --merge-alleles.
Removed 0 SNPs with missing values.
Removed 0 SNPs with INFO <= 0.9.
Removed 27921 SNPs with MAF <= 0.01.
Removed 0 SNPs with out-of-bounds p-values.
Removed 87 variants that were not SNPs or were strand-ambiguous.
1136322 SNPs remain.
Removed 0 SNPs with duplicated rs numbers (1136322 SNPs remain).
Removed 0 SNPs with N < 257688.666667 (1136322 SNPs remain).
Median value of OR was 1.0, which seems sensible.
Removed 120 SNPs whose alleles did not match --merge-alleles (1136202 SNPs remain).

ERROR converting summary statistics:

Traceback (most recent call last):
  File "/raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py", line 707, in munge_sumstats
    dat = allele_merge(dat, merge_alleles, log)
  File "/raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py", line 445, in allele_merge
    dat.loc[~jj, [i for i in dat.columns if i != 'SNP']] = float('nan')
  File "/tools/anaconda/3-4.4.0/envs/py27_anaconda3_PT170705/lib/python2.7/site-packages/pandas/core/indexing.py", line 188, in __setitem__
    indexer = self._get_setitem_indexer(key)
  File "/tools/anaconda/3-4.4.0/envs/py27_anaconda3_PT170705/lib/python2.7/site-packages/pandas/core/indexing.py", line 166, in _get_setitem_indexer
    return self._convert_tuple(key, is_setter=True)
  File "/tools/anaconda/3-4.4.0/envs/py27_anaconda3_PT170705/lib/python2.7/site-packages/pandas/core/indexing.py", line 247, in _convert_tuple
    idx = self._convert_to_indexer(k, axis=i, is_setter=is_setter)
  File "/tools/anaconda/3-4.4.0/envs/py27_anaconda3_PT170705/lib/python2.7/site-packages/pandas/core/indexing.py", line 1327, in _convert_to_indexer
    .format(mask=objarr[mask]))
KeyError: '[-2 -2 -2 ... -1 -1 -2] not in index'

###################################### CMDS ######################################

### BMI_Yengo2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/BMI_HEIGHT_Yengo2018/Meta-analysis_Locke_et_al+UKBiobank_2018.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--a1 Tested_Allele \
--a2 Other_Allele \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/BMI_Yengo2018


### HEIGHT_Yengo2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/BMI_HEIGHT_Yengo2018/Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--a1 Tested_Allele \
--a2 Other_Allele \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/HEIGHT_Yengo2018

### EA3_Lee2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/EA3_Lee2018/GWAS_EduYears_excl23andMe.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 766345 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/EA3_Lee2018


### SCZ_Ripke2014
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/SCZ_Ripke2014/ckqny.scz2snpres.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 150064 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/SCZ_Ripke2014


### LIPIDS_HDL_Willer2013
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/LIPIDS_Willer2013/jointGwasMc_HDL.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/LIPIDS_HDL_Willer2013


### RA_Okada2014
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/RA_Okada2014/RA_GWASmeta_European_v2.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 58284 \
--signed-sumstats "OR(A1),1" \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/RA_Okada2014


### ASD_iPSYCH_PGC_Grove2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/ASD_iPSYCH_PGC_Grove2018/iPSYCH-PGC_ASD_Nov2017.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 46350 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/ASD_iPSYCH_PGC_Grove2018

################## NEW ##################

### WHR_adjBMI_Shungin2015
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/WHR_Shungin2015/GIANT_2015_WHRadjBMI_COMBINED_EUR.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/WHR_adjBMI_Shungin2015

### WHR_Shungin2015
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/WHR_Shungin2015/GIANT_2015_WHR_COMBINED_EUR.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/WHR_Shungin2015


###
GWAS=AN_PGC_Duncan2017
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/AN_PGC_Duncan2017/pgc.ed.freeze1.summarystatistics.July2017.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 14477 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


### 
GWAS=ADHD_PGC_Demontis2017
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/ADHD_PGC_Demontis2017/adhd_eur_jun2017.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 53293 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


### 
GWAS=MDD_PGC_Wray2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/MDD_PGC_Wray2018/MDD2018_ex23andMe.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 153780 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


### BLOOD COUNTS
for GWAS in blood_EOSINOPHIL_COUNT blood_PLATELET_COUNT; do
	python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
	--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/alkesgroup_UKBB/${GWAS}.sumstats.gz \
	--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
	--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS} &
done
# blood_PLATELET_COUNT ---> ValueError: Columns ['P'] are expected to be numeric


### NULL GWAS
for i in {1..10}; do
	GWAS=1KG_phase3_EUR_null_gwas_P${i}
	echo $GWAS
	python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
	--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/NULL_GWAS/plink_out.1KG_phase3_EUR_null_gwas.P${i}.qassoc_with_alleles \
	--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
	--N 503 \
	--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS} &
done


###################################### NEW 2019 ######################################

GWAS=INSOMNIA_Jansen2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/INSOMNIA_Jansen2018/Insomnia_sumstats_Jansenetal.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 386533 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}





###################################### DOCS munge_sumstats.py ######################################

--no-alleles: Don't require alleles. Useful if only unsigned summary statistics are available and the goal is h2 /partitioned h2 estimation rather than rg estimation
--no-alleles and --merge-alleles are not compatible.

--a1: A1 -- first allele (effect allele)
--a2: A2-- second allele (other allele)


  --signed-sumstats SIGNED_SUMSTATS : Name of signed sumstat column, comma null value (e.g., Z,0 or OR,1). NB: case insensitive.

