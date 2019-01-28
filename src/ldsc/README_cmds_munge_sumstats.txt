
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

### BMI_UPDATE_Yengo2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/BMI_HEIGHT_Yengo2018/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--a1 Tested_Allele \
--a2 Other_Allele \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/BMI_UPDATE_Yengo2018

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


### SCZ_EUR_Ripke2014 (EUR ONLY)
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/SCZ_Ripke2014/daner_PGC_SCZ49.sh2_mds10_1000G-frq_2.OR-FILTER.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 77096 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/SCZ_EUR_Ripke2014


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

GWAS=BIP_PGC2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/BIP_SCZ_PGC2018/BDvsCONT.sumstats.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 41653 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=SCZ_Pardinas2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/SCZ_Pardinas2018/clozuk_pgc2.meta.sumstats.rs_snps_only.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 105318 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=EA2_Okbay2016
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/EA2_Okbay2016/EduYears_Main.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 217568 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=AD_Jansen2019
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/AD_Jansen2019/AD_sumstats_Jansenetal.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--signed-sumstats Z,0 \
--N-col Nsum \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=AD_Jansen2019
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/AD_Jansen2019/AD_sumstats_Jansenetal.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--signed-sumstats Z,0 \
--N-col Nsum \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=NEUROTICISM_Nagel2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/NEUROTICISM_Nagel2018/sumstats_neuroticism_ctg_format.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--ignore SNP \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=DEPRESSION_Nagel2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/NEUROTICISM_Nagel2018/sumstats_depression_ctg_format.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--ignore SNP \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=DEPRESSED_AFFECT_Nagel2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/NEUROTICISM_Nagel2018/sumstats_depressed_affect_ctg_format.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--signed-sumstats Z,0 \
--ignore SNP,MAF \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=WORRY_Nagel2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/NEUROTICISM_Nagel2018/sumstats_worry_ctg_format.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--signed-sumstats Z,0 \
--ignore SNP,MAF \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=WORRY_Nagel2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/NEUROTICISM_Nagel2018/sumstats_worry_ctg_format.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--signed-sumstats Z,0 \
--ignore SNP,MAF \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=INTELLIGENCE_Sniekers2017
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/INTELLIGENCE_Sniekers2017/sumstats.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--a1 ref \
--a2 alt \
--ignore Zscore \
--N 78308 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=INTELLIGENCE_Savage2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/INTELLIGENCE_Savage2018/SavageJansen_2018_intelligence_metaanalysis.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N-col N_analyzed \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=RB_Linner_2019
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/RB_Linner_2019/RISK_GWAS_MA_UKB+replication.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 466571 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=BMI_Locke2015
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/BMI_Locke2015/EUR_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=HEIGHT_Wood2014
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/HEIGHT_Wood2014/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=T2DadjBMI_DIAMANTE_Mahajan2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/T2D_DIAMANTE_Mahajan2018/Mahajan.NatGenet2018b.T2Dbmiadj.European.w_snpsnap_rsid.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N-col Neff \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}
# Number of SNPs in GWAS data: 21635866
# Number of SNPs in GWAS data *NOT found* in SNP chr pos mapping file: 13965922
# Percent SNPs not found: 64.55 %

GWAS=T2D_DIAMANTE_Mahajan2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/T2D_DIAMANTE_Mahajan2018/Mahajan.NatGenet2018b.T2D.European.w_snpsnap_rsid.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N-col Neff \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}
# Number of SNPs in GWAS data: 23465132
# Number of SNPs in GWAS data *NOT found* in SNP chr pos mapping file: 15792766
# Percent SNPs not found: 67.30 %

GWAS=T2D_UKBB_DIAMANTE_Mahajan2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/T2D_DIAMANTE_Mahajan2018/Mahajan.NatGenet2018b.UKBB.HRC.T2D.European.w_snpsnap_rsid.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}
# Number of SNPs in GWAS data: 13583103
# Number of SNPs in GWAS data *NOT found* in SNP chr pos mapping file: 5931856
# Percent SNPs not found: 43.67 %


GWAS=T2D_Xue2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/T2D_Xue2018/Xue_et_al_T2D_META_Nat_Commun_2018.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=HBA1C_MAGIC_Wheeler2017
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/HBA1C_MAGIC_Wheeler2017/HbA1c_METAL_European.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 123665 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=FG_Male_Lagou2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/FI_FG_Lagou2018/FG_STAGE1_2_3_SEX_GWAS_2018.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--signed-sumstats male_beta,0 \
--p male_pvalue \
--N 67506 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=FG_Female_Lagou2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/FI_FG_Lagou2018/FG_STAGE1_2_3_SEX_GWAS_2018.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--signed-sumstats female_beta,0 \
--p female_pvalue \
--N 73089 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=FI_Male_Lagou2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/FI_FG_Lagou2018/FI_STAGE1_2_3_SEX_GWAS_2018.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--signed-sumstats male_beta,0 \
--p male_pvalue \
--N 47806 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=FI_Female_Lagou2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/FI_FG_Lagou2018/FI_STAGE1_2_3_SEX_GWAS_2018.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--signed-sumstats female_beta,0 \
--p female_pvalue \
--N 50404 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=HBA1C_MAGIC_Wheeler2017
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/HBA1C_MAGIC_Wheeler2017/HbA1c_METAL_European.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 159940 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=WHRadjBMI_Pulit2019
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/WHR_Pulit2019/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.rs_snps_clean.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--a1 Tested_Allele \
--a2 Other_Allele \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=WHR_Pulit2019
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/WHR_Pulit2019/whr.giant-ukbb.meta-analysis.combined.23May2018.rs_snps_clean.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--a1 Tested_Allele \
--a2 Other_Allele \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=BMI_Pulit2019
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/WHR_Pulit2019/bmi.giant-ukbb.meta-analysis.combined.23May2018.rs_snps_clean.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--a1 Tested_Allele \
--a2 Other_Allele \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=BMI_female_Pulit2019
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/WHR_Pulit2019/bmi.giant-ukbb.meta-analysis.females.23May2018.rs_snps_clean.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--a1 Tested_Allele \
--a2 Other_Allele \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=BMI_male_Pulit2019
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/WHR_Pulit2019/bmi.giant-ukbb.meta-analysis.males.23May2018.rs_snps_clean.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--a1 Tested_Allele \
--a2 Other_Allele \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=MDD_Howard2019
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/MDD_Howard2019/PGC_UKB_depression_genome-wide.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 500199 \
--signed-sumstats LogOR,0 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

###################################### alkesgroup_UKBB/UKBB_460k ######################################
### OBS: Neff used for N


GWAS=BMI_UKBB_Loh2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/alkesgroup_UKBB/UKBB_460k/body_BMIz.sumstats.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 532396 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=HEIGHT_UKBB_Loh2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/alkesgroup_UKBB/UKBB_460k/body_HEIGHTz.sumstats.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 673879 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}



GWAS=WHRadjBMI_UKBB_Loh2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/alkesgroup_UKBB/UKBB_460k/body_WHRadjBMIz.sumstats.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 502773 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=DIASTOLICadjMED_UKBB_Loh2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/alkesgroup_UKBB/UKBB_460k/bp_DIASTOLICadjMEDz.sumstats.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 490470 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=SYSTOLICadjMED_UKBB_Loh2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/alkesgroup_UKBB/UKBB_460k/bp_SYSTOLICadjMEDz.sumstats.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 469767 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=CARDIOVASCULAR_UKBB_Loh2018
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/alkesgroup_UKBB/UKBB_460k/disease_CARDIOVASCULAR.sumstats.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 477807 \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=LIPIDS_HDL_Willer2013
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/LIPIDS_Willer2013/jointGwasMc_HDL.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=LIPIDS_LDL_Willer2013
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/LIPIDS_Willer2013/jointGwasMc_LDL.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=LIPIDS_TC_Willer2013
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/LIPIDS_Willer2013/jointGwasMc_TC.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=LIPIDS_TG_Willer2013
python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/LIPIDS_Willer2013/jointGwasMc_TG.txt.gz \
--merge-alleles /raid5/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--out /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}



###################################### SYMLINK ALKES sumstats_formatted ######################################

PATH_SOURCE=/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/alkesgroup/sumstats_formatted
PATH_DEST=/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection
ln -s $PATH_SOURCE/Type_1_Diabetes.sumstats.gz $PATH_DEST/T1D_Bradfield2011.sumstats.gz
ln -s $PATH_SOURCE/HDL.sumstats.gz $PATH_DEST/LIPIDS_HDL_Teslovich2010.sumstats.gz
ln -s $PATH_SOURCE/LDL.sumstats.gz $PATH_DEST/LIPIDS_LDL_Teslovich2010.sumstats.gz
ln -s $PATH_SOURCE/Triglycerides.sumstats.gz $PATH_DEST/LIPIDS_TG_Teslovich2010.sumstats.gz
ln -s $PATH_SOURCE/Celiac.sumstats.gz $PATH_DEST/CELIAC_Dubois2010.sumstats.gz
ln -s $PATH_SOURCE/Ulcerative_Colitis.sumstats.gz $PATH_DEST/UC_Jostins2012.sumstats.gz
ln -s $PATH_SOURCE/Crohns_Disease.sumstats.gz $PATH_DEST/CROHNS_Jostins2012.sumstats.gz
ln -s $PATH_SOURCE/IBD.sumstats.gz $PATH_DEST/IBD_Jostins2012.sumstats.gz
ln -s $PATH_SOURCE/Primary_biliary_cirrhosis.sumstats.gz $PATH_DEST/PBC_Cordell2015.sumstats.gz
ln -s $PATH_SOURCE/Lupus.sumstats.gz $PATH_DEST/LUPUS_2015.sumstats.gz
ln -s $PATH_SOURCE/Rheumatoid_Arthritis.sumstats.gz $PATH_DEST/RA_Okada2014.sumstats.gz
ln -s $PATH_SOURCE/Coronary_Artery_Disease.sumstats.gz $PATH_DEST/CAD_Schunkert2011.sumstats.gz

ln -s $PATH_SOURCE/Alzheimer.sumstats.gz $PATH_DEST/AD_Lambert2013.sumstats.gz
ln -s $PATH_SOURCE/Multiple_sclerosis.sumstats.gz $PATH_DEST/MS_Patsopoulos2011.sumstats.gz

ln -s $PATH_SOURCE/SWB.sumstats.gz $PATH_DEST/DS_Okbay2016.sumstats.gz
ln -s $PATH_SOURCE/DS.sumstats.gz $PATH_DEST/NEUROTICISM_OKBAY2016.sumstats.gz
ln -s $PATH_SOURCE/Neuroticism.sumstats.gz $PATH_DEST/SWB_Okbay2016.sumstats.gz


###################################### DOCS munge_sumstats.py ######################################

--no-alleles: Don't require alleles. Useful if only unsigned summary statistics are available and the goal is h2 /partitioned h2 estimation rather than rg estimation
--no-alleles and --merge-alleles are not compatible.

--a1: A1 -- first allele (effect allele)
--a2: A2-- second allele (other allele)


  --signed-sumstats SIGNED_SUMSTATS : Name of signed sumstat column, comma null value (e.g., Z,0 or OR,1). NB: case insensitive.

