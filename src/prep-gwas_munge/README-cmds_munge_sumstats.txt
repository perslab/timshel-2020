
###################################### SETUP ######################################

### Activate
conda activate ldsc_munge

### Check
python -V
---> Python 2.7.16 :: Anaconda, Inc.

### Notes
The environment was created from this file:
CELLECT/ldsc/environment_munge.yml


###################################### NEW 2020 [with ldsc_munge env] ######################################

GWAS=HIPPOVOL_Hilbar2017
python /projects/timshel/sc-genetics/CELLECT/ldsc/munge_sumstats.py \
--sumstats /nfsdata/projects/timshel/data/gwas_sumstats_raw/HIPPOVOL_Hilbar2017/CHARGE-ENIGMA-HV-METAANALYSIS-201311141.TBL.FINAL.gz \
--merge-alleles /projects/timshel/sc-genetics/CELLECT/data/ldsc/w_hm3.snplist \
--ignore Weight,MarkerName \
--snp RSNUMBERS \
--N-col N \
--out /projects/timshel/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


###################################### CMDS ######################################

GWAS=BMI_Yengo2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/BMI_HEIGHT_Yengo2018/Meta-analysis_Locke_et_al+UKBiobank_2018.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--a1 Tested_Allele \
--a2 Other_Allele \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=BMI_UPDATE_Yengo2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/BMI_HEIGHT_Yengo2018/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--a1 Tested_Allele \
--a2 Other_Allele \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=HEIGHT_Yengo2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/BMI_HEIGHT_Yengo2018/Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--a1 Tested_Allele \
--a2 Other_Allele \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=EA3_Lee2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/EA3_Lee2018/GWAS_EduYears_excl23andMe.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 766345 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=SCZ_EUR_Ripke2014
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/SCZ_Ripke2014/daner_PGC_SCZ49.sh2_mds10_1000G-frq_2.OR-FILTER.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 77096 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=ASD_iPSYCH_PGC_Grove2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/ASD_iPSYCH_PGC_Grove2018/iPSYCH-PGC_ASD_Nov2017.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 46350 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=WHRadjBMI_Shungin2015
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/WHR_Shungin2015/GIANT_2015_WHRadjBMI_COMBINED_EUR.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=WHR_Shungin2015
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/WHR_Shungin2015/GIANT_2015_WHR_COMBINED_EUR.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=AN_PGC_Duncan2017
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/AN_PGC_Duncan2017/pgc.ed.freeze1.summarystatistics.July2017.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 14477 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=ADHD_PGC_Demontis2017
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/ADHD_PGC_Demontis2017/adhd_eur_jun2017.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 53293 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=MDD_PGC_Wray2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/MDD_PGC_Wray2018/MDD2018_ex23andMe.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 153780 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=INSOMNIA_Jansen2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/INSOMNIA_Jansen2018/Insomnia_sumstats_Jansenetal.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 386533 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=BIP_PGC2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/BIP_SCZ_PGC2018/BDvsCONT.sumstats.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 41653 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=SCZ_Pardinas2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/SCZ_Pardinas2018/clozuk_pgc2.meta.sumstats.rs_snps_only.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 105318 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=EA2_Okbay2016
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/EA2_Okbay2016/EduYears_Main.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 217568 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=AD_Jansen2019
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/AD_Jansen2019/AD_sumstats_Jansenetal.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--signed-sumstats Z,0 \
--N-col Nsum \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=AD_Jansen2019
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/AD_Jansen2019/AD_sumstats_Jansenetal.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--signed-sumstats Z,0 \
--N-col Nsum \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=NEUROTICISM_Nagel2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/NEUROTICISM_Nagel2018/sumstats_neuroticism_ctg_format.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--ignore SNP \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=DEPRESSION_Nagel2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/NEUROTICISM_Nagel2018/sumstats_depression_ctg_format.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--ignore SNP \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=DEPRESSED_AFFECT_Nagel2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/NEUROTICISM_Nagel2018/sumstats_depressed_affect_ctg_format.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--signed-sumstats Z,0 \
--ignore SNP,MAF \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=WORRY_Nagel2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/NEUROTICISM_Nagel2018/sumstats_worry_ctg_format.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--signed-sumstats Z,0 \
--ignore SNP,MAF \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=WORRY_Nagel2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/NEUROTICISM_Nagel2018/sumstats_worry_ctg_format.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--signed-sumstats Z,0 \
--ignore SNP,MAF \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=INTELLIGENCE_Sniekers2017
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/INTELLIGENCE_Sniekers2017/sumstats.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--a1 ref \
--a2 alt \
--ignore Zscore \
--N 78308 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=INTELLIGENCE_Savage2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/INTELLIGENCE_Savage2018/SavageJansen_2018_intelligence_metaanalysis.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N-col N_analyzed \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=RB_Linner_2019
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/RB_Linner_2019/RISK_GWAS_MA_UKB+replication.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 466571 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=BMI_Locke2015
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/BMI_Locke2015/EUR_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=HEIGHT_Wood2014
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/HEIGHT_Wood2014/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=T2DadjBMI_DIAMANTE_Mahajan2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/T2D_DIAMANTE_Mahajan2018/Mahajan.NatGenet2018b.T2Dbmiadj.European.w_snpsnap_rsid.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N-col Neff \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}
# Number of SNPs in GWAS data: 21635866
# Number of SNPs in GWAS data *NOT found* in SNP chr pos mapping file: 13965922
# Percent SNPs not found: 64.55 %

GWAS=T2D_DIAMANTE_Mahajan2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/T2D_DIAMANTE_Mahajan2018/Mahajan.NatGenet2018b.T2D.European.w_snpsnap_rsid.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N-col Neff \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}
# Number of SNPs in GWAS data: 23465132
# Number of SNPs in GWAS data *NOT found* in SNP chr pos mapping file: 15792766
# Percent SNPs not found: 67.30 %

GWAS=T2D_UKBB_DIAMANTE_Mahajan2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/T2D_DIAMANTE_Mahajan2018/Mahajan.NatGenet2018b.UKBB.HRC.T2D.European.w_snpsnap_rsid.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}
# Number of SNPs in GWAS data: 13583103
# Number of SNPs in GWAS data *NOT found* in SNP chr pos mapping file: 5931856
# Percent SNPs not found: 43.67 %


GWAS=T2D_Xue2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/T2D_Xue2018/Xue_et_al_T2D_META_Nat_Commun_2018.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=HBA1C_MAGIC_Wheeler2017
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/HBA1C_MAGIC_Wheeler2017/HbA1c_METAL_European.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 123665 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=FG_Male_Lagou2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/FI_FG_Lagou2018/FG_STAGE1_2_3_SEX_GWAS_2018.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--signed-sumstats male_beta,0 \
--p male_pvalue \
--N 67506 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=FG_Female_Lagou2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/FI_FG_Lagou2018/FG_STAGE1_2_3_SEX_GWAS_2018.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--signed-sumstats female_beta,0 \
--p female_pvalue \
--N 73089 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=FI_Male_Lagou2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/FI_FG_Lagou2018/FI_STAGE1_2_3_SEX_GWAS_2018.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--signed-sumstats male_beta,0 \
--p male_pvalue \
--N 47806 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=FI_Female_Lagou2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/FI_FG_Lagou2018/FI_STAGE1_2_3_SEX_GWAS_2018.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--signed-sumstats female_beta,0 \
--p female_pvalue \
--N 50404 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=HBA1C_MAGIC_Wheeler2017
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/HBA1C_MAGIC_Wheeler2017/HbA1c_METAL_European.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 159940 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=WHRadjBMI_Pulit2019
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/WHR_Pulit2019/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.rs_snps_clean.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--a1 Tested_Allele \
--a2 Other_Allele \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=WHR_Pulit2019
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/WHR_Pulit2019/whr.giant-ukbb.meta-analysis.combined.23May2018.rs_snps_clean.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--a1 Tested_Allele \
--a2 Other_Allele \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=BMI_Pulit2019
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/WHR_Pulit2019/bmi.giant-ukbb.meta-analysis.combined.23May2018.rs_snps_clean.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--a1 Tested_Allele \
--a2 Other_Allele \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=BMI_female_Pulit2019
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/WHR_Pulit2019/bmi.giant-ukbb.meta-analysis.females.23May2018.rs_snps_clean.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--a1 Tested_Allele \
--a2 Other_Allele \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=BMI_male_Pulit2019
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/WHR_Pulit2019/bmi.giant-ukbb.meta-analysis.males.23May2018.rs_snps_clean.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--a1 Tested_Allele \
--a2 Other_Allele \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=MDD_Howard2019
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/MDD_Howard2019/PGC_UKB_depression_genome-wide.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 500199 \
--signed-sumstats LogOR,0 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=AD_Lambert2013
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/AD_Lambert2013/IGAP_stage_1.txt \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 54162 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=RA_Okada2014
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/RA_Okada2014/RA_GWASmeta_European_v2.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 58284 \
--signed-sumstats "OR(A1)",1 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=DS_Okbay2016
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/SWB_Neuroticism_DS_Okbay2016/DS_Full.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 161460 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=NEUROTICISM_Okbay2016
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/SWB_Neuroticism_DS_Okbay2016/Neuroticism_Full.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 170911 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=SWB_Okbay2016
# *OBS*: N is subtracted 23andMe cohort size
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/SWB_Neuroticism_DS_Okbay2016/SWB_Full.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 204966 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=LIPIDS_HDL_Willer2013
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/LIPIDS_Willer2013/jointGwasMc_HDL.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=LIPIDS_LDL_Willer2013
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/LIPIDS_Willer2013/jointGwasMc_LDL.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=LIPIDS_TC_Willer2013
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/LIPIDS_Willer2013/jointGwasMc_TC.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=LIPIDS_TG_Willer2013
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/LIPIDS_Willer2013/jointGwasMc_TG.txt.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}



GWAS=LIPIDS_HDL_Teslovich2010
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/LIPIDS_Teslovich2010/HDL_ONE_Europeans.tbl.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N-col Weight \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=LIPIDS_LDL_Teslovich2010
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/LIPIDS_Teslovich2010/LDL_ONE_Europeans.tbl.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N-col Weight \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=LIPIDS_TG_Teslovich2010
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/LIPIDS_Teslovich2010/TG_ONE_Europeans.tbl.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N-col Weight \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=LIPIDS_TC_Teslovich2010
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/LIPIDS_Teslovich2010/TC_ONE_Europeans.tbl.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N-col Weight \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}




###################################### alkesgroup_UKBB/UKBB_460k ######################################
### OBS: Neff used for N

GWAS=BMI_UKBB_Loh2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/alkesgroup_UKBB/UKBB_460k/body_BMIz.sumstats.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 532396 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=HEIGHT_UKBB_Loh2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/alkesgroup_UKBB/UKBB_460k/body_HEIGHTz.sumstats.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 673879 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=WHRadjBMI_UKBB_Loh2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/alkesgroup_UKBB/UKBB_460k/body_WHRadjBMIz.sumstats.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 502773 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=DIASTOLICadjMED_UKBB_Loh2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/alkesgroup_UKBB/UKBB_460k/bp_DIASTOLICadjMEDz.sumstats.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 490470 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=SYSTOLICadjMED_UKBB_Loh2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/alkesgroup_UKBB/UKBB_460k/bp_SYSTOLICadjMEDz.sumstats.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 469767 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}


GWAS=CARDIOVASCULAR_UKBB_Loh2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/alkesgroup_UKBB/UKBB_460k/disease_CARDIOVASCULAR.sumstats.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 477807 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}

GWAS=T2D_UKBB_Loh2018
python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
--sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/alkesgroup_UKBB/UKBB_460k/disease_T2D.sumstats.gz \
--merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
--N 468298 \
--out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}



###################################### LOOPS ######################################

### BLOOD COUNTS
for GWAS in blood_EOSINOPHIL_COUNT blood_PLATELET_COUNT; do
  python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
  --sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/alkesgroup_UKBB/${GWAS}.sumstats.gz \
  --merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
  --out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS} &
done
# blood_PLATELET_COUNT ---> ValueError: Columns ['P'] are expected to be numeric


### NULL GWAS
for i in {1..10}; do
  GWAS=1KG_phase3_EUR_null_gwas_P${i}
  echo $GWAS
  python2 /projects/timshel/sc-genetics/ldsc/ldsc/munge_sumstats.py \
  --sumstats /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/NULL_GWAS/plink_out.1KG_phase3_EUR_null_gwas.P${i}.qassoc_with_alleles \
  --merge-alleles /projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist \
  --N 503 \
  --out /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS} &
done


###################################### SYMLINK ALKES sumstats_formatted ######################################
### TRAITS I DON'T HAVE

PATH_SOURCE=/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/alkesgroup/sumstats_formatted
PATH_DEST=/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection

GWAS=T1D_Bradfield2011
ln -s $PATH_SOURCE/Type_1_Diabetes.sumstats.gz $PATH_DEST/${GWAS}.sumstats.gz
GWAS=CELIAC_Dubois2010
ln -s $PATH_SOURCE/Celiac.sumstats.gz $PATH_DEST/${GWAS}.sumstats.gz
GWAS=UC_Jostins2012
ln -s $PATH_SOURCE/Ulcerative_Colitis.sumstats.gz $PATH_DEST/${GWAS}.sumstats.gz
GWAS=CROHNS_Jostins2012
ln -s $PATH_SOURCE/Crohns_Disease.sumstats.gz $PATH_DEST/${GWAS}.sumstats.gz
GWAS=IBD_Jostins2012
ln -s $PATH_SOURCE/IBD.sumstats.gz $PATH_DEST/${GWAS}.sumstats.gz
GWAS=PBC_Cordell2015
ln -s $PATH_SOURCE/Primary_biliary_cirrhosis.sumstats.gz $PATH_DEST/${GWAS}.sumstats.gz
GWAS=LUPUS_Bentham2015
ln -s $PATH_SOURCE/Lupus.sumstats.gz $PATH_DEST/${GWAS}.sumstats.gz
GWAS=CAD_Schunkert2011
ln -s $PATH_SOURCE/Coronary_Artery_Disease.sumstats.gz $PATH_DEST/${GWAS}.sumstats.gz
GWAS=MS_Patsopoulos2011
ln -s $PATH_SOURCE/Multiple_sclerosis.sumstats.gz $PATH_DEST/${GWAS}.sumstats.gz



