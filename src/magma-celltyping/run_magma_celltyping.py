#!/usr/bin/env python2.7

import os
import sys
import time

import glob

import subprocess

###################################### CMDS manual ######################################


# time Rscript magma_celltyping.R EA2_Okbay2016 &> log.tmp.magma_celltyping.EA2_Okbay2016.out.txt
# time Rscript magma_celltyping.R SCZ_Ripke2014 &> log.tmp.magma_celltyping.SCZ_Ripke2014.out.txt
# time Rscript magma_celltyping.R HEIGHT_Wood2014 &> log.tmp.magma_celltyping.HEIGHT_Wood2014.out.txt


###################################### MAIN ######################################


### GWAS
list_gwas = ["BMI_Yengo2018",
"EA3_Lee2018",
"SCZ_Ripke2014",
"HEIGHT_Yengo2018",
"blood_EOSINOPHIL_COUNT"] 

### GWAS
# list_gwas = ["ADHD_PGC_Demontis2017",
# "AN_PGC_Duncan2017",
# "ASD_iPSYCH_PGC_Anney2017",
# "ASD_iPSYCH_PGC_Grove2018",
# "BIP_SCZ_BDvsCONT_PGC2018",
# "BIP_SCZ_SCZvsBD_PGC2018",
# "BMI_Yengo2018",
# "HEIGHT_Yengo2018",
# "EA2_Okbay2016",
# "EA3_Lee2018",
# "FG_male_Lagou2018",
# "FI_male_Lagou2018",
# "HEIGHT_Wood2014",
# "INSOMNIA_Hammerschlag2017",
# "LIPIDS_HDL_Willer2013",
# "LIPIDS_LDL_Willer2013",
# "LIPIDS_TC_Willer2013",
# "MDD_PGC_Wray2018",
# "RA_Okada2014",
# "SCZ_Pardinas2018",
# "SCZ_Ripke2014",
# "T2D_70kforT2D_Guarch2018",
# "WHR_adjBMI_Shungin2015",
# "WHR_Shungin2015",
# # "Daytime_dozing_sleeping_UKBB2018", # Rosa
# # "Nap_during_day_UKBB2018", # Rosa
# "blood_EOSINOPHIL_COUNT"] 

### NULL GWAS (1...10)
# list_gwas = ["1KG_phase3_EUR_null_gwas_P{}".format(i) for i in range(1,11)] # ['1KG_phase3_EUR_null_gwas_P1', '1KG_phase3_EUR_null_gwas_P2', '1KG_phase3_EUR_null_gwas_P3', '1KG_phase3_EUR_null_gwas_P4', '1KG_phase3_EUR_null_gwas_P5', '1KG_phase3_EUR_null_gwas_P6', '1KG_phase3_EUR_null_gwas_P7', '1KG_phase3_EUR_null_gwas_P8', '1KG_phase3_EUR_null_gwas_P9', '1KG_phase3_EUR_null_gwas_P10']


### LDSC gwas
# list_gwas = ["AgeFirstBirth",
# "Alzheimer",
# "Anorexia",
# "Autism",
# "Bipolar_Disorder",
# "BMI1",
# "Celiac",
# "Coronary_Artery_Disease",
# "Crohns_Disease",
# "DS",
# "ENIGMA2_MeanPutamen",
# "Ever_Smoked",
# "FastingGlucose_Manning",
# "HbA1C",
# "HDL",
# "Height1",
# "IBD",
# "LDL",
# "Lupus",
# "Menarche2017",
# "Neuroticism",
# "NumberChildrenEverBorn",
# "Primary_biliary_cirrhosis",
# "Rheumatoid_Arthritis",
# "Schizophrenia",
# "SWB",
# "Triglycerides",
# "Type_2_Diabetes",
# "Ulcerative_Colitis",
# "Years_of_Education1",
# "Years_of_Education2"]


# list_cmds = ["Rscript magma_celltyping.R {gwas_name} &> log.tmp.magma_celltyping.{gwas_name}.out.txt".format(gwas_name=gwas_name) for gwas_name in list_gwas] # standard cell-typing
list_cmds = ["Rscript magma_celltyping-non_quantiles.R {gwas_name} &> log.tmp.magma_celltyping-non-quantiles.{gwas_name}.out.txt".format(gwas_name=gwas_name) for gwas_name in list_gwas]
# list_cmds = ["sleep 10"] * 10 # ['sleep 3', 'sleep 3', 'sleep 3']

BATCH_SIZE = 25 # number of jobs to run parallel

list_of_processes = []
batch = 1
for i, cmd in enumerate(list_cmds, start=1):
	print "batch = {} | i = {} | Running command: {}".format(batch, i, cmd)
	print cmd
	p = subprocess.Popen(cmd, shell=True)
	list_of_processes.append(p)
	print "batch = {} | i = {} | list_of_processes:".format(batch, i)
	print list_of_processes
	if i % BATCH_SIZE == 0:
		batch += 1
		for process in list_of_processes:
			print "=========== Waiting for process: {} ===========".format(process.pid)
			process.wait()
		list_of_processes = [] # 'reset' list

### wait for the rest for the rest of the processes
print "*********** DONE IN MAIN LOOP*************."
for process in list_of_processes:
	print "=========== Waiting for process: {} ===========".format(process.pid)
	process.wait()
			


print "SCRIPT ENDED"

