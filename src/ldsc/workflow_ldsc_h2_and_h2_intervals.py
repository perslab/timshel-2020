
from __future__ import print_function
import argparse
import glob
import os
import subprocess
import re
import sys


import collections


###################################### USAGE ######################################
# Compatibility: Python 2 and 3

### Run in unbuffered mode
# time python -u workflow_ldsc_cts.py |& tee workflow_ldsc_cts.UNNAMED.out.txt

###################################### DOCS ######################################


###################################### Pipeline description ######################################

# 1. run quantile_M.pl for your annotation. [Only need to do once for each annotation]
# 2. run LDSC --h2 with GWAS --> output jackknife and .results with coefficients 
# 3. quantile_h2g.r --> compute h2 based on quantiles
# 4. (parse results)


###################################### UTILS - ALL GENES ######################################


def get_all_genes_ref_ld_chr_name(dataset):
	""" Function to get the ref_ld_chr_name for 'all genes annotation' for ldsc.py --h2/--h2-cts command """
	# *IMPORTANT*: ldsc_all_genes_ref_ld_chr_name MUST be full file path PLUS trailing "."
	dict_dataset_all_genes_path_prefix = {"mousebrain":"/scratch/sc-ldsc/control.all_genes_in_dataset/per_annotation/control.all_genes_in_dataset__all_genes_in_dataset.mousebrain.",
										"tabula_muris":"/scratch/sc-ldsc/control.all_genes_in_dataset/per_annotation/control.all_genes_in_dataset__all_genes_in_dataset.tabula_muris.",
										"campbell":"/scratch/sc-ldsc/control.all_genes_in_dataset/per_annotation/control.all_genes_in_dataset__all_genes_in_dataset.campbell.",
										 }
	if not dataset in dict_dataset_all_genes_path_prefix:
		raise KeyError("dataset={} is not found in dict_dataset_all_genes_path_prefix.".format(dataset))
	ldsc_all_genes_ref_ld_chr_name = dict_dataset_all_genes_path_prefix[dataset]
	# some obnoxious validation of the matches
	files_ldscore = glob.glob("{}*l2.ldscore.gz".format(ldsc_all_genes_ref_ld_chr_name)) # get ldscore files for all chromosomes. glob() returns full file paths.
	if not len(files_ldscore) == 22: # we must have ldscore files for every chromosome, so the length 
		raise ValueError("dataset={} only has n={} matching {}*l2.ldscore.gz files. Expected 22 files. Check the ldscore file directory or update the dict_dataset_all_genes_path_prefix inside this function.".format(dataset, len(files_ldscore), ldsc_all_genes_ref_ld_chr_name))
	return(ldsc_all_genes_ref_ld_chr_name)

###################################### Job scheduler ######################################

def job_scheduler(list_cmds, n_parallel_jobs):
	""" Schedule parallel jobs with at most n_parallel_jobs parallel jobs."""
	list_of_processes = []
	batch = 1
	for i, cmd in enumerate(list_cmds, start=1):
		print("job schedule batch = {} | i = {} | Running command: {}".format(batch, i, cmd))
		## p = subprocess.Popen(cmd, shell=True, bufsize=0 if FLAG_UNBUFFERED else -1, stdout=FNULL, stderr=subprocess.STDOUT)
		### You need to keep devnull open for the entire life of the Popen object, not just its construction. 
		### FNULL = open(os.devnull, 'w') # devnull filehandle does not need to be closed?
		p = subprocess.Popen(cmd, shell=True, bufsize=0 if FLAG_UNBUFFERED else -1)
		list_of_processes.append(p)
		print("job schedule batch = {} | i = {} | PIDs of running jobs (list_of_processes):".format(batch, i))
		print(" ".join([str(p.pid) for p in list_of_processes])) # print PIDs
		if i % n_parallel_jobs == 0: # JOB BATCH SIZE
			batch += 1
			for p in list_of_processes:
				print("=========== Waiting for process: {} ===========".format(p.pid))
				sys.stdout.flush()
				p.wait()
				print("Returncode = {}".format(p.returncode))
			list_of_processes = [] # 'reset' list

	### wait for the rest for the rest of the processes
	for p in list_of_processes:
		print("=========== Waiting for process: {} ===========".format(p.pid))
		p.wait()

	return list_of_processes

##################################################################################################
############################################ CONSTANTS ###########################################
##################################################################################################



PYTHON3_EXEC = "/tools/anaconda/3-4.4.0/envs/py3_anaconda3_PT180510/bin/python3"
PYTHON2_EXEC = "/tools/anaconda/3-4.4.0/envs/py27_anaconda3_PT170705/bin/python2"

PATH_LDSC_SCRIPT = "/projects/timshel/sc-genetics/ldsc/ldsc-timshel/ldsc.py" 
PATH_LDSC_QUANTILE_PERL_SCRIPT="/projects/timshel/sc-genetics/ldsc/ldsc-timshel/ContinuousAnnotations/quantile_M_fixed_non_zero_quantiles.pl"
PATH_LDSC_H2_RSCRIPT_SCRIPT="/projects/timshel/sc-genetics/ldsc/ldsc-timshel/ContinuousAnnotations/quantile_h2g.r"
PATH_LDSC_DATA_MAIN="/projects/timshel/sc-genetics/ldsc/data"

LIST_CHROMOSOMES = list(range(1,23)) # 1..22 | only used to validate that all annot file exists

FLAG_UNBUFFERED = True

N_PARALLEL_QUANTILE_PERL_JOBS = 50
N_PARALLEL_LDSC_REGRESSION_JOBS = 10
N_PARALLEL_H2_RSCRIPT_JOBS = 50


### PRIM LIST GWAS [ALL RUN on h2 + h2_quantile]
# list_gwas = ["BMI_UKBB_Loh2018",
# 			 "BMI_UPDATE_Yengo2018",
# 			 "T2D_DIAMANTE_Mahajan2018",
# 			 "T2D_UKBB_DIAMANTE_Mahajan2018",
# 			 "T2D_UKBB_Loh2018",
# 			 "HEIGHT_UKBB_Loh2018",
# 			 "HEIGHT_Yengo2018",
# 			 "LIPIDS_LDL_Teslovich2010",
# 			 "RA_Okada2014",
# 			]


### Meta-analysis traits [n=39]
list_gwas = ["AD_Jansen2019",
"ADHD_PGC_Demontis2017",
"AN_PGC_Duncan2017",
"ASD_iPSYCH_PGC_Grove2018",
"BIP_PGC2018",
"BMI_UKBB_Loh2018",
"CAD_Schunkert2011",
"CARDIOVASCULAR_UKBB_Loh2018",
"CELIAC_Dubois2010",
"CROHNS_Jostins2012",
"DEPRESSION_Nagel2018",
"DIASTOLICadjMED_UKBB_Loh2018",
"EA3_Lee2018",
"FG_Female_Lagou2018",
"FI_Female_Lagou2018",
"HBA1C_MAGIC_Wheeler2017",
"HEIGHT_UKBB_Loh2018",
"IBD_Jostins2012",
"INSOMNIA_Jansen2018",
"INTELLIGENCE_Savage2018",
"LIPIDS_HDL_Teslovich2010",
"LIPIDS_LDL_Teslovich2010",
"LIPIDS_TG_Teslovich2010",
"LUPUS_Bentham2015",
"MDD_Howard2019",
"MS_Patsopoulos2011",
"NEUROTICISM_Nagel2018",
"PBC_Cordell2015",
"RA_Okada2014",
"RB_Linner_2019",
"SCZ_Pardinas2018",
"SWB_Okbay2016",
"SYSTOLICadjMED_UKBB_Loh2018",
"T1D_Bradfield2011",
"T2D_UKBB_Loh2018",
"UC_Jostins2012",
"WHR_Pulit2019",
"WHRadjBMI_UKBB_Loh2018",
"WORRY_Nagel2018"]

##################################################################################################
############################################ PARAMS ##############################################
##################################################################################################

FLAG_QUANTILE_FIXED = True # True is used for publication
# FLAG_QUANTILE_FIXED = False


################## MOUSEBRAIN ##################
dict_annotations = collections.defaultdict(dict)
list_annotations = ["TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5","MEGLU10","TEGLU17","MEGLU11","TEGLU4","DEGLU4","TEINH12"]
for annotation in list_annotations:
	dict_annotations[annotation]["name_context"] = "mousebrain_all.{}.sem_mean".format(annotation)
	dict_annotations[annotation]["file_path_prefix"] = "/scratch/sc-ldsc/celltypes.mousebrain.all/per_annotation/celltypes.mousebrain.all__{}".format(dict_annotations[annotation]["name_context"])
# ALT using dict comprehension: dict_annotations = {key:{"name_context":"mousebrain_all.{}.sem_mean".format(key)} for key in list_annotations}

dict_annotations_mb = dict_annotations


################## TABULA MURIS ##################
dict_annotations = collections.defaultdict(dict)
### LIST for h2_quantile
list_annotations = ["Brain_Non-Myeloid.neuron","Brain_Non-Myeloid.oligodendrocyte_precursor_cell","Brain_Non-Myeloid.oligodendrocyte"]
### FULL LIST for h2
# list_annotations = ["Brain_Non-Myeloid.neuron","Brain_Non-Myeloid.oligodendrocyte_precursor_cell","Brain_Non-Myeloid.oligodendrocyte","Brain_Non-Myeloid.astrocyte",
# "Pancreas.type_B_pancreatic_cell","Pancreas.pancreatic_A_cell","Pancreas.pancreatic_D_cell","Pancreas.pancreatic_PP_cell",
# "Liver.hepatocyte", # LDL
# "Spleen.T_cell", # RA
# "Trachea.mesenchymal_cell", # Height
# "Limb_Muscle.mesenchymal_stem_cell", # Height
# ]

for annotation in list_annotations:
	dict_annotations[annotation]["name_context"] = "tabula_muris.{}.sem_mean".format(annotation)
	dict_annotations[annotation]["file_path_prefix"] = "/scratch/sc-ldsc/celltypes.tabula_muris.all/per_annotation/celltypes.tabula_muris.all__{}".format(dict_annotations[annotation]["name_context"])
dict_annotations_tm = dict_annotations

################## WGCNA ##################
# /scratch/sc-ldsc/wgcna.mousebrain-190111.fdr_sign_celltypes.continuous/per_annotation/wgcna.mousebrain-190111.fdr_sign_celltypes.continuous__dodgerblue.7.l2.M
dict_annotations = collections.defaultdict(dict)
list_annotations = ["dodgerblue","darkgoldenrod3","wheat3"]
for annotation in list_annotations:
	dict_annotations[annotation]["name_context"] = "{}".format(annotation)
	dict_annotations[annotation]["file_path_prefix"] = "/scratch/sc-ldsc/wgcna.mousebrain-190111.fdr_sign_celltypes.continuous/per_annotation/wgcna.mousebrain-190111.fdr_sign_celltypes.continuous__{}".format(dict_annotations[annotation]["name_context"])
dict_annotations_wgcna = dict_annotations



################## COMBINE ##################
# Output name: <runname/key_prim>__<annotation_id>__<gwas>
dict_run = {"celltypes.mousebrain":
						{"dataset":"mousebrain",
						"dict_annotations":dict_annotations_mb},
			"celltypes.tabula_muris":
						{"dataset":"tabula_muris",
						"dict_annotations":dict_annotations_tm},
			"wgcna.mousebrain-190111":
						{"dataset":"mousebrain",
						"dict_annotations":dict_annotations_wgcna},
			}

print(dict_run)


#########################################################################################
###################################### DRY-RUN (TEST) ######################################
#########################################################################################

### Make sure that all_genes annotation is present, so program does not fail later on.
for run_name in dict_run.keys():
	param_dict = dict_run[run_name]
	ldsc_all_genes_ref_ld_chr_name = get_all_genes_ref_ld_chr_name(param_dict["dataset"])


#########################################################################################
########################### QUANTILE PREPROCESS ANNOTATIONS #############################
#########################################################################################


if FLAG_QUANTILE_FIXED:
	out_suffix = "qfixed"
else: 
	out_suffix = "q5_exclude_zero"


### Create job commands
list_cmds_quantile = []
# for FLAG_QUANTILE_FIXED in [True, False]:
# 	if FLAG_QUANTILE_FIXED:
# 		out_suffix = "qfixed"
# 	else: 
# 		out_suffix = "q5_exclude_zero"
for run_name, param_dict in dict_run.items():
	ldsc_all_genes_ref_ld_chr_name = get_all_genes_ref_ld_chr_name(param_dict["dataset"])
	for annotation_id in param_dict["dict_annotations"]:
		### [slightly overkill to check all chromosomes - one should be enough] check that all .annot files exists
		for chromosome in LIST_CHROMOSOMES:
			file_annot = "{}.{}.annot.gz".format(param_dict["dict_annotations"][annotation_id]["file_path_prefix"], chromosome)
			if not os.path.exists(file_annot):
				raise ValueError("file_annot does not exist: {}".format(file_annot))
		### set output filenames and check for exisistance
		file_out_q = "{}.{}.M".format(param_dict["dict_annotations"][annotation_id]["file_path_prefix"], out_suffix)
		file_out_log = "{}.{}.log".format(param_dict["dict_annotations"][annotation_id]["file_path_prefix"], out_suffix)
		if os.path.exists(file_out_q):
			print("file exists: {}. Will not make a new".format(file_out_q))
			continue
		### Command
		### NOTE: I'm note sure it is necessary to add the 'all genes', but I think it is since quantile_h2g.r relies on the .results and .q.M. files to have the same annotations (in the same order)
		cmd = """perl {PATH_LDSC_QUANTILE_PERL_SCRIPT} \
        --ref-annot-chr {PATH_LDSC_DATA_MAIN}/baseline_v1.1_thin_annot/baseline.,{ldsc_all_genes_ref_ld_chr_name},{file_annot_prefix}. \
        --frqfile-chr {PATH_LDSC_DATA_MAIN}/1000G_Phase3_frq/1000G.EUR.QC. \
        --annot-header {name_annot_context} \
        --nb-quantile 5 \
        --maf 0.05 \
        --thin-annot \
        {flag_dependent_argument} \
        --out {file_out} |& tee {file_out_log}""".format( # See here about tee and buffering : https://stackoverflow.com/questions/11337041/force-line-buffering-of-stdout-when-piping-to-tee
		PATH_LDSC_QUANTILE_PERL_SCRIPT=PATH_LDSC_QUANTILE_PERL_SCRIPT,
		PATH_LDSC_DATA_MAIN=PATH_LDSC_DATA_MAIN,
		ldsc_all_genes_ref_ld_chr_name=ldsc_all_genes_ref_ld_chr_name,
		file_annot_prefix=param_dict["dict_annotations"][annotation_id]["file_path_prefix"],  # this is a file prefix. quantile_M_fixed_non_zero_quantiles will add .<CHR>.annot.gz to the file
		name_annot_context=param_dict["dict_annotations"][annotation_id]["name_context"],
		flag_dependent_argument="--fixed-quantiles" if FLAG_QUANTILE_FIXED else "--exclude0",
		file_out=file_out_q,
		file_out_log=file_out_log)
		list_cmds_quantile.append(cmd)


### Call scheduler
job_scheduler(list_cmds=list_cmds_quantile, n_parallel_jobs=N_PARALLEL_QUANTILE_PERL_JOBS)

################## OUTPUT quantile_M_fixed_non_zero_quantiles.pl ##################

# ### The file contains number of SNPs (or annotation size='sum of annotation' for continuos annotation) for each annotation loading in the --ref-annot-chr argument.
# First column = N_SNPs/annotation_size in lowest quantile (of the annotation stratified on)
# Last column = N_SNPs/annotation_size in highest quantile (of the annotation stratified on)


### snippet qfixed | /scratch/sc-ldsc/celltypes.mousebrain.all/per_annotation/celltypes.mousebrain.all__mousebrain_all.TEINH12.sem_mean.qfixed.M
# N       4081518 412332  630022  365058  219945  252284
# base    4081518 412332  630022  365058  219945  252284
# Coding_UCSC.bed 36656   9055    15431   11046   6187    6720
# Coding_UCSC.extend.500.bed      165519  42589   69895   45784   27669   28042
# Conserved_LindbladToh.bed       100200  10548   17478   10365   6648    7874
# ...
# WeakEnhancer_Hoffman.bed        69748   11857   18109   11375   6738    7281
# WeakEnhancer_Hoffman.extend.500.bed     297928  49986   76299   47203   27705   30150
# all_genes_in_dataset.mousebrain 2085845 412332  630022  365058  219945  252284
# mousebrain_all.TEINH12.sem_mean 0       40918.9212836192        192911.49975945 178369.311397958        153360.932629917        230021.339404443

### snippet q5_exclude_zero | /scratch/sc-ldsc/celltypes.mousebrain.all/per_annotation/celltypes.mousebrain.all__mousebrain_all.TEINH12.sem_mean.q5_exclude_zero.M
# N       376102  376046  376194  377725  373574
# base    376102  376046  376194  377725  373574
# Coding_UCSC.bed 8153    8873    10093   10983   10337
# Coding_UCSC.extend.500.bed      38308   40848   44197   46535   44091
# Conserved_LindbladToh.bed       9821    10410   10155   10900   11627
# ...
# WeakEnhancer_Hoffman.bed        10642   10621   11515   11574   11008
# WeakEnhancer_Hoffman.extend.500.bed     45196   45151   47584   48031   45381
# all_genes_in_dataset.mousebrain 376102  376046  376194  377725  373574
# mousebrain_all.TEINH12.sem_mean 34073.4491382769        95629.6952472526        139792.10504084 206291.635054323        319795.119994766

#########################################################################################
###################################### RUN LDSC h2 ######################################
#########################################################################################


### Create job commands
list_cmds_ldsc_prim = []
for run_name, param_dict in dict_run.items():
	ldsc_all_genes_ref_ld_chr_name = get_all_genes_ref_ld_chr_name(param_dict["dataset"])
	for annotation_id in param_dict["dict_annotations"]:
		for gwas in list_gwas:
			# Output name: <runname/key_prim>__<annotation_id>__<gwas>
			fileout_prefix_ldsc_h2 = "/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc_h2/{run_name}__{annotation_id}__{gwas}".format(run_name=run_name, annotation_id=annotation_id, gwas=gwas)
			if os.path.exists("{}.results".format(fileout_prefix_ldsc_h2)):
				print("GWAS={}, run_name={},  annotation_id={} | LDSC outout file exists: {}. Will skip this LDSC regression...".format(gwas, run_name,  annotation_id, fileout_prefix_ldsc_h2))
				continue
			### I'm 90% sure that ldsc.py ONLY runs on python2 - and not python3.
			### *OBS*: we are runnin ldsc python script with UNBUFFERED stdout and stderr
			### REF: https://stackoverflow.com/questions/230751/how-to-flush-output-of-print-function
			### python -u: Force the stdout and stderr streams to be unbuffered. THIS OPTION HAS NO EFFECT ON THE STDIN STREAM [or writing of other files, e.g. the ldsc .log file]. See also PYTHONUNBUFFERED.
			cmd = """{PYTHON2_EXEC} {flag_unbuffered} {script} --h2 /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/{gwas}.sumstats.gz \
            --ref-ld-chr {PATH_LDSC_DATA_MAIN}/baseline_v1.1_thin_annot/baseline.,{ldsc_all_genes_ref_ld_chr_name},{file_annot_prefix}. \
            --frqfile-chr {PATH_LDSC_DATA_MAIN}/1000G_Phase3_frq/1000G.EUR.QC. \
            --w-ld-chr {PATH_LDSC_DATA_MAIN}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
            --overlap-annot \
            --thin-annot \
            --print-cov --print-coefficients --print-delete-vals \
            --out {fileout_prefix_ldsc_h2}""".format(
				PYTHON2_EXEC=PYTHON2_EXEC,
				flag_unbuffered="-u" if FLAG_UNBUFFERED else "",
				script=PATH_LDSC_SCRIPT,
				gwas=gwas,
				run_name=run_name,
				ldsc_all_genes_ref_ld_chr_name=ldsc_all_genes_ref_ld_chr_name,
				PATH_LDSC_DATA_MAIN=PATH_LDSC_DATA_MAIN,
				file_annot_prefix=param_dict["dict_annotations"][annotation_id]["file_path_prefix"],  # this is a file prefix. LDSC will add .<CHR>.l2.ldsc.gz to the file
				fileout_prefix_ldsc_h2=fileout_prefix_ldsc_h2
				)
            ### NOT sure why "--print-cov" is needed. It outputs ".cov" file.
			list_cmds_ldsc_prim.append(cmd)


### Call scheduler
job_scheduler(list_cmds=list_cmds_ldsc_prim, n_parallel_jobs=N_PARALLEL_LDSC_REGRESSION_JOBS)

###################################### OUTPUT --h2 ######################################

### Example
# celltypes.mousebrain__TEGLU4__BMI_UKBB_Loh2018.cov
# celltypes.mousebrain__TEGLU4__BMI_UKBB_Loh2018.delete
# celltypes.mousebrain__TEGLU4__BMI_UKBB_Loh2018.log
# celltypes.mousebrain__TEGLU4__BMI_UKBB_Loh2018.part_delete
# celltypes.mousebrain__TEGLU4__BMI_UKBB_Loh2018.results

#########################################################################################
################################ RUN RSCRIPT h2 QUANTILE ################################
#########################################################################################

# PATH_LDSC_H2_RSCRIPT_SCRIPT

### WIKI quantile_h2g.r
# Arg1: the file containing the sum of each annotation by quantile of the continuous annotation (e.g. .q5.M file) ['annotfile' inside the script]
# Arg2: Specify the prefix of the ldsc result outputs (.results and .part_delete files) ['resultfile' inside the script]
# Arg3: Specify the output filename ['outfile' inside the script]


### Create job commands
list_cmds_h2_rscript = []
for run_name, param_dict in dict_run.items():
	for annotation_id in param_dict["dict_annotations"]:
		for gwas in list_gwas:
			### Cache function
			fileout_rscript =  "/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc_h2/{run_name}__{annotation_id}__{gwas}.results_quantile_{quantile_mode}_h2".format(run_name=run_name, 
																																								  annotation_id=annotation_id, 
																																								  gwas=gwas,
																																								  quantile_mode=out_suffix)
			if os.path.exists(fileout_rscript):
				print("GWAS={}, run_name={},  annotation_id={} | Rscript output exists: {}. Will skip this LDSC regression...".format(gwas, run_name,  annotation_id, fileout_rscript))
				continue
			### Make sure that LDSC h2 results exists
			fileout_prefix_ldsc_h2 = "/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc_h2/{run_name}__{annotation_id}__{gwas}".format(run_name=run_name, annotation_id=annotation_id, gwas=gwas)
			if not os.path.exists("{}.results".format(fileout_prefix_ldsc_h2)):
				raise ValueError("GWAS={}, run_name={},  annotation_id={} | LDSC h2 .results.txt outoup file does NOT exist: {}. Cannot run Rscript".format(gwas, run_name,  annotation_id, fileout_prefix_ldsc_h2))
			### File q
			file_out_q = "{}.{}.M".format(param_dict["dict_annotations"][annotation_id]["file_path_prefix"], out_suffix)
			### CMD
			cmd = """/tools/R/3.4.3/bin/Rscript {PATH_LDSC_H2_RSCRIPT_SCRIPT} \
			{file_out_q} \
			{fileout_prefix_ldsc_h2} \
			{fileout_rscript}""".format(PATH_LDSC_H2_RSCRIPT_SCRIPT=PATH_LDSC_H2_RSCRIPT_SCRIPT,
										file_out_q=file_out_q,
										fileout_prefix_ldsc_h2=fileout_prefix_ldsc_h2,
										fileout_rscript=fileout_rscript)
			list_cmds_h2_rscript.append(cmd)


### Call scheduler
job_scheduler(list_cmds=list_cmds_h2_rscript, n_parallel_jobs=N_PARALLEL_H2_RSCRIPT_JOBS)


###################################### OUTPUT quantile_h2g.R ######################################
# This output file has one row for each quantile (starting with lowest values) and column summarizing the heritability explained by each quantile, its enrichment and corresponding standard error and P value.
# h2g     h2g_se  prop_h2g        prop_h2g_se     enr     enr_se  enr_pval
# 0.0468467691665829      0.00348227231769801     0.283751125243022       0.017433979820968       1.41876013045356        0.0871701758501825      3.7199370574134e-06
# 0.0404121018424434      0.00179110903328923     0.244776311690852       0.00620168879845983     1.22388201557721        0.0310084555740338      1.27746090080597e-11
# 0.0354534348018723      0.00113298619206756     0.214741639556595       0.00132708163840661     1.07370559107702        0.00663539208285362     2.7158960368993e-27
# 0.0292348058573656      0.0013806156620049      0.177075371596941       0.00708277648692335     0.885378842131785       0.0354139617978289      0.00148371877436259
# 0.0131509795779308      0.00262523459741717     0.0796555519125901      0.0162441807658585      0.398276420747926       0.0812206308043398      4.77095371863701e-12


print("Script is done!")



