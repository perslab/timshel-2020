
from __future__ import print_function
import argparse
import glob
import os
import subprocess

###################################### USAGE ######################################
# Compatibility: Python 2 and 3

###################################### TODO ######################################

### ARGS
# GWAS
# CTS file
# Outdir
# Flag to run without BASELINE
# Flag to run with/without all_genes --> give all genes file prefix


# ### Run the regression
# GWAS=BMI_Yengo2018
# CTS=novo_bulk_sema_lira
# python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/ldsc.py \
#     --h2-cts /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}.sumstats.gz \
#     --ref-ld-chr /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_baseline/baseline. \
#     --out /raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/wgcna_modules/${CTS}_${GWAS} \
#     --ref-ld-chr-cts /XXX/${CTS}.ldcts \
#     --w-ld-chr /raid5/projects/timshel/sc-genetics/ldsc/data/weights_hm3_no_hla/weights.

###################################### DESCRIPTION ######################################


### Output
# This script will call ldsc.py to do CTS regression.
# The following output files will be written to the --prefix_annot_files:
# <OUT>.cell_type_results.txt
# <OUT>.log


###################################### CONSTANTS ######################################


python_exec = "/tools/anaconda/2-4.4.0/bin/python2" # runs on python2

###################################### ARGS ######################################



# CTS = "wgcna.maca"
# CTS = "wgcna.mousebrain"
# CTS = "wgcna.hypothalamus_mette_thesis"


###################################### FUNCTIONS ######################################
def preprocess(prefix_out, file_multi_gene_set):
	### Make annot
	cmd = """{python_exec} make_annot_from_geneset_all_chr.py \
	--file_multi_gene_set {file_multi_gene_set} \
	--file_gene_coord /raid5/projects/timshel/sc-genetics/ldsc/data/gene_coords/gene_annotation.hsapiens_all_genes.GRCh37.ens_v91.LDSC_fmt.txt \
	--windowsize 100000 \
	--bimfile_basename /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC \
	--out_dir /scratch/sc-ldsc/{prefix_out} \
	--out_prefix {prefix_out} \
	--flag_wgcna \
	--flag_mouse""".format(python_exec=python_exec, file_multi_gene_set=file_multi_gene_set, prefix_out=prefix_out)

	print("Running command: {}".format(cmd))
	p = subprocess.Popen(cmd, shell=True)
	p.wait()
	print("Return code: {}".format(p.returncode))


	### compute LD scores
	cmd="{python_exec} wrapper_compute_ldscores.py --prefix_annot_files /scratch/sc-ldsc/{prefix_out}/ --n_parallel_jobs 4".format(python_exec=python_exec, prefix_out=prefix_out)
	print("Running command: {}".format(cmd))
	p = subprocess.Popen(cmd, shell=True)
	p.wait()
	print("Return code: {}".format(p.returncode))
	# RUNTIME ----> ~6 h for ~500 modules with --n_parallel_jobs=4


	### split LD scores
	cmd="{python_exec} split_ldscores.py --prefix_ldscore_files /scratch/sc-ldsc/{prefix_out}/".format(python_exec=python_exec, prefix_out=prefix_out)
	print("Running command: {}".format(cmd))
	p = subprocess.Popen(cmd, shell=True)
	p.wait()
	print("Return code: {}".format(p.returncode))
	# RUNTIME ----> ~10 min


	### make cts file
	cmd="{python_exec} make_cts_file.py --prefix_ldscore_files /scratch/sc-ldsc/{prefix_out}/per_annotation/ --cts_outfile /raid5/projects/timshel/sc-genetics/sc-genetics/src/ldsc/cts_files/wgcna.{prefix_out}.ldcts.txt".format(python_exec=python_exec, prefix_out=prefix_out)
	# ^*OBS***:DIRTY USING wgcna. as prefix in  wgcna.{prefix_out}.ldcts.txt. FIX THIS.
	print("Running command: {}".format(cmd))
	p = subprocess.Popen(cmd, shell=True)
	p.wait()
	print("Return code: {}".format(p.returncode))
	# RUNTIME ----> 0 min

###################################### CALL FUNCTION ######################################


dict_CTS = {"mousebrain_Ependymal":"/raid5/projects/jonatan/tmp-mousebrain/tables/mousebrain_Ependymal_ClusterName_2_cell_cluster_module_genes.csv",
"mousebrain_Immune":"/raid5/projects/jonatan/tmp-mousebrain/tables/mousebrain_Immune_ClusterName_2_cell_cluster_module_genes.csv",
"mousebrain_PeripheralGlia":"/raid5/projects/jonatan/tmp-mousebrain/tables/mousebrain_PeripheralGlia_ClusterName_2_cell_cluster_module_genes.csv",
"mousebrain_Vascular":"/raid5/projects/jonatan/tmp-mousebrain/tables/mousebrain_Vascular_ClusterName_2_cell_cluster_module_genes.csv",
"mousebrain_Astrocytes":"/raid5/projects/jonatan/tmp-mousebrain/tables/mousebrain_Astrocytes_ClusterName_2_cell_cluster_module_genes.csv",
"mousebrain_Oligos":"/raid5/projects/jonatan/tmp-mousebrain/tables/mousebrain_Oligos_ClusterName_2_cell_cluster_module_genes.csv"}

for prefix_out, file_multi_gene_set in dict_CTS.items():
	preprocess(prefix_out, file_multi_gene_set)


###################################### CALL LDSC ######################################

path_ldsc_script = "/raid5/projects/timshel/sc-genetics/ldsc/ldsc/ldsc.py" 

n_parallel_jobs = 4

list_GWAS = ["BMI_Yengo2018",
"EA3_Lee2018",
"SCZ_Ripke2014",
"HEIGHT_Yengo2018"]

list_cmds = []
for CTS in dict_CTS:
	for GWAS in list_GWAS:
		cmd = """{python_exec} {script} --h2-cts /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/{GWAS}.sumstats.gz \
	    --ref-ld-chr /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_baseline/baseline. \
	    --out /raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/wgcna.{CTS}.{GWAS} \
	    --ref-ld-chr-cts /raid5/projects/timshel/sc-genetics/sc-genetics/src/ldsc/cts_files/wgcna.{CTS}.ldcts.txt \
	    --w-ld-chr /raid5/projects/timshel/sc-genetics/ldsc/data/weights_hm3_no_hla/weights.""".format(
			python_exec=python_exec,
			script=path_ldsc_script,
			GWAS=GWAS,
			CTS=CTS
			)
		list_cmds.append(cmd)


### You need to keep devnull open for the entire life of the Popen object, not just its construction. 
# FNULL = open(os.devnull, 'w') # devnull filehandle does not need to be closed?

list_of_processes = []
batch = 1
for i, cmd in enumerate(list_cmds, start=1):
	print("batch = {} | i = {} | Running command: {}".format(batch, i, cmd))
	## p = subprocess.Popen(cmd, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
	p = subprocess.Popen(cmd, shell=True)
	list_of_processes.append(p)
	print("batch = {} | i = {} | PIDs of running jobs (list_of_processes):".format(batch, i))
	print(" ".join([str(p.pid) for p in list_of_processes])) # print PIDs
	if i % n_parallel_jobs == 0: # JOB BATCH SIZE
		batch += 1
		for p in list_of_processes:
			print("=========== Waiting for process: {} ===========".format(p.pid))
			p.wait()
			print("Returncode = {}".format(p.returncode))
		list_of_processes = [] # 'reset' list

### wait for the rest for the rest of the processes
for p in list_of_processes:
	print("=========== Waiting for process: {} ===========".format(p.pid))
	p.wait()



print("Script is done!")



