
from __future__ import print_function
import argparse
import glob
import os
import subprocess
import re
import configparser

prefix_genomic_annot="control.all_genes_in_dataset"

#########################################################################################
###################################### PRE-PROCESS ######################################
#########################################################################################
if __name__ == "__main__":


	#########################################################################################
	###################################### PARSE VARIABLES ##################################
	#########################################################################################
	
	## Command line options
	parser = argparse.ArgumentParser()

	# parser.add_argument("--esmetric", default='sem_mean', help="Produce a binary annotation file to calculate LDSC with.")
	# parser.add_argument("--celltypes", default='all', help="Produce a binary annotation file to calculate LDSC with.")
	parser.add_argument("--binary", action='store_true', help="Produce a binary annotation file to calculate LDSC with.")
	parser.add_argument("--wgcna", action='store_true', help="Tells the script that the input multi-geneset file is from WGCNA.")
	parser.add_argument("--windowsize", type=int, default=100, help="Specify an alternate window size in kb, default is 100.")
	parser.add_argument("--buffered", action='store_true', help="Set the output to be buffered rather than the default unbuffered.")
	# parser.add_argument("--dataset", action='store_true', help="Set the output to be buffered rather than the default unbuffered.")

	args = parser.parse_args()
	
	FLAG_BINARY = args.binary
	FLAG_WGCNA = args.wgcna
	WINDOW_SIZE_KB = args.windowsize
	FLAG_UNBUFFERED = not args.buffered

	## Config file
	config = configparser.SafeConfigParser(os.environ)
	config.read('workflow_ldsc_config.ini')

	# Paths
	PYTHON2_EXEC = config['PATHS']['PYTHON2_EXEC']
	PYTHON3_EXEC = config['PATHS']['PYTHON3_EXEC']
	PATH_LDSC_SCRIPT = config['PATHS']['PATH_LDSC_SCRIPT']
	OUTPUT_DIR = config['PATHS']['OUTPUT_DIR']
	GENE_COORDINATES = config['PATHS']['GENE_COORDINATES']
	BIMFILES_BASENAME = config['PATHS']['BIMFILES_BASENAME']
	MULTIGENESET_DIRECTORY = config['PATHS']['MULTIGENESET_DIRECTORY']
	FILE_MULTIGENESET = MULTIGENESET_DIRECTORY + '/multi_geneset.all_genes_in_dataset.txt'


	### Make annot
	###  *RESOURCE NOTES*: if you have many modules (~3000-6000) then set --n_parallel_jobs to ~2-5 (instead of 22). Otherwise the script will up all the MEMORY on yggdrasil and fail.
	cmd = """{PYTHON3_EXEC} make_annot_from_geneset_all_chr.py \
	--file_multi_gene_set {file_multi_gene_set} \
	--file_gene_coord {file_gene_coordinates} \
	--windowsize {window_size} \
	--bimfile_basename {bimfiles_basename} \
	{flag_binary} \
	{flag_wgcna} \
	--out_dir {output_dir}/pre-computation/{prefix_genomic_annot} \
	--out_prefix {prefix_genomic_annot}
	""".format(PYTHON3_EXEC=PYTHON3_EXEC, 
		file_multi_gene_set=FILE_MULTIGENESET, 
		prefix_genomic_annot=prefix_genomic_annot, 
		flag_wgcna="--flag_wgcna --flag_mouse" if FLAG_WGCNA else "",
		flag_binary="--flag_encode_as_binary_annotation" if FLAG_BINARY else "",
		window_size = WINDOW_SIZE_KB * 1000,
		file_gene_coordinates = GENE_COORDINATES,
		output_dir = OUTPUT_DIR,
		bimfiles_basename = BIMFILES_BASENAME
		) 
	# --n_parallel_jobs 11

	print("Running command: {}".format(cmd))
	p = subprocess.Popen(cmd, shell=True)
	p.wait()
	print("Return code: {}".format(p.returncode))
	if not p.returncode == 0:
		raise Exception("Got non zero return code")


	### compute LD scores
	### *RESOURCE NOTES*: this script uses a lot of CPU. Never run more than 4 parallel jobs. 4 parallel jobs will use ~220% CPU
	cmd="{PYTHON3_EXEC} wrapper_compute_ldscores.py --prefix_annot_files {output_dir}/pre-computation/{prefix_genomic_annot}/ --n_parallel_jobs 2".format(PYTHON3_EXEC=PYTHON3_EXEC, prefix_genomic_annot=prefix_genomic_annot, output_dir = OUTPUT_DIR)
	print("Running command: {}".format(cmd))
	p = subprocess.Popen(cmd, shell=True)
	p.wait()
	print("Return code: {}".format(p.returncode))
	# RUNTIME ----> ~6 h for ~500 modules with --n_parallel_jobs=4
	if not p.returncode == 0:
		raise Exception("Got non zero return code")

	### split LD scores
	### This script will read 1 ".COMBINED_ANNOT.$CHR.l2.ldscore.gz" file  (N_SNPs x N_Modules) per parallel process.
	###  *RESOURCE NOTES*: this script may use quiet a lot of memory for many modules? Not sure
	cmd="{PYTHON3_EXEC} split_ldscores.py --prefix_ldscore_files {output_dir}/pre-computation/{prefix_genomic_annot}/ --n_parallel_jobs 22".format(PYTHON3_EXEC=PYTHON3_EXEC, prefix_genomic_annot=prefix_genomic_annot, output_dir = OUTPUT_DIR)
	print("Running command: {}".format(cmd))
	p = subprocess.Popen(cmd, shell=True)
	p.wait()
	print("Return code: {}".format(p.returncode))
	# RUNTIME ----> ~10 min
	if not p.returncode == 0:
		raise Exception("Got non zero return code")

	###################################### XXXX ######################################

	print("Script is done!")



