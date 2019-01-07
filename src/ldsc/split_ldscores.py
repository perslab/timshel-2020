
from __future__ import print_function
import argparse
import glob
import os
import subprocess


### PT adds
import re
import pandas as pd
import multiprocessing
import os


import pdb

###################################### USAGE ######################################
# Compatibility: Python 2 and 3

# time python split_ldscores.py --prefix_ldscore_files /scratch/sc-ldsc/nn_lira_sema-combined/nn_lira_sema
# time python split_ldscores.py --prefix_ldscore_files /scratch/sc-ldsc/maca/
# time python split_ldscores.py --prefix_ldscore_files /scratch/sc-ldsc/mousebrain/
# time python split_ldscores.py --prefix_ldscore_files /scratch/sc-ldsc/hypothalamus_mette_thesis/

###################################### TODO ######################################

# Fix reliance on 'COMBINED_ANNOT' in file name.

###################################### DESCRIPTION ######################################

# Purpose: split ldscore files (l2.ldscore.gz) per annotation.


# N output files = 3000*22*3= ~200k!!!
# put in DIR per chromosome (22)
# pit in DIR per annotation (~3000) containing 22*3 files

### Steps
# get list of ldscore files
# PER CHR
	# read ldscore file. Warn if file does not exists.
	# write out
# NB: make_cts will check that all chromosomes are present for a given annotation.
###################################### USAGE ######################################


def split_M_files(file_ldscore, chromosome, prefix_genomic_annot, list_annotations):
	### Both file_M and file_M_5_50 contains a single line with the same number of fields as the number of annotations in the ldscore file
	# nn_lira_sema.COMBINED_ANNOT.6.l2.ldscore.gz
	# nn_lira_sema.COMBINED_ANNOT.6.l2.M
	# nn_lira_sema.COMBINED_ANNOT.6.l2.M_5_50
	print("CHR={} | Writing .M and .M_5_50 files".format(chromosome))
	file_ldscore_base = re.sub(r"\.l2\.ldscore\.gz$", "", file_ldscore) # .e.g /scratch/sc-ldsc/nn_lira_sema/nn_lira_sema.COMBINED_ANNOT.9.
	file_M = "{}.l2.M".format(file_ldscore_base)
	file_M_5_50 = "{}.l2.M_5_50".format(file_ldscore_base)
	with open(file_M, "r") as fh_M, open(file_M_5_50, "r") as fh_M_5_50:
		list_M =  fh_M.readline().rstrip().split()
		list_M_5_50 = fh_M_5_50.readline().rstrip().split()
	assert(len(list_M) == len(list_M_5_50) == len(list_annotations))
	for i in range(len(list_annotations)):
		annotation = list_annotations[i]
		M = list_M[i]
		M_5_50 = list_M_5_50[i]
		file_out_M = "{}/{}__{}.{}.l2.M".format(out_dir, prefix_genomic_annot, annotation, chromosome)
		file_out_M_5_50 = "{}/{}__{}.{}.l2.M_5_50".format(out_dir, prefix_genomic_annot, annotation, chromosome)
		with open(file_out_M, "w") as fh_out_M, open(file_out_M_5_50, "w") as fh_out_M_5_50:
			fh_out_M.write(M)
			fh_out_M_5_50.write(M_5_50)
	print("CHR={} | DONE writing .M and .M_5_50 files".format(chromosome))

def all_files_exist(chromosome, prefix_genomic_annot, annotations_clean):
	print("CHR={} | Checking if all files already exists...".format(chromosome))
	# generate filenames for all files to be generated (for this chromosome):
	set_all_files_to_be_generated = set()
	for annotation in annotations_clean:
		# NB: we only generate the basefilename - not the full filepath - because of how os.listdir() works. See below.
		file_out_M = "{}__{}.{}.l2.M".format(prefix_genomic_annot, annotation, chromosome)
		file_out_M_5_50 = "{}__{}.{}.l2.M_5_50".format(prefix_genomic_annot, annotation, chromosome)
		file_out_ldscore = "{}__{}.{}.l2.ldscore.gz".format(prefix_genomic_annot, annotation, chromosome)
		set_all_files_to_be_generated.update([file_out_M, file_out_M_5_50, file_out_ldscore]) # .update(): add multiple elements to set
	# get all/any existing filenames (across all chromosome, but that does not matter - we only loose a bit of speed):
	bool_all_files_exists = set_all_files_to_be_generated.issubset(os.listdir(out_dir)) # test if set_all_files_to_be_generated is a subset of all files in out_dir. Returns True if it is a subset (i.e. all files already exists)
	# ^ OBS: os.listdir() only gives the basename of the files in the dirs - not their full paths. E.g. nn_lira_sema.springgreen.1.l2.ldscore.gz.
	return bool_all_files_exists

def split_ldscore_file_per_annotation(file_ldscore):
	print("Processing file_ldscore {}".format(file_ldscore))
	m = re.search(r"(.*)\.COMBINED_ANNOT\.(\d{1,2})\.l2.ldscore.gz$", os.path.basename(file_ldscore)) # file_ldscore e.g. "/scratch/sc-ldsc/nn_lira_sema/nn_lira_sema.COMBINED_ANNOT.9.l2.ldscore.gz"
	prefix_genomic_annot = m.groups()[0] # e.g nn_lira_sema
	chromosome = m.groups()[1]
	df = pd.read_csv(file_ldscore, sep="\t") # no index
	# CHR     SNP     BP      blueL2  cornsilkL2 ...
	annotations_header = df.columns[3:].tolist() # skip the first three columns (CHR, SNP, BP)
	annotations_clean = [re.sub(r"L2$", "", x) for x in annotations_header] # remove trailing L2 in name. Do not use x.rstrip("L2")
	if all_files_exist(chromosome, prefix_genomic_annot, annotations_clean):
		print("CHR={} | All files already exists! Will do nothing for this chromosome".format(chromosome))
		return None # do nothing
	else:
		split_M_files(file_ldscore, chromosome, prefix_genomic_annot, annotations_clean)
		for counter, annotation in enumerate(annotations_header):
			if counter % 25 == 0:
				print("CHR={} #{}/#{}| Writing annot files".format(chromosome, counter, len(annotations_header)))
			annotation_clean = re.sub(r"L2$", "", annotation)
			file_out_ldscore = "{}/{}__{}.{}.l2.ldscore.gz".format(out_dir, prefix_genomic_annot, annotation_clean, chromosome)
			df[["CHR", "SNP", "BP", annotation]].to_csv(file_out_ldscore, sep="\t", index=False, compression="gzip") # no index




###################################### MAIN ######################################

parser = argparse.ArgumentParser()

parser = argparse.ArgumentParser()
parser.add_argument('--prefix_ldscore_files', type=str, help="""File path prefix to LDSC ldscore files to split per annotation. 
	Will split ldscores files matching <prefix_ldscore_files>*.l2.ldscore.gz
	Notice that if prefix_ldscore_files is a directory (with a trailing '/'), then all ldscore file in the directory will be used.
	Example with prefix_genomic_annot as input: /scratch/sc-ldsc/ldsc_analysis_dir/files_prefix
	Example with dir as input: /scratch/sc-ldsc/ldsc_analysis_dir/ [OBS the trailing slash is important!]""")
parser.add_argument('--n_parallel_jobs', type=int, default=22, help='Number of processes. Default 22 (the number of chromosomes)')


args = parser.parse_args()


main_dir = os.path.dirname(args.prefix_ldscore_files)
out_dir = os.path.join(main_dir, "per_annotation")

### Make out_dir
if not os.path.exists(out_dir):
    print("Making output dir {}".format(out_dir))
    os.makedirs(out_dir)


file_ldscores = glob.glob("{}*.l2.ldscore.gz".format(args.prefix_ldscore_files)) # e.g. <SOMEPATH>/nn_lira_sema.COMBINED_ANNOT.22.l2.ldscore.gz


# ### Check for existing annot files
# pool = multiprocessing.Pool(processes=len(list_chromosomes))
# list_chromosomes_to_run = pool.map(check_annot_file, list_chromosomes) # check_annot_file returns None if annot file exists and is ok
# list_chromosomes_to_run = [x for x in list_chromosomes_to_run if x is not None] # filter away None

# print("========================== CHROMOSOMES TO RUN ====================")
# print("N chromosomes = {}".format(len(list_chromosomes_to_run)))
# print("Chromosomes: {}".format(",".join(map(str, list_chromosomes_to_run))))
# print("=========================================================================")


print("Starting pool...")
pool = multiprocessing.Pool(processes=min(args.n_parallel_jobs, len(file_ldscores)))
pool.map(split_ldscore_file_per_annotation, file_ldscores)
# split_ldscore_file_per_annotation(file_ldscores[0]) # for debugging

print("Script is done!")







