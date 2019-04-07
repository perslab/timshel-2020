#!/usr/bin/env python2.7

import ConfigParser
import argparse
import subprocess
import time

###################################### Notes ######################################
### USAGE (should run in python2.7 environment because runs on python2.7)
# module unload anaconda
# module load anaconda/2-4.4.0
# python run_depict_from_template_cfg.py --cfg_template_file /projects/timshel/DEPICT/reproductive_behavior/TEMPLATE_reproductive_behavior.cfg \
# --gwas_summary_statistics_file /data/genetic-gwas/reproductive_behaviour_2018/chr1to22/AFB/meta_AFB_POOLED_chr1to22_pos_snpid.txt.gz \
# --label_for_output_files AFB_pooled_chr1to22_1e-5_maca_celltype \
# --tissue_expression_file /projects/tp/tmp-bmi-brain/data/maca/maca.per_celltype.celltype_expr.avg_expr.hsapiens_orthologs_std.tab

### DESCRIPTION
# 1) settings from cfg_template_file will be updated based on arguments passed to this script
# 2) an updated config file will be written to current working directory
# 3) DEPICT will be run with the updated config file.

### Assumptions
# - analysis path not updated from config. So use label_for_output_files instead
# - will assume that GWAS files are formatted identically. That is, template config file should fit the GWAS input file.


###################################### CONSTANTS ######################################
script_depict = "/scratch/tmp-depict-1.0.174/src/python/depict.py"

###################################### SCRIPT ######################################


# Read path to config file
parser = argparse.ArgumentParser()
parser.add_argument('--cfg_template_file', type=str, help='DEPICT TEMPLATE configuration file')
parser.add_argument('--gwas_summary_statistics_file', type=str, help='DEPICT TEMPLATE configuration file')
parser.add_argument('--label_for_output_files', type=str, help='DEPICT TEMPLATE configuration file')
parser.add_argument('--tissue_expression_file', type=str, help='DEPICT TEMPLATE configuration file')
# parser.add_argument('--association_pvalue_cutoff', type=str, help='DEPICT TEMPLATE configuration file')
args = parser.parse_args()


# Parse the config file
cfg = ConfigParser.ConfigParser()
cfg.read(args.cfg_template_file)

### Set/update template config file with inputs from this script
cfg.set("GWAS FILE SETTINGS",'gwas_summary_statistics_file', args.gwas_summary_statistics_file)
cfg.set("GWAS FILE SETTINGS",'label_for_output_files', args.label_for_output_files)
cfg.set("MISC SETTINGS",'tissue_expression_file', args.tissue_expression_file)
# cfg.set("GWAS FILE SETTINGS",'association_pvalue_cutoff', args.association_pvalue_cutoff)

### Write config file
cfg_file_out = "config-{}.cfg".format(args.label_for_output_files)
with open(cfg_file_out, 'w') as fh_out:
	cfg.write(fh_out)

### Run DEPICT
cmd = "python {} {}".format(script_depict, cfg_file_out)
print "Making call: {}".format(cmd)
p = subprocess.Popen(cmd, shell=True)
p.wait()
print "SCRIPT DONE. Returncode={}".format(p.returncode)



