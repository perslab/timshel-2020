
import os
import sys
import time
import glob

import subprocess


###################################### USAGE ######################################

# time python wrapper_run_rolypoly_rscript.py |& tee wrapper_run_rolypoly_rscript.out.txt

###################################### XXX ######################################


SCRIPT2RUN = "/raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/run_rolypoly_multi_expr_datasets-v2.R"
BATCH_SIZE = 3 # maximal number of processes to run cuncurrently



DIR_GWAS = "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats"
OUTDIR = "/raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2"


#DIR_GWAS = "/raid5/projects/timshel/sc-genetics/sc-genetics/src/preparation_steps"
list_files = glob.glob("{}/*.rolypoly_fmt.tab.gz".format(DIR_GWAS))
for x in list_files: print x

# time Rscript run_rolypoly_multi_expr_datasets.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_WHR_Shungin2015.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_WHR_Shungin2015 --run_name squared_tss_10kb --n_cores 10 |& tee log.body_WHR_Shungin2015.out.txt
# Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file {} --gwas_name {} --run_name {} --outdir {} --window_position tss --window_size_kb 10 --pos_transformation square --n_cores 10 > log.{}.out.txt


###################################### XXXX ######################################

list_of_processes = []
batch_counter = 1
for i, gwas_file in enumerate(list_files, start=1):
    gwas_name = os.path.basename(gwas_file).split(".")[0] # e.g. mental_NEUROTICISM
        # /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_alkes_UKBB/mental_NEUROTICISM.gwassumstats.rolypoly_fmt.tab.gz --> mental_NEUROTICISM
    print "batch_counter = {} | i = {} / {} | Running gwas_name: {}".format(batch_counter, i, len(list_files), gwas_name)

    cmd = "Rscript {} --gwas_file {} --gwas_name {} --run_name tss.10kb.square.all_genes --outdir {} --window_position tss --window_size_kb 10 --pos_transformation square --n_cores 10 &> log.{}.out.txt".format(SCRIPT2RUN,
        gwas_file,
        gwas_name,
        OUTDIR,
        gwas_name)

    print "Running command: {}".format(cmd)
    
    FNULL = open(os.devnull, 'w')
    p = subprocess.Popen(cmd, stdout=FNULL, stderr=subprocess.STDOUT, shell=True)
    list_of_processes.append(p)
    print "batch_counter = {} | i = {} / {} | list_of_processes:".format(batch_counter, i, len(list_files))
    print list_of_processes
    if i % BATCH_SIZE == 0:
        batch_counter += 1
        for process in list_of_processes:
            print "=========== Waiting for process: {} ===========".format(process.pid)
            process.wait()
        list_of_processes = [] # 'reset' list

### wait for the rest for the rest of the processes
print "*********** DONE IN MAIN LOOP*************."
for process in list_of_processes:
    print "=========== Waiting for process: {} ===========".format(process.pid)
    process.wait()
            


print "PYTHON SCRIPT DONE"



###################################### OLD METHOD: single jobs - works ######################################

# for gwas_file in list_files:
#     gwas_name = os.path.basename(gwas_file).split(".")[0] # e.g. mental_NEUROTICISM
#         # /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_alkes_UKBB/mental_NEUROTICISM.gwassumstats.rolypoly_fmt.tab.gz --> mental_NEUROTICISM
    
#     cmd = "Rscript {} --gwas_file {} --gwas_name {} --run_name tss.10kb.square.all_genes --outdir {} --window_position tss --window_size_kb 10 --pos_transformation square --n_cores 10 > log.{}.out.txt".format(SCRIPT2RUN,
#         gwas_file,
#         gwas_name,
#         outdir,
#         gwas_name)
#     print "Running command: {}".format(cmd)
#     p = subprocess.Popen(cmd, shell=True)
#     print "Waiting..."
#     p.wait()
