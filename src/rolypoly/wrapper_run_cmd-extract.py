
import os
import sys
import time
import glob

import subprocess


###################################### USAGE ######################################

# time python wrapper_run_cmd.py |& tee wrapper_run_cmd.out.txt

###################################### XXX ######################################

script2run = "/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/rolypoly_extract_inference_results-BACKUP.R"


files_in = glob.glob("/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v1/*.final.RData")
# files_in = glob.glob("/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v1.clean/*.inference.RData")

list_cmds = ["Rscript {} --input_file {}".format(script2run, filepath) for filepath in files_in]

# time Rscript rolypoly_extract_inference_results.R --input_file /projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v1/out.rolypoly_objs.disease_HI_CHOL_SELF_REP.squared_tss_10kb.final.RData
# time Rscript rolypoly_extract_inference_results.R --input_file /projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v1/out.rolypoly_objs.body_BMI_Locke2015.squared_tss_10kb.final.RData

print list_cmds
print len(list_cmds)

# list_cmds = ["sleep 10"] * 10 # ['sleep 3', 'sleep 3', 'sleep 3']



list_of_processes = []
batch = 1
for i, cmd in enumerate(list_cmds, start=1):
    print "batch = {} | i = {} | Running command: {}".format(batch, i, cmd)
    p = subprocess.Popen(cmd, shell=True)
    list_of_processes.append(p)
    print "batch = {} | i = {} | list_of_processes:".format(batch, i)
    print list_of_processes
    if i % 20 == 0:
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
            


print "PYTHON SCRIPT DONE"