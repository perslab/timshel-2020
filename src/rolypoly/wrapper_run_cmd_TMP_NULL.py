
import os
import sys
import time
import glob

import subprocess



###################################### USAGE ######################################

# time python wrapper_run_cmd.py |& tee wrapper_run_cmd.out.txt


################## INFERENCE ##################

N_CORES=8
LOG_PREFIX="null_gwas"
OUTDIR = "/scratch/sc-genetics"
if not os.path.exists(OUTDIR): # need to make it for the log files to go the same place.
    os.makedirs(OUTDIR)
list_cmds = [
"Rscript run_rolypoly_multi_expr_datasets-v3.R --gwas_name {gwas_name} --run_name {run_name} --outdir {outdir}/out.rolypoly_objs-v3.NULL_GWAS --gwas_linked_file /scratch/sc-genetics/out.rolypoly_objs-v3.NULL_GWAS/rolypoly_objs.{gwas_name}.tss.10kb.all_genes.gwas_linked.RData --expr_data_list /raid5/projects/timshel/sc-genetics/sc-genetics/data/list_expr_data.atlas.txt --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_bootstrap 100 --n_cores {n_cores} &> {outdir}/log.{log_prefix}.{gwas_name}-{run_name}.out.txt".format(gwas_name="NULL_GWAS_1KG_phase3_EUR_N1", run_name="tss.10kb.all_genes", outdir=OUTDIR, n_cores=N_CORES, log_prefix=LOG_PREFIX),
"Rscript run_rolypoly_multi_expr_datasets-v3.R --gwas_name {gwas_name} --run_name {run_name} --outdir {outdir}/out.rolypoly_objs-v3.NULL_GWAS --gwas_linked_file /scratch/sc-genetics/out.rolypoly_objs-v3.NULL_GWAS/rolypoly_objs.{gwas_name}.tss.10kb.all_genes.gwas_linked.RData --expr_data_list /raid5/projects/timshel/sc-genetics/sc-genetics/data/list_expr_data.atlas.txt --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_bootstrap 100 --n_cores {n_cores} &> {outdir}/log.{log_prefix}.{gwas_name}-{run_name}.out.txt".format(gwas_name="NULL_GWAS_1KG_phase3_EUR_N2", run_name="tss.10kb.all_genes", outdir=OUTDIR, n_cores=N_CORES, log_prefix=LOG_PREFIX),
]



################## GWAS LINKING ##################
# N_CORES=6
# LOG_PREFIX="null_gwas"
# OUTDIR = "/scratch/sc-genetics"
# if not os.path.exists(OUTDIR): # need to make it for the log files to go the same place.
# 	os.makedirs(OUTDIR)
# list_cmds = [
# # "Rscript run_rolypoly_multi_expr_datasets-v3.R --gwas_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_gwas_null/{gwas_name}.gwassumstats.rolypoly_fmt.tab.gz --gwas_name {gwas_name} --run_name {run_name} --outdir {outdir}/out.rolypoly_objs-v3.NULL_GWAS --expr_data_list /raid5/projects/timshel/sc-genetics/sc-genetics/data/list_expr_data.atlas.txt --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_bootstrap 100 --gwas_linking_only --n_cores {n_cores} &> {outdir}/log.{log_prefix}.{gwas_name}-{run_name}.out.txt".format(gwas_name="NULL_GWAS_1KG_phase3_EUR_N1", run_name="tss.10kb.all_genes", outdir=OUTDIR, n_cores=N_CORES, log_prefix=LOG_PREFIX),
# # "Rscript run_rolypoly_multi_expr_datasets-v3.R --gwas_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_gwas_null/{gwas_name}.gwassumstats.rolypoly_fmt.tab.gz --gwas_name {gwas_name} --run_name {run_name} --outdir {outdir}/out.rolypoly_objs-v3.NULL_GWAS --expr_data_list /raid5/projects/timshel/sc-genetics/sc-genetics/data/list_expr_data.atlas.txt --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_bootstrap 100 --gwas_linking_only --n_cores {n_cores} &> {outdir}/log.{log_prefix}.{gwas_name}-{run_name}.out.txt".format(gwas_name="NULL_GWAS_1KG_phase3_EUR_N2", run_name="tss.10kb.all_genes", outdir=OUTDIR, n_cores=N_CORES, log_prefix=LOG_PREFIX),
# "Rscript run_rolypoly_multi_expr_datasets-v3.R --gwas_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_gwas_null/{gwas_name}.gwassumstats.rolypoly_fmt.tab.gz --gwas_name {gwas_name} --run_name {run_name} --outdir {outdir}/out.rolypoly_objs-v3.NULL_GWAS --expr_data_list /raid5/projects/timshel/sc-genetics/sc-genetics/data/list_expr_data.atlas.txt --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_bootstrap 100 --gwas_linking_only --n_cores {n_cores} &> {outdir}/log.{log_prefix}.{gwas_name}-{run_name}.out.txt".format(gwas_name="NULL_GWAS_1KG_phase3_EUR_N3", run_name="tss.10kb.all_genes", outdir=OUTDIR, n_cores=N_CORES, log_prefix=LOG_PREFIX),
# "Rscript run_rolypoly_multi_expr_datasets-v3.R --gwas_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_gwas_null/{gwas_name}.gwassumstats.rolypoly_fmt.tab.gz --gwas_name {gwas_name} --run_name {run_name} --outdir {outdir}/out.rolypoly_objs-v3.NULL_GWAS --expr_data_list /raid5/projects/timshel/sc-genetics/sc-genetics/data/list_expr_data.atlas.txt --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_bootstrap 100 --gwas_linking_only --n_cores {n_cores} &> {outdir}/log.{log_prefix}.{gwas_name}-{run_name}.out.txt".format(gwas_name="NULL_GWAS_1KG_phase3_EUR_N4", run_name="tss.10kb.all_genes", outdir=OUTDIR, n_cores=N_CORES, log_prefix=LOG_PREFIX),
# "Rscript run_rolypoly_multi_expr_datasets-v3.R --gwas_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_gwas_null/{gwas_name}.gwassumstats.rolypoly_fmt.tab.gz --gwas_name {gwas_name} --run_name {run_name} --outdir {outdir}/out.rolypoly_objs-v3.NULL_GWAS --expr_data_list /raid5/projects/timshel/sc-genetics/sc-genetics/data/list_expr_data.atlas.txt --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_bootstrap 100 --gwas_linking_only --n_cores {n_cores} &> {outdir}/log.{log_prefix}.{gwas_name}-{run_name}.out.txt".format(gwas_name="NULL_GWAS_1KG_phase3_EUR_N5", run_name="tss.10kb.all_genes", outdir=OUTDIR, n_cores=N_CORES, log_prefix=LOG_PREFIX),
# "Rscript run_rolypoly_multi_expr_datasets-v3.R --gwas_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_gwas_null/{gwas_name}.gwassumstats.rolypoly_fmt.tab.gz --gwas_name {gwas_name} --run_name {run_name} --outdir {outdir}/out.rolypoly_objs-v3.NULL_GWAS --expr_data_list /raid5/projects/timshel/sc-genetics/sc-genetics/data/list_expr_data.atlas.txt --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_bootstrap 100 --gwas_linking_only --n_cores {n_cores} &> {outdir}/log.{log_prefix}.{gwas_name}-{run_name}.out.txt".format(gwas_name="NULL_GWAS_1KG_phase3_EUR_N6", run_name="tss.10kb.all_genes", outdir=OUTDIR, n_cores=N_CORES, log_prefix=LOG_PREFIX),
# "Rscript run_rolypoly_multi_expr_datasets-v3.R --gwas_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_gwas_null/{gwas_name}.gwassumstats.rolypoly_fmt.tab.gz --gwas_name {gwas_name} --run_name {run_name} --outdir {outdir}/out.rolypoly_objs-v3.NULL_GWAS --expr_data_list /raid5/projects/timshel/sc-genetics/sc-genetics/data/list_expr_data.atlas.txt --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_bootstrap 100 --gwas_linking_only --n_cores {n_cores} &> {outdir}/log.{log_prefix}.{gwas_name}-{run_name}.out.txt".format(gwas_name="NULL_GWAS_1KG_phase3_EUR_N7", run_name="tss.10kb.all_genes", outdir=OUTDIR, n_cores=N_CORES, log_prefix=LOG_PREFIX),
# "Rscript run_rolypoly_multi_expr_datasets-v3.R --gwas_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_gwas_null/{gwas_name}.gwassumstats.rolypoly_fmt.tab.gz --gwas_name {gwas_name} --run_name {run_name} --outdir {outdir}/out.rolypoly_objs-v3.NULL_GWAS --expr_data_list /raid5/projects/timshel/sc-genetics/sc-genetics/data/list_expr_data.atlas.txt --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_bootstrap 100 --gwas_linking_only --n_cores {n_cores} &> {outdir}/log.{log_prefix}.{gwas_name}-{run_name}.out.txt".format(gwas_name="NULL_GWAS_1KG_phase3_EUR_N8", run_name="tss.10kb.all_genes", outdir=OUTDIR, n_cores=N_CORES, log_prefix=LOG_PREFIX),
# "Rscript run_rolypoly_multi_expr_datasets-v3.R --gwas_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_gwas_null/{gwas_name}.gwassumstats.rolypoly_fmt.tab.gz --gwas_name {gwas_name} --run_name {run_name} --outdir {outdir}/out.rolypoly_objs-v3.NULL_GWAS --expr_data_list /raid5/projects/timshel/sc-genetics/sc-genetics/data/list_expr_data.atlas.txt --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_bootstrap 100 --gwas_linking_only --n_cores {n_cores} &> {outdir}/log.{log_prefix}.{gwas_name}-{run_name}.out.txt".format(gwas_name="NULL_GWAS_1KG_phase3_EUR_N9", run_name="tss.10kb.all_genes", outdir=OUTDIR, n_cores=N_CORES, log_prefix=LOG_PREFIX),
# "Rscript run_rolypoly_multi_expr_datasets-v3.R --gwas_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_gwas_null/{gwas_name}.gwassumstats.rolypoly_fmt.tab.gz --gwas_name {gwas_name} --run_name {run_name} --outdir {outdir}/out.rolypoly_objs-v3.NULL_GWAS --expr_data_list /raid5/projects/timshel/sc-genetics/sc-genetics/data/list_expr_data.atlas.txt --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_bootstrap 100 --gwas_linking_only --n_cores {n_cores} &> {outdir}/log.{log_prefix}.{gwas_name}-{run_name}.out.txt".format(gwas_name="NULL_GWAS_1KG_phase3_EUR_N10", run_name="tss.10kb.all_genes", outdir=OUTDIR, n_cores=N_CORES, log_prefix=LOG_PREFIX),
# ]



# list_cmds = ["sleep 10"] * 10 # ['sleep 3', 'sleep 3', 'sleep 3']


list_of_processes = []
batch = 1
for i, cmd in enumerate(list_cmds, start=1):
    print "batch = {} | i = {} | Running command: {}".format(batch, i, cmd)
    print cmd
    p = subprocess.Popen(cmd, shell=True)
    list_of_processes.append(p)
    print "batch = {} | i = {} | list_of_processes:".format(batch, i)
    print list_of_processes
    if i % 1 == 0:
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
            


# print "PYTHON SCRIPT DONE"