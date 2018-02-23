
import os
import sys
import time
import glob

import subprocess


###################################### USAGE ######################################

# time python wrapper_run_cmd.py |& tee wrapper_run_cmd.out.txt

###################################### XXX ######################################

# ### BMI
# list_cmds = ["Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_BMI_Locke2015 --run_name gene.10kb.squared.protein_coding_only --window_position gene --window_size_kb 10 --pos_transformation square --protein_coding_only --n_cores 10 &> log.body_BMI_Locke2015-gene.10kb.squared.protein_coding_only.out.txt",
# "Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_BMI_Locke2015 --run_name gene.10kb.squared.all_genes --window_position gene --window_size_kb 10 --pos_transformation square --n_cores 10 &> log.body_BMI_Locke2015-gene.10kb.squared.all_genes.out.txt",
# "Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_BMI_Locke2015 --run_name tss.10kb.squared.all_genes --window_position tss --window_size_kb 10 --pos_transformation square --n_cores 10 &> log.body_BMI_Locke2015-tss.10kb.squared.all_genes.out.txt",
# "Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_BMI_Locke2015 --run_name tss.20kb.squared.all_genes --window_position tss --window_size_kb 20 --pos_transformation square --n_cores 10 &> log.body_BMI_Locke2015-tss.20kb.squared.all_genes.out.txt",
# "Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_BMI_Locke2015 --run_name tss.10kb.pos_only.all_genes --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_cores 10 &> log.body_BMI_Locke2015-tss.10kb.pos_only.all_genes.out.txt",
# "Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_BMI_Locke2015 --run_name tss.10kb.abs.all_genes --window_position tss --window_size_kb 10 --pos_transformation abs --n_cores 10 &> log.body_BMI_Locke2015-tss.10kb.abs.all_genes.out.txt"
# ]

# ### SCZ
# list_cmds = ["Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/mental_SCZ_Ripke2014.gwassumstats.rolypoly_fmt.tab.gz --gwas_name mental_SCZ_Ripke2014 --run_name gene.10kb.squared.protein_coding_only --window_position gene --window_size_kb 10 --pos_transformation square --protein_coding_only --n_cores 15 &> log.mental_SCZ_Ripke2014-gene.10kb.squared.protein_coding_only.out.txt",
# "Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/mental_SCZ_Ripke2014.gwassumstats.rolypoly_fmt.tab.gz --gwas_name mental_SCZ_Ripke2014 --run_name gene.10kb.squared.all_genes --window_position gene --window_size_kb 10 --pos_transformation square --n_cores 15 &> log.mental_SCZ_Ripke2014-gene.10kb.squared.all_genes.out.txt",
# "Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/mental_SCZ_Ripke2014.gwassumstats.rolypoly_fmt.tab.gz --gwas_name mental_SCZ_Ripke2014 --run_name tss.10kb.squared.all_genes --window_position tss --window_size_kb 10 --pos_transformation square --n_cores 15 &> log.mental_SCZ_Ripke2014-tss.10kb.squared.all_genes.out.txt",
# "Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/mental_SCZ_Ripke2014.gwassumstats.rolypoly_fmt.tab.gz --gwas_name mental_SCZ_Ripke2014 --run_name tss.20kb.squared.all_genes --window_position tss --window_size_kb 20 --pos_transformation square --n_cores 15 &> log.mental_SCZ_Ripke2014-tss.20kb.squared.all_genes.out.txt",
# "Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/mental_SCZ_Ripke2014.gwassumstats.rolypoly_fmt.tab.gz --gwas_name mental_SCZ_Ripke2014 --run_name tss.10kb.pos_only.all_genes --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_cores 15 &> log.mental_SCZ_Ripke2014-tss.10kb.pos_only.all_genes.out.txt",
# "Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/mental_SCZ_Ripke2014.gwassumstats.rolypoly_fmt.tab.gz --gwas_name mental_SCZ_Ripke2014 --run_name tss.10kb.abs.all_genes --window_position tss --window_size_kb 10 --pos_transformation abs --n_cores 15 &> log.mental_SCZ_Ripke2014-tss.10kb.abs.all_genes.out.txt"
# ]

# ### XXXX
list_cmds = [
"Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/{gwas_name}.gwassumstats.rolypoly_fmt.tab.gz --gwas_name {gwas_name} --run_name {run_name} --outdir /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2.pos_only --window_position tss --window_size_kb 10 --pos_transformation none --n_cores 5 &> log.{gwas_name}-{run_name}.out.txt".format(gwas_name="body_BMI_Locke2015", run_name="tss.10kb.none.all_genes"),
"Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/{gwas_name}.gwassumstats.rolypoly_fmt.tab.gz --gwas_name {gwas_name} --run_name {run_name} --outdir /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2.pos_only --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_cores 5 &> log.{gwas_name}-{run_name}.out.txt".format(gwas_name="body_BMI_Locke2015", run_name="tss.10kb.pos_only.all_genes"),
"Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/{gwas_name}.gwassumstats.rolypoly_fmt.tab.gz --gwas_name {gwas_name} --run_name {run_name} --outdir /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2.pos_only --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_cores 5 &> log.{gwas_name}-{run_name}.out.txt".format(gwas_name="mental_SCZ_Ripke2014", run_name="tss.10kb.pos_only.all_genes"),
"Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/{gwas_name}.gwassumstats.rolypoly_fmt.tab.gz --gwas_name {gwas_name} --run_name {run_name} --outdir /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2.pos_only --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_cores 5 &> log.{gwas_name}-{run_name}.out.txt".format(gwas_name="disease_T2D", run_name="tss.10kb.pos_only.all_genes"),
"Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/{gwas_name}.gwassumstats.rolypoly_fmt.tab.gz --gwas_name {gwas_name} --run_name {run_name} --outdir /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2.pos_only --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_cores 5 &> log.{gwas_name}-{run_name}.out.txt".format(gwas_name="body_WHR_Shungin2015", run_name="tss.10kb.pos_only.all_genes"),
"Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/{gwas_name}.gwassumstats.rolypoly_fmt.tab.gz --gwas_name {gwas_name} --run_name {run_name} --outdir /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2.pos_only --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_cores 5 &> log.{gwas_name}-{run_name}.out.txt".format(gwas_name="lipids_TC", run_name="tss.10kb.pos_only.all_genes"),
"Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/{gwas_name}.gwassumstats.rolypoly_fmt.tab.gz --gwas_name {gwas_name} --run_name {run_name} --outdir /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2.pos_only --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_cores 5 &> log.{gwas_name}-{run_name}.out.txt".format(gwas_name="blood_EOSINOPHIL_COUNT", run_name="tss.10kb.pos_only.all_genes"),
"Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_BMI_Locke2015 --run_name {run_name} --outdir /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2.pos_only --window_position tss --window_size_kb 100 --pos_transformation pos_only --n_cores 3 &> log.body_BMI_Locke2015-{run_name}.out.txt".format(run_name="tss.100kb.pos_only.all_genes"),
"Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_BMI_Locke2015 --run_name {run_name} --outdir /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2.pos_only --window_position tss --window_size_kb 500 --pos_transformation pos_only --n_cores 3 &> log.body_BMI_Locke2015-{run_name}.out.txt".format(run_name="tss.100kb.pos_only.all_genes")
]

 
### with --gwas_linked_file
# "Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/{gwas_name}.gwassumstats.rolypoly_fmt.tab.gz --gwas_name {gwas_name} --run_name {run_name} --outdir /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2.pos_only --gwas_linked_file /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2/rolypoly_objs.{gwas_name}.tss.10kb.square.all_genes.gwas_linked.RData --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_cores 5 &> log.{gwas_name}-{run_name}.out.txt".format(gwas_name="body_BMI_Locke2015", run_name="tss.10kb.pos_only.all_genes"),
### ALL WRITTEN OUT
# Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_BMI_Locke2015 --run_name tss.10kb.none.all_genes --outdir /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2.pos_only --gwas_linked_file /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2/rolypoly_objs.body_BMI_Locke2015.tss.10kb.square.all_genes.gwas_linked.RData --window_position tss --window_size_kb 10 --pos_transformation none --n_cores 5 &> log.body_BMI_Locke2015-tss.10kb.none.all_genes.out.txt
# Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_BMI_Locke2015 --run_name tss.10kb.pos_only.all_genes --outdir /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2.pos_only --gwas_linked_file /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2/rolypoly_objs.body_BMI_Locke2015.tss.10kb.square.all_genes.gwas_linked.RData --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_cores 5 &> log.body_BMI_Locke2015-tss.10kb.pos_only.all_genes.out.txt
# Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/mental_SCZ_Ripke2014.gwassumstats.rolypoly_fmt.tab.gz --gwas_name mental_SCZ_Ripke2014 --run_name tss.10kb.pos_only.all_genes --outdir /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2.pos_only --gwas_linked_file /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2/rolypoly_objs.mental_SCZ_Ripke2014.tss.10kb.square.all_genes.gwas_linked.RData --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_cores 5 &> log.mental_SCZ_Ripke2014-tss.10kb.pos_only.all_genes.out.txt
# Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/disease_T2D.gwassumstats.rolypoly_fmt.tab.gz --gwas_name disease_T2D --run_name tss.10kb.pos_only.all_genes --outdir /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2.pos_only --gwas_linked_file /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2/rolypoly_objs.disease_T2D.tss.10kb.square.all_genes.gwas_linked.RData --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_cores 5 &> log.disease_T2D-tss.10kb.pos_only.all_genes.out.txt
# Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_WHR_Shungin2015.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_WHR_Shungin2015 --run_name tss.10kb.pos_only.all_genes --outdir /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2.pos_only --gwas_linked_file /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2/rolypoly_objs.body_WHR_Shungin2015.tss.10kb.square.all_genes.gwas_linked.RData --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_cores 5 &> log.body_WHR_Shungin2015-tss.10kb.pos_only.all_genes.out.txt
# Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/lipids_TC.gwassumstats.rolypoly_fmt.tab.gz --gwas_name lipids_TC --run_name tss.10kb.pos_only.all_genes --outdir /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2.pos_only --gwas_linked_file /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2/rolypoly_objs.lipids_TC.tss.10kb.square.all_genes.gwas_linked.RData --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_cores 5 &> log.lipids_TC-tss.10kb.pos_only.all_genes.out.txt
# Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/blood_EOSINOPHIL_COUNT.gwassumstats.rolypoly_fmt.tab.gz --gwas_name blood_EOSINOPHIL_COUNT --run_name tss.10kb.pos_only.all_genes --outdir /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2.pos_only --gwas_linked_file /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2/rolypoly_objs.blood_EOSINOPHIL_COUNT.tss.10kb.square.all_genes.gwas_linked.RData --window_position tss --window_size_kb 10 --pos_transformation pos_only --n_cores 5 &> log.blood_EOSINOPHIL_COUNT-tss.10kb.pos_only.all_genes.out.txt
# Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_BMI_Locke2015 --run_name tss.100kb.pos_only.all_genes --outdir /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2.pos_only --window_position tss --window_size_kb 100 --pos_transformation pos_only --n_cores 3 &> log.body_BMI_Locke2015-tss.100kb.pos_only.all_genes.out.txt
# Rscript run_rolypoly_multi_expr_datasets-v2.R --gwas_file /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz --gwas_name body_BMI_Locke2015 --run_name tss.100kb.pos_only.all_genes --outdir /raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2.pos_only --window_position tss --window_size_kb 500 --pos_transformation pos_only --n_cores 3 &> log.body_BMI_Locke2015-tss.100kb.pos_only.all_genes.out.txt



# list_cmds = ["sleep 10"] * 10 # ['sleep 3', 'sleep 3', 'sleep 3']

script2run = "/raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/run_rolypoly_multi_expr_datasets-v2.R"


list_of_processes = []
batch = 1
for i, cmd in enumerate(list_cmds, start=1):
    print "batch = {} | i = {} | Running command: {}".format(batch, i, cmd)
    print cmd
    p = subprocess.Popen(cmd, shell=True)
    list_of_processes.append(p)
    print "batch = {} | i = {} | list_of_processes:".format(batch, i)
    print list_of_processes
    if i % 10 == 0:
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