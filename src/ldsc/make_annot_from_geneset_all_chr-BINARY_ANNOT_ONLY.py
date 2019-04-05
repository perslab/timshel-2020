#!/usr/bin/env python2.7
from __future__ import print_function
import pandas as pd
import numpy as np
import argparse
from pybedtools import BedTool
# import gzip


### PT adds
import multiprocessing
import os
import sys
import subprocess

import pdb

###################################### USAGE ######################################
# Must run on python2.7, because of issue with pybedtools in python3.

### test
# python2 make_annot_from_geneset_all_chr.py \
# --file_multi_gene_set /projects/timshel/sc-genetics/ldsc/ldsc/test_file_multi_gene_set_wgcna200.csv \
# --file_gene_coord /projects/timshel/sc-genetics/ldsc/data/gene_coords/gene_annotation.hsapiens_all_genes.GRCh37.ens_v91.LDSC_fmt.txt \
# --windowsize 100000 \
# --bimfile_basename /projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC \
# --out_dir ./tmp_tmp_test_make_annot \
# --out_prefix test_xxx1 \
# --flag_wgcna \
# --flag_mouse


### lira sema (scratch)
# time python2 make_annot_from_geneset_all_chr.py \
# --file_multi_gene_set /projects/timshel/sc-genetics/sc-genetics/out/out.wgcna/nn_lira_sema/tables/nn_lira_sema_per_brain_area_run1_cell_cluster_module_genes.csv \
# --file_gene_coord /projects/timshel/sc-genetics/ldsc/data/gene_coords/gene_annotation.hsapiens_all_genes.GRCh37.ens_v91.LDSC_fmt.txt \
# --windowsize 100000 \
# --bimfile_basename /projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC \
# --out_dir /scratch/sc-ldsc/nn_lira_sema \
# --out_prefix nn_lira_sema \
# --flag_wgcna \
# --flag_mouse


### MACA
# time python2 make_annot_from_geneset_all_chr.py \
# --file_multi_gene_set /projects/jonatan/tmp-maca/tables/maca_tissue_cell_type_kME_cell_cluster_module_genes.csv \
# --file_gene_coord /projects/timshel/sc-genetics/ldsc/data/gene_coords/gene_annotation.hsapiens_all_genes.GRCh37.ens_v91.LDSC_fmt.txt \
# --windowsize 100000 \
# --bimfile_basename /projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC \
# --out_dir /scratch/sc-ldsc/maca \
# --out_prefix maca_tissue_cell_type \
# --flag_wgcna \
# --flag_mouse

### Mousebrain Neurons_ClusterName (n=3833)
### 11 processes takes ~500 GB
# time python2 make_annot_from_geneset_all_chr.py \
# --file_multi_gene_set /projects/jonatan/tmp-mousebrain/tables/mousebrain_Neurons_ClusterName_2_cell_cluster_module_genes.csv \
# --file_gene_coord /projects/timshel/sc-genetics/ldsc/data/gene_coords/gene_annotation.hsapiens_all_genes.GRCh37.ens_v91.LDSC_fmt.txt \
# --windowsize 100000 \
# --bimfile_basename /projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC \
# --out_dir /scratch/sc-ldsc/mousebrain \
# --out_prefix Neurons_ClusterName \
# --flag_wgcna \
# --flag_mouse \
# --n_parallel_jobs 3


### Mette thesis hypothalamus
# time python2 make_annot_from_geneset_all_chr.py \
# --file_multi_gene_set /projects/timshel/sc-genetics/sc-genetics/data/gene_lists/mludwig_thesis_hypothalamus_wgcna_modules.csv \
# --file_gene_coord /projects/timshel/sc-genetics/ldsc/data/gene_coords/gene_annotation.hsapiens_all_genes.GRCh37.ens_v91.LDSC_fmt.txt \
# --windowsize 100000 \
# --bimfile_basename /projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC \
# --out_dir /scratch/sc-ldsc/hypothalamus_mette_thesis \
# --out_prefix hypothalamus_mette_thesis \
# --flag_wgcna \
# --flag_mouse


### Mousebrain Ependymal_ClusterName (n=XXX)
# time python2 make_annot_from_geneset_all_chr.py \
# --file_multi_gene_set /projects/jonatan/tmp-mousebrain/tables/mousebrain_Ependymal_ClusterName_2_cell_cluster_module_genes.csv \
# --file_gene_coord /projects/timshel/sc-genetics/ldsc/data/gene_coords/gene_annotation.hsapiens_all_genes.GRCh37.ens_v91.LDSC_fmt.txt \
# --windowsize 100000 \
# --bimfile_basename /projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC \
# --out_dir /scratch/sc-ldsc/mousebrain \
# --out_prefix Ependymal_ClusterName \
# --flag_wgcna \
# --flag_mouse \
# --n_parallel_jobs 3


###################################### TMP ANNOTBED - delete ######################################

#  (Pdb) bed_for_annot.head(1000) ---> len=~66
# chr1    51619935        52185000
#  chr1   70410488        70871303
#  chr1   85584164        86243933
#  chr1   202948059       203355877
#  chr10  43851792        44270066
#  chr10  75681524        76110821
#  chr10  76769912        77191206
#  chr10  120663598       121138345
#  chr11  118030300       118469926
#  chr12  21454715        21871342
#  chr13  47145391        47571367
#  chr14  23646017        24048981
#  chr14  62253803        62768431
#  chr14  75545477        76327532
#  chr14  77693018        78124295
#  chr14  90662846        91074605
#  chr15  73652355        74126475
#  chr17  11724141        12247147
#  chr17  16950141        17384607
#  chr17  29399031        29824557
#  chr17  33133009        33616338
#  chr17  43038067        43447407
#  chr17  45818872        46225654
#  chr18  40647843        41057615
#  chr2   9524101 9971143
#  chr2   30169807        30583399
#  chr2   73856086        74300786
#  chr2   172440880       173064766
#  chr20  17722241        18149623
#  chr20  42924862        43350750
#  chr21  34604792        35052318
#  chr3   5029331 5461642
#  chr3   44756749        45217677
#  chr3   49260379        49666759
#  chr3   112080556       112504424
#  chr3   118987785       119413555
#  chr3   150059781       150521015
#  chr3   196473214       196895931
#  chr4   157797209       158293242
#  chr4   166048775       166464312
#  chr5   61499799        62124409
#  chr5   133506870       133927683
#  chr5   141137893       141569856
#  chr5   146570374       147089619
#  chr5   149888002       150338671
#  chr5   158490089       158913044
#  chr6   88099839        88577169
#  chr6   99116420        99595849
#  chr7   38562563        39171994
#  chr7   129946089       130553598
#  chr7   154889486       155301945
#  chr8   11453082        11896818
#  chr8   41147915        41568499
#  chr8   42795556        43257998
#  chr8   80323049        81343467
#  chr8   96057147        96481429
#  chr8   98587285        99065241
#  chr8   101069288       101548446
#  chr8   130651839       131229375
#  chr9   17379080        17997127
#  chr9   81986688        82541658
#  chr9   123377774       123805262
#  chrX   9925024 10405700
#  chrX   77185245        77595203
#  chrX   99683667        100094988
#  chrX   119361682       119803220

###################################### FILE SNIPPETS ######################################


# @@@@@@@@@@@@@@@@@@ file_multi_gene_set @@@@@@@@@@@@@@@@@@

### WGCNA input (flag_wgcna = True). 
### We use the 'module' column as annotation name, since they should be unique
# cell_cluster,module,ensembl,hgcn,pkME
# Aorta_endothelial cell,antiquewhite3,ENSMUSG00000017639,Rab11fip4,0.888954070670456
# Aorta_endothelial cell,antiquewhite3,ENSMUSG00000004233,Wars2,0.860974901276187
# Aorta_endothelial cell,antiquewhite3,ENSMUSG00000031487,Brf2,0.855919687669775
# Aorta_endothelial cell,antiquewhite3,ENSMUSG00000032997,Chpf,0.824264829666081

### non WGCNA input.
### NO HEADER. Delim=csv. 
### Col1=annotation_name, Col2=EnsemblID (mouse or human)
# brain_cortex,ENSMUSG00000017639
# brain_cortex,ENSMUSG00000004233
# brain_cortex,ENSMUSG00000031487
# brain_cortex,ENSMUSG00000032997
# ...
# brain_hypothlamus,ENSMUSG00000032997

def map_ensembl_genes_mouse_to_human(df):
    """
    DESCRIPTION
        function to map mouse ensembl genes to human ortholog genes.
    INPUT
       df: a dataframe with two columns: "annotation" and "gene_input". "gene_input" column should contain Ensembl mouse gene_input IDs to be mapped.
    OUTPUT
       df: input dataframe with mapped genes. Mouse genes that could not be mapped are removed.
       file_mapping_summary: a summary file with mapping stats

    REMARKS
        We assume the file_mapping contains only 1-1 mapping (which is true for gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz).
        Otherwise the .map() function might fail
    """
    file_out_mapping_stats = "{}/log.{}.make_annotation_mapping_stats.txt".format(out_dir, out_prefix)
    
    file_mapping = "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz"
    ### SNIPPET
    # ensembl_gene_id chromosome_name start_position  end_position    mmusculus_homolog_ensembl_gene  mmusculus_homolog_orthology_confidence
    # ENSG00000138593 15      49280673        49338760        ENSMUSG00000035093      1
    # ENSG00000166351 21      14982498        15013906        ENSMUSG00000095294      0
    # ENSG00000168675 18      13217497        13652754        ENSMUSG00000024544      1
    
    df_mapping = pd.read_csv(file_mapping, delim_whitespace=True, usecols=["ensembl_gene_id", "mmusculus_homolog_ensembl_gene"], index_col=1) # index is MOUSE genes
    
    ### map ortholog genes
    df["gene"] = df["gene_input"].map(df_mapping["ensembl_gene_id"]) # df["gene_input"] is mouse genes. df_mapping["ensembl_gene_id"] is human.
    # ^ .map() returns NaN values for genes not mapped. 

    ### make summary of mapping
    df_summary = df.groupby("annotation")["gene"].agg({'n_genes_input': lambda x: len(x),
                                                         'n_genes_output': lambda x: len(x)-sum(pd.isnull(x)), 
                                                        'n_genes_not_mapped' : lambda x: sum(pd.isnull(x)),
                                                        'pct_genes_not_mapped': lambda x: "{:.2f}".format(sum(pd.isnull(x))/float(len(x))*100)})
    # ^ we use pd.isnull() instead of np.isnan() because of the issue described here: https://stackoverflow.com/a/36001191/6639640
    df_summary.sort_values(by=['n_genes_output'], inplace=True) 
    df_summary.to_csv(file_out_mapping_stats, sep="\t")
    
    ### final processing
    df.dropna(axis=0, inplace=True) # remove non-mapped genes
    return df


def read_multi_gene_set_file(file_multi_gene_set, flag_wgcna, flag_mouse):
    if flag_wgcna:
        df_multi_gene_set = pd.read_csv(file_multi_gene_set) 
        df_multi_gene_set.rename(columns={"module":"annotation", "ensembl":"gene_input"}, inplace=True) # df.rename(columns={'oldName1': 'newName1', 'oldName2': 'newName2'})
        # ^ 'module' = first column; will later be renamed to 'annotation'
        # ^ 'ensembl' = second column; will later be renamed to 'gene'
    else:
        df_multi_gene_set = pd.read_csv(file_multi_gene_set, header=None)
        df_multi_gene_set.columns = ["annotation", "gene_input"] # setting or renaming column names
    print("========================== STATS file_multi_gene_set ====================")
    print("Number of gene sets: {}".format(df_multi_gene_set["annotation"].nunique()))
    print("=========================================================================")
    if flag_mouse:
        df_multi_gene_set = map_ensembl_genes_mouse_to_human(df_multi_gene_set) # adds "gene" column
    else:
        df_multi_gene_set["gene"] = df_multi_gene_set["gene_input"] # copy
    ### write_multi_geneset_file
    file_out_multi_geneset = "{}/log.{}.multi_geneset.txt".format(out_dir, out_prefix)
    df_multi_gene_set.to_csv(file_out_multi_geneset, sep="\t", index=False)
    return df_multi_gene_set


def multi_gene_sets_to_dict_of_beds(df_multi_gene_set, df_gene_coord, windowsize):
    """ 
    INPUT
        df_multi_gene_set: two columns "annotation" and "gene". Gene is human Ensembl gene names.
    OUTPUT
        dict_of_beds: returns a dict of beds. Keys are annotation names from df_multi_gene_set.
    """
    print('making gene set bed files')
    dict_of_beds = {}
    for name_annotation, df_group in df_multi_gene_set.groupby("annotation"):
        print(name_annotation)
        df = pd.merge(df_gene_coord, df_group, left_on="GENE", right_on="gene", how = "inner") 
        df['START'] = np.maximum(0, df['START'] - windowsize)
        df['END'] = df['END'] + windowsize
        iter_df = [['chr'+(str(x1).lstrip('chr')), x2, x3] for (x1,x2,x3) in np.array(df[['CHR', 'START', 'END']])] # notice that only 3 columns from df is used.
        bed_for_annot = BedTool(iter_df).sort().merge() # PT NOTE: .merge(): Merge overlapping features together. https://daler.github.io/pybedtools/autodocs/pybedtools.bedtool.BedTool.merge.html#pybedtools.bedtool.BedTool.merge
        dict_of_beds[name_annotation] = bed_for_annot
    return dict_of_beds


def get_annot_file_path(chromosome):
    file_out_annot_combined = "{}/{}.{}.{}.annot.gz".format(out_dir, out_prefix, "COMBINED_ANNOT", chromosome) # set output filename. ${prefix}.${chr}.annot.gz
    return file_out_annot_combined
    


def check_annot_file(chromosome):
    ### check for existing annot file (and it's integrity)
    file_out_annot_combined = get_annot_file_path(chromosome)
    if os.path.exists(file_out_annot_combined):
        print("CHR={} | file_out_annot_combined exists: {}. Will check for integrity...".format(chromosome, file_out_annot_combined))
        cmd_gzip_test = "gzip -t {}".format(file_out_annot_combined)
        with open(os.devnull, 'w') as fnull:
            p = subprocess.Popen(cmd_gzip_test, shell=True, stdout=fnull, stderr=subprocess.STDOUT)
            # gzip -t doesn't have any output, other than the return code, if it's a correct gzip compressed file. But if it is a corrupted file, it will output something:
            # CORRUPT FILE --> gzip -t Neurons_ClusterName.COMBINED_ANNOT.17.annot.gz --> "gzip: Neurons_ClusterName.COMBINED_ANNOT.17.annot.gz: unexpected end of file"
            p.wait()
        if p.returncode == 0: # gz file ok
            print("CHR={} | file_out_annot_combined integrity OK: {}. Will not make new file".format(chromosome, file_out_annot_combined))
            return None
        else: # file not ok.
            return chromosome
    else: # file does not exist
        return chromosome



def make_annot_file_per_chromosome(chromosome):
    """ 
    Input
        chromosome: integer (1..22)
    
    *OBS* this function RELIES on MANY GLOBAL scope VARIABLES
    """
    # TODO: parse variables to function
    
    ### make annot file
    print('making annot files for chromosome {}'.format(chromosome))
    bimfile = "{}.{}.bim".format(args.bimfile_basename, chromosome) # e.g. <LONGPATH/1000G.EUR.QC>.15.bim for chromosome 15
    df_bim = pd.read_csv(bimfile, delim_whitespace=True, usecols = [0,1,2,3], names = ['CHR','SNP','CM','BP'])
    # (Pdb) df_bim.head()
    #    CHR          SNP        CM       BP
    # 0   21  rs146134162 -0.908263  9412099
    # 1   21  rs578050168 -0.908090  9412377
    # 2   21  rs527616997 -0.907297  9413645
    # 3   21  rs544748596 -0.906578  9414796
    # 4   21  rs528236937 -0.906500  9414921
    iter_bim = [['chr'+str(x1), x2, x2] for (x1, x2) in np.array(df_bim[['CHR', 'BP']])]
    bimbed = BedTool(iter_bim)
    counter = 1 # just to print status message
    list_df_annot = []
    for name_annotation in sorted(dict_of_beds): # we sort to make output more consistent.
        print("CHR={} | annotation={}, #{}/#{}".format(chromosome, name_annotation, counter, len(dict_of_beds)))
        bed_for_annot = dict_of_beds[name_annotation] # get bed
        # (Pdb)len(bed_for_annot)
        # 66
        #  (Pdb) bed_for_annot.head()
        # chr1    51619935        52185000
        #  chr1   70410488        70871303
        #  chr1   85584164        86243933
        #  chr1   202948059       203355877
        #  chr10  43851792        44270066
        #  chr10  75681524        76110821
        #  chr10  76769912        77191206
        #  chr10  120663598       121138345
        #  chr11  118030300       118469926
        #  chr12  21454715        21871342
        annotbed = bimbed.intersect(bed_for_annot) # PT NOTE: this finds SNPs in bim file that OVERLAP with the annotation bed (gene)? 
        # (Pdb) annotbed.head() | when chromosome=21
        # chr21   34605531        34605531
        #  chr21  34605604        34605604
        #  chr21  34605644        34605644
        #  chr21  34605778        34605778
        #  chr21  34606634        34606634
        #  chr21  34606840        34606840
        #  chr21  34607223        34607223
        #  chr21  34607358        34607358
        # SEE https://daler.github.io/pybedtools/intersections.html
        # SEE https://daler.github.io/pybedtools/autodocs/pybedtools.bedtool.BedTool.intersect.html#pybedtools.bedtool.BedTool.intersect
        bp = [x.start for x in annotbed] # PT NOTE: make list of all bp positions for the overlapping SNPs | All features, no matter what the file type, have chrom, start, stop, name, score, and strand attributes.
        df_int = pd.DataFrame({'BP': bp, name_annotation:1}) # FINUCANE ORIG: df_int = pd.DataFrame({'BP': bp, 'ANNOT':1})
        #             BP  blue
        # 0     34605531     1
        # 1     34605604     1
        # 2     34605644     1
        # 3     34605778     1
        # 4     34606634     1
        # 5     34606840     1
        # 6     34607223     1
        df_annot = pd.merge(df_bim, df_int, how='left', on='BP') # merges will ALL snps from the bim file.
        # (Pdb) df_annot.head()
        #    CHR          SNP        CM       BP  blue
        # 0   21  rs146134162 -0.908263  9412099   NaN
        # 1   21  rs578050168 -0.908090  9412377   NaN
        # 2   21  rs527616997 -0.907297  9413645   NaN
        # 3   21  rs544748596 -0.906578  9414796   NaN
        # 4   21  rs528236937 -0.906500  9414921   NaN
        df_annot.fillna(0, inplace=True)
        df_annot = df_annot[[name_annotation]].astype(int) # *OBS*: KEEPS ONLY ONE COLUMN. FINUCANE ORIG: df_annot = df_annot[['ANNOT']].astype(int)
        list_df_annot.append(df_annot)
        if args.flag_annot_file_per_geneset: # write annot file per annotation per chromosome
            file_out_annot = "{}/{}.{}.{}.annot.gz".format(out_dir, out_prefix, name_annotation, chromosome) # set output filename. ${prefix}.${chr}.annot.gz
            df_annot.to_csv(file_out_annot, sep="\t", index=False, compression="gzip")
        counter += 1
    print("CHR={} | Joining and writing annotations...".format(chromosome))
    df_annot_combined = pd.concat(list_df_annot, axis=1)
    file_out_annot_combined = get_annot_file_path(chromosome) 
    df_annot_combined.to_csv(file_out_annot_combined, sep="\t", index=False, compression="gzip")
    return None



###################################### MAIN ######################################

parser = argparse.ArgumentParser()
parser.add_argument('--file_multi_gene_set', type=str, help='a file of gene names, one line per gene.')
parser.add_argument('--file_gene_coord', type=str, help='a file with columns GENE, CHR, START, and END, where START and END are base pair coordinates of TSS and TES. This file can contain more genes than are in the gene set. We provide ENSG_coord.txt as a default.')
parser.add_argument('--windowsize', type=int, help='how many base pairs to add around the transcribed region to make the annotation? Finucane uses 100 kb.')
parser.add_argument('--bimfile_basename', type=str, help='plink bim BASENAME for the dataset you will use to compute LD scores. If argument is "1000G.EUR.QC", then the files "1000G.EUR.QC.1.bim", "1000G.EUR.QC.2.bim", ..., "1000G.EUR.QC.22.bim" will be loaded')
parser.add_argument('--out_dir', type=str, help='output directory to write annot files. Relative or absolute. Dir be created if it does not exist. ')
parser.add_argument('--out_prefix', type=str, help='Prefix for output files. Outputfiles will be <out_dir>/<out_prefix>.<name_annotation>.<chromosome>.annot.gz')
parser.add_argument('--flag_annot_file_per_geneset', action='store_true', help='set flag to write one annot file per gene set per chromosome. Default is to write a combined annot file per chromosome containing all gene sets. NB: the annotation files are always split per chromosome.')
parser.add_argument('--flag_wgcna', action='store_true', help='set flag if file_multi_gene_set input is from WGCNA pipeline')
parser.add_argument('--flag_mouse', action='store_true', help='set flag if ile_multi_gene_set input has mouse genes (instead of human)')
parser.add_argument('--n_parallel_jobs', type=int, default=22, help='Number of processes. Default 22 (the number of chromosomes)')

#parser.add_argument('--annot-file', type=str, help='the name of the annot file to output.')


args = parser.parse_args()

### TODO
# [ ] check that ALL bim files exist. (so pool of workers does not fail)

out_prefix = args.out_prefix
out_dir = args.out_dir


### Make out_dir
if not os.path.exists(out_dir):
    print("Making output dir {}".format(out_dir))
    os.makedirs(out_dir)


list_chromosomes = range(1,23) # 1...22
### Check for existing annot files
pool = multiprocessing.Pool(processes=len(list_chromosomes))
list_chromosomes_to_run = pool.map(check_annot_file, list_chromosomes) # check_annot_file returns None if annot file exists and is ok
list_chromosomes_to_run = [x for x in list_chromosomes_to_run if x is not None] # filter away None

print("========================== CHROMOSOMES TO RUN ====================")
print("N chromosomes = {}".format(len(list_chromosomes_to_run)))
print("Chromosomes: {}".format(",".join(map(str, list_chromosomes_to_run))))
print("=========================================================================")
if len(list_chromosomes_to_run) == 0:
    print("Quitting...")
    sys.exit()

### read coord file
df_gene_coord = pd.read_csv(args.file_gene_coord, delim_whitespace = True)
### read gene list file (a list of gene list files)
df_multi_gene_set = read_multi_gene_set_file(args.file_multi_gene_set, flag_wgcna=args.flag_wgcna, flag_mouse=args.flag_mouse)
### make beds
dict_of_beds = multi_gene_sets_to_dict_of_beds(df_multi_gene_set, df_gene_coord, args.windowsize)

print("Starting pool...")
pool = multiprocessing.Pool(processes=args.n_parallel_jobs)
pool.map(make_annot_file_per_chromosome, list_chromosomes_to_run)
# make_annot_file_per_chromosome(21)

print("Script is done!")



