#!/usr/bin/python


import sys
import os
import random
from os import listdir
from os.path import isfile, join
from joblib import Parallel, delayed

inputdir = "/Bmo/data/projects/human_gut_meta_NCBI/"
outputdir = "/Iceking/dantipov/human_gut/other_datasets/"
#inputdir = sys.argv[1]
sra_pref = "/home/dantipov/other_tools/sratoolkit.2.10.7-ubuntu64/bin/"
other_datasets = "/Bmo/dantipov/data/500_random_datasets/"
list = "/home/dantipov/scripts/human_gut_virome/all_updated.list"
suff = ["_1.fastq.gz", ".fastq.gz", "_2.fastq.gz"]
uncompressed_suff = ["_1.fastq", "_2.fastq", ".fastq"]
def download_sample(id, outdir):
    for ready_dir in [other_datasets, outdir]:
        for cur_suff in suff:
            if isfile(join(ready_dir, id + cur_suff)):
                print id + " found "
                return
    pr_line = sra_pref + "prefetch --max-size 40000000 " + id
    print pr_line
    os.system(pr_line)
    fq_dump_line = sra_pref + "fasterq-dump -t /Bmo/dantipov/tmp/ " + id + " -O " + outdir
    
    print fq_dump_line
    os.system(fq_dump_line)
    for s in uncompressed_suff:
        
        file = join (outdir, id + s)
        if isfile(file):
            os.system ("gzip " + file)

def cleanup (inputlist, fasta_dir):
    for line in open(inputlist,"r"):
        rand = random.randint(0, 10)
#        if rand != 0:
#            continue
        id = line.split()[0]
        strat = line.split()[1]
        tech = line.split()[2]
        if isfile(join(fasta_dir, id+".sra")) and isfile(join(outputdir, id + "_1.fastq.gz")):
            os.remove(join(fasta_dir, id+".sra"))

def process_list(inputlist, outdir ):
    random.seed(239)
    ids = []
    for line in open(inputlist,"r"):
        rand = random.randint(0, 10)
#        if rand != 0:
#            continue
        id = line.split()[0]
        strat = line.split()[1]
        tech = line.split()[2]
        if strat != "WGS" or tech != "ILLUMINA":
            print ("Line is wrong " + line.strip())
            continue
        ids.append(id)    
    Parallel(n_jobs=5)(delayed(download_sample)(id, outdir)
    for id  in ids)
#        exit()
def create_download_list(inputdir):
    names = []
    for f in listdir(inputdir):
        arr = f.split('_')
        if len(arr) < 3:
            continue
        if arr[2][0:3] == "SRR" or arr[2][0:3] == "ERR":
            names.append([arr[2],f])
    random.shuffle(names)
    for name in names:
#        print name[0]+"\t" + name[1]
        print name[0]    

#process_list(sys.argv[1], sys.argv[2])
#create_download_list(inputdir)
process_list(list, outputdir)
#cleanup(list, "/Iceking/dantipov/ncbi/sra/")
#create_download_list(inputdir)
