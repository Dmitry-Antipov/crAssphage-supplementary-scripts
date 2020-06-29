#!/usr/bin/python


import sys
import os
import random
from os import listdir
from os.path import isfile, join

from joblib import Parallel, delayed
import csv
from StringIO import StringIO
import subprocess


#constants
low_coverage = 10
samtools = "samtools"
bcftools = "~/other_tools/bcftools/bcftools"
minimap2 = "~/other_tools/minimap2/minimap2"
ref = "/Bmo/dantipov/gut_pipeline/june_abund/all_genomes.fa"
#ref = "/Bmo/dantipov/metagenomes_Koonin/study/crass/crassphage.fasta"
#workdir = "/Iceking/dantipov/human_gut/depth_reports/"
#workdir = "/Iceking/dantipov/human_gut/crass_reports/"
#workdir = "/Iceking/dantipov/human_gut/check/"
workdir = "/Bmo/dantipov/gut_pipeline/june_abund/kraken_res/"
kraken_bin = "/home/dantipov/other_tools/kraken2/kraken/kraken2"
kraken_db ="/Bmo/dantipov/gut_pipeline/kraken_viral_db/"
sra_tools =  "/home/dantipov/other_tools/sratoolkit.2.10.7-ubuntu64/bin/"
list = "/home/dantipov/scripts/human_gut_virome/sra_only.list"
length_list = "/Bmo/dantipov/gut_pipeline/june_abund/all_genomes.length"
inputdir = "/Bmo/dantipov/data/500_random_datasets/"
#classified = "/Bmo/dantipov/gut_pipeline/abundancy_check/596_crass_related_gut_contigs.tsv"
classified = "/Bmo/dantipov/gut_pipeline/june_abund/table_1.tsv"

class reference_stats:
    def __init__(self, name):
        self.name = name
        self.samples = 0
        self.total_coverage = 0
        self.length = 0
        self.additional = ""
        self.intid = 0

def run_sample(sample_descr):
    srr = sample_descr[0]
#    if srr != "SRR11771929":
#        return
    gisaid = sample_descr[-1]
    workdir = join(work_pref, srr)
    if not os.path.isdir (workdir):
        os.mkdir(workdir)
    if os.path.exists(join(workdir, srr+".report")):
        print srr + " exists"
        return
    
    res = extract_nextstrain_id(srr, gisaid, workdir)
    if res != 0:       
        print "reference extraction failed " + srr
        return
    print "Processing " + srr
    process_sample(sample_descr, workdir)


def get_kraken_str(srr_id, inputdir, workdir):
    read1 = join(inputdir, srr_id+ "_1.fastq.gz")
    read2 = join(inputdir, srr_id+ "_2.fastq.gz")
#./kraken/kraken2 /Bmo/dantipov/data/500_random_datasets//ERR688506_1.fastq.gz /Bmo/dantipov/data/500_random_datasets//ERR688506_2.fastq
#.gz --threads 20 --db /Bmo/dantipov/gut_pipeline/kraken_viral_db/ --report ERR688506_upd.report>/dev/null
    res = kraken_bin + " " + read1 + " " + read2  + " --threads 15 --db " + kraken_db + " --report " + join(workdir,srr_id+".report") + " >/dev/null"
    return res



def kraken_sample(inputdir, srr_id, workdir):
    if  not os.path.isdir (workdir):
        os.mkdir(workdir)
#    if sample_descr[-2] == "OXFORD_NANOPORE":
#TODO elif ionTORRENT    
#        minimap_str = get_minimap_nanopore(srr_id, workdir)
#    elif sample_descr[2] == "PAIRED":
#        minimap_str = get_minimap_illumina_paired(srr_id, workdir)
#    else:
#        minimap_str = get_minimap_illumina_single(srr_id, workdir)
#    print minimap_str
    if not os.path.isfile (join(inputdir, srr_id + "_1.fastq.gz")):
        print (srr_id + " is not present")
        return
    if not os.path.isfile (join(inputdir, srr_id + "_2.fastq.gz")):
        print (srr_id + " is not paired illumina")
        return
    kraken_str = get_kraken_str(srr_id, inputdir, workdir)
    print kraken_str
    os.system(kraken_str)

def run_all(list, inputdir, workdir):
    ids = []
    for line in open (list, "r"):
        ids.append(line.strip())
#    for id in ids:
#        kraken_sample(inputdir, id, workdir)

    Parallel(n_jobs=5)(delayed(kraken_sample)(inputdir, id, workdir)
    for id  in ids)

def count_absolute(workdir, length_list):
    ref_len = {}
    total_amount = {}
    stats = {}
    ordered = []
    intid = 1
    for line in open(length_list, 'r'):
        arr = line.split()
        ordered.append(arr[0])
        stats[arr[0]] = reference_stats(arr[0])
        stats[arr[0]].length = int(arr[1])
        stats[arr[0]].intid = intid
        intid += 1
    additional_header = ""
    for line in open(classified, "r"):
        arr = line.strip().split('\t')
        if arr[0][0] == "#":
            additional_header =  "\t".join(arr[1:])
        else:
            stats[arr[0]].additional = "\t".join(arr[1:])
    tst = 0
    for f in listdir(workdir):
        if f.split('.')[-1] == "depth":
#            print f
            samples = set()
            for line in open (join(workdir, f), 'r'): 
                arr = line.split()
                name = arr[0]
                if not name in samples:
                    samples.add(name)
                    stats[name].samples += 1
                stats[name].total_coverage += int(arr[2])
            tst +=1
#            if tst == 10:
#                break
    print("id\tname\tsamples_present\ttotal_covered\taverage_depth\t" + additional_header)
    for r in ordered:
        print (str(stats[r].intid) + "\t" + r + "\t" + str(stats[r].samples) + "\t" + str(stats[r].total_coverage) + "\t" + str(stats[r].total_coverage/stats[r].length)+ "\t" + stats[r].additional)
                   
#    for id in ids:
#        map_sample(inputdir, id, workdir)


run_all(list, inputdir, workdir)

#count_absolute(workdir, length_list)
