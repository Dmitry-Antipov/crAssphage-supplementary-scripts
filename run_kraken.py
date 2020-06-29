#!/usr/bin/python3


import sys
import os
import random
from os import listdir
from os.path import isfile, join

from joblib import Parallel, delayed
import csv
import subprocess
import taxonomy_scripts

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
braken_bin = "/Nancy/mrayko/Libs/Bracken-2.5/bracken"
kraken_dir = "/home/dantipov/other_tools/kraken2/kraken/"

kraken_db ="/Bmo/dantipov/gut_pipeline/kraken_viral_db/"
sra_tools =  "/home/dantipov/other_tools/sratoolkit.2.10.7-ubuntu64/bin/"
list = "/home/dantipov/scripts/human_gut_virome/sra_only.list"
patched_list = "/home/dantipov/scripts/human_gut_virome/500_patched.list"
length_list = "/Bmo/dantipov/gut_pipeline/june_abund/all_genomes.length"
inputdir = "/Bmo/dantipov/data/500_random_datasets/"
#classified = "/Bmo/dantipov/gut_pipeline/abundancy_check/596_crass_related_gut_contigs.tsv"
classified = "/Bmo/dantipov/gut_pipeline/june_abund/table_1.tsv"


class bracken_stats:
    def __init__ (self, str):
#100.00  443943  0       D       10239     Viruses
        arr = str.strip().split('\t')
        self.percentage = float(arr[0])
        self.reads_clade = int(arr[1])
        self.reads_self = int(arr[2])
        self.level = arr[3]
        self.id = int(arr[4])
        self.name = arr[5]

def run_sample(sample_descr):
    srr = sample_descr[0]
#    if srr != "SRR11771929":
#        return
    gisaid = sample_descr[-1]
    workdir = join(work_pref, srr)
    if not os.path.isdir (workdir):
        os.mkdir(workdir)
    if os.path.exists(join(workdir, srr+".report")):
        print (srr + " exists")
        return
    
    res = extract_nextstrain_id(srr, gisaid, workdir)
    if res != 0:       
        print ("reference extraction failed " + srr)
        return
    print ("Processing " + srr)
    process_sample(sample_descr, workdir)


def get_kraken_str(srr_id, inputdir, workdir):
    read1 = join(inputdir, srr_id+ "_1.fastq.gz")
    read2 = join(inputdir, srr_id+ "_2.fastq.gz")
#./kraken/kraken2 /Bmo/dantipov/data/500_random_datasets//ERR688506_1.fastq.gz /Bmo/dantipov/data/500_random_datasets//ERR688506_2.fastq
#.gz --threads 20 --db /Bmo/dantipov/gut_pipeline/kraken_viral_db/ --report ERR688506_upd.report>/dev/null
    res = kraken_bin + " " + read1 + " " + read2  + " --threads 15 --db " + kraken_db + " --report " + join(workdir,srr_id+".report") + " >/dev/null"
    return res


def get_bracken_str(srr_id, length, workdir):
#-d /Bmo/dantipov/gut_pipeline/kraken_viral_db/  -i  ERR688506_upd.report -o ERR688506_nodes.bracken -r 100
    infile = join(workdir, srr_id+".report")
    outfile = join(workdir, srr_id+"_bracken.nodes")
    if not os.path.exists(infile):
        print (srr_id + " not found")
        return ""

    if os.path.exists(outfile):
        print (srr_id + " processed")
        return ""
    res =[]
    if not os.path.exists(join(kraken_db, "database{}mers.kmer_distrib".format(length))):
#/Bmo/dantipov/gut_pipeline/kraken_viral_db/database101mers.kmer_distrib
        res.append(braken_bin+ "-build  -d " + kraken_db + " -t 20  -l " + length + " -x " + kraken_dir + " > bracken.log")
    else:
        print  ("db for read length {} already constructed, skipping".format(length))
    res.append( braken_bin + " -d " + kraken_db + " -i " +infile + " -o " + outfile + " -r " + length + " > bracken.log" )
#./bracken-build -d ${KRAKEN_DB} -t ${THREADS} -k ${KMER_LEN} -l ${READ_LEN} -x ${KRAKEN_INSTALLATION}

    return res

def bracken_sample(srr_id, length, workdir):
    bracken_str = get_bracken_str(srr_id, length, workdir)
    if bracken_str != "":
        print (bracken_str)
        for line in bracken_str:
            os.system(line)

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
    print (kraken_str)
    os.system(kraken_str)

def run_all_kraken(list, inputdir, workdir):
    ids = []
    for line in open (list, "r"):
        ids.append(line.strip())
    Parallel(n_jobs=5)(delayed(kraken_sample)(inputdir, id, workdir)
    for id  in ids)

def run_all_bracken (list, workdir):
    ids = []
    for line in open (list, "r"):
        arr = line.strip().split()
        if  not arr[1].isdecimal():
            continue
        bracken_sample(arr[0], arr[1], workdir)        

def merge_brackens(workdir):
    nodes = taxonomy_scripts.read_nodes(taxonomy_scripts.kraken_names, taxonomy_scripts.kraken_nodes)
    bracken_all = {}
    for f in workdir:
        if f.split("_")[1] == "bracken.report":
            for line in open(join(workdir,f),'r'):
                br = bracken_stats(line)
                if not br.id in bracken_all:
                    bracken_all[br.id] = br
                else:
                    bracken_all[br.id].reads_clade += br.reads_clade
                    bracken_all[br.id].reads_self += br.reads_self
    for br in bracken_all:
        bracken_all[br].percentage = 100.00 * bracken_all[br].reads_clade / bracken_all[1].reads_clade
#merge_brackens(workdir)
run_all_bracken(patched_list, workdir)
#run_all(list, inputdir, workdir)

#count_absolute(workdir, length_list)
