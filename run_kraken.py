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
workdir = "/Bmo/dantipov/gut_pipeline/june_abund/kraken_res2/"
kraken_bin = "/home/dantipov/other_tools/kraken2/kraken/kraken2"
braken_bin = "/Nancy/mrayko/Libs/Bracken-2.5/bracken"
kraken_dir = "/home/dantipov/other_tools/kraken2/kraken/"
tmp_dir = "/Bmo/dantipov/tmp/"
trimmomatic_jar = "/Nancy/mrayko/Libs/Trimmomatic-0.39/trimmomatic-0.39.jar"
adapters_fasta = "/home/dantipov/scripts/human_gut_virome/adapters.fa"
trim_in_suff = ["_1.fastq.gz", ".fastq.gz", "_2.fastq.gz"]
trim_out_suff = ["_1P.fastq", "_2P.fastq", "_1U.fastq", "_2U.fastq", "_S.fastq"]

kraken_db ="/Bmo/dantipov/gut_pipeline/kraken_viral_db/"
sra_tools =  "/home/dantipov/other_tools/sratoolkit.2.10.7-ubuntu64/bin/"
#list = "/home/dantipov/scripts/human_gut_virome/500_updated.list"
#ID Strategy Tech reads_length reads_number
list = "/Bmo/dantipov/data/500_random_datasets/500_updated.list"
#patched_list = "/home/dantipov/scripts/human_gut_virome/all_length.list"
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
        self.kraken_clade = 0
        self.kraken_self = 0
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


def get_kraken_str(srr_id, tmp_dir, workdir):
#    read1 = join(inputdir, srr_id+ "_1.fastq.gz")
#    read2 = join(inputdir, srr_id+ "_2.fastq.gz")
    reads = " "
    for input_suffix in trim_out_suff:
        f = join(tmp_dir, srr_id + input_suffix)
        if isfile (f):
            reads += " " + f

#./kraken/kraken2 /Bmo/dantipov/data/500_random_datasets//ERR688506_1.fastq.gz /Bmo/dantipov/data/500_random_datasets//ERR688506_2.fastq
#.gz --threads 20 --db /Bmo/dantipov/gut_pipeline/kraken_viral_db/ --report ERR688506_upd.report>/dev/null
    res = kraken_bin + " " + reads  + " --threads 15 --db " + kraken_db + " --report " + join(workdir,srr_id+".report") + " >/dev/null"
    return res

def get_trimmomatic_str(srr_id, input_dir, outdir):
    reads = " "
    count = 0
    for input_suffix in trim_in_suff:
        f = join(input_dir, srr_id + input_suffix)
        if isfile (f):
            reads += " " + f
            count += 1
    if count == 2:
        mode = " PE "
        out_opt = " -baseout "+ join(outdir, srr_id)+".fastq" 
    else:
        mode = " SE "
        out_opt = join(outdir, srr_id)+"_S.fastq" 
    str = "java -jar "+ trimmomatic_jar +mode + reads + " -threads 20 " + out_opt  + "  ILLUMINACLIP:" + adapters_fasta + ":2:30:10 >/dev/null"
    return str
#java -jar  /Nancy/mrayko/Libs/Trimmomatic-0.39/trimmomatic-0.39.jar  PE -threads 30 -basein /Bmo/dantipov/data/500_random_datasets/ERR525702_1.fastq.gz  -baseout ERR525702 ILLUMINACLIP:/Nancy/mrayko/Libs/Trimmomatic-0.39/adapters/total.fa:2:30:10

def clear_temporary_trim(srr_id, tmp_dir):
    reads = " "
    for input_suffix in trim_out_suff:
        f = join(tmp_dir, srr_id + input_suffix)
        if isfile (f):
            os.remove(f)
    

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
   # print (join(inputdir, srr_id + "_1.fastq.gz"))
    if os.path.isfile (join(workdir, srr_id +".report")):
        print (srr_id + " is already processed")
        return
    if not os.path.isfile (join(inputdir, srr_id + "_1.fastq.gz")) and not os.path.isfile(join(inputdir, srr_id + ".fastq.gz")) :
        print (srr_id + " is not present")
        return
    trimmomatic_str = get_trimmomatic_str (srr_id, inputdir, tmp_dir)
#    print (trimmomatic_str)
    os.system(trimmomatic_str)
    kraken_str = get_kraken_str(srr_id,tmp_dir, workdir)
#    print (kraken_str)
    os.system(kraken_str)  
    clear_temporary_trim(srr_id, tmp_dir )

def run_all_kraken(list, inputdir, workdir):
    ids = []
    for line in open (list, "r"):
        arr = line.split()
        if len(arr) > 0 and arr[1] == "WGS" and arr[2] == "ILLUMINA":
           ids.append(arr[0])
    Parallel(n_jobs=7)(delayed(kraken_sample)(inputdir, id, workdir)
    for id  in ids)

def run_all_bracken (list, workdir):
    ids = []
    for line in open (list, "r"):
        arr = line.strip().split()
        try:
            length = float(arr[3])
            length = int(round(length))

        except:    
            print (arr[0] + " is not illumina " + arr[1])
            continue
        bracken_sample(arr[0], str(length), workdir)        
        
def merge_brackens(workdir):
    nodes = taxonomy_scripts.read_nodes(taxonomy_scripts.kraken_names, taxonomy_scripts.kraken_nodes)
    bracken_all = {}
    kracken_all = {}
    add_kraken = True
    total_reads = 0
    for f in os.listdir(workdir):
        arr = f.split("_")
        
        if len(arr) > 1 and arr[-1] == "bracken.report":
            tmp  = {}
            for line in open(join(workdir,f),'r'):
                br = bracken_stats(line)
                tmp[br.id] = br
                if not br.id in bracken_all:
                    bracken_all[br.id] = br
                else:
                    bracken_all[br.id].reads_clade += br.reads_clade
                    bracken_all[br.id].reads_self += br.reads_self
            kr_file = arr[0] + ".report"
            for line in open(join(workdir, kr_file),'r'):
                kr = bracken_stats(line)
                if kr.id <= 1 :
                    total_reads += kr.reads_clade
                if kr.id in bracken_all:
                    bracken_all[kr.id].kraken_clade += kr.reads_clade
                    bracken_all[kr.id].kraken_self += kr.reads_self
                    if kr.id not in tmp:
                        continue
#                    tmp[kr.id].kraken_clade += kr.reads_clade
#                if kr.id == 10239 and tmp[kr.id].reads_clade < kr.reads_clade:
#                    print ("{}\t{}\t{}".format(tmp[kr.id].reads_clade, kr.reads_clade,  arr[0]))
                    
#                    exit()
    start_node = 1
#    print(str(total_reads) + " total reads" )
#root in podoviridae
#    start_node = 10744
    for br in bracken_all:
        bracken_all[br].percentage = 100.00 * bracken_all[br].reads_clade / bracken_all[start_node].reads_clade
    for br in bracken_all:
        bracken_all[br].childs = {}
    for br in bracken_all:
        if br != 1:
            bracken_all[nodes[br].parent_id].childs[br] = bracken_all[br].reads_clade
    def print_childs(id):
#100.00  443943  0       D       10239     Viruses
        info = bracken_all[id]
        res = "{:.2f}\t{}\t{}\t{}\t{}\t{}".format(info.percentage, info.reads_clade, info.reads_self, info.level, info.id, info.name)
# + "\t" + str(info.reads_clade) + "\t" + str(info.reads_self) +  "\t" + info.level + "\t" + str(info.id) + "\t" + info.name
        if add_kraken:
            res += "\t{}\t{}\t{:.4f}".format(info.kraken_clade, info.kraken_self, info.kraken_clade/info.reads_clade)
#        print ("{:.2f}".format(info.percentage) + "\t" + str(info.reads_clade) + "\t" + str(info.reads_self) +  "\t" + info.level + "\t" + str(info.id) + "\t" + info.name)
        print (res)
        for child in sorted(info.childs, key=info.childs.get, reverse=True):
            print_childs(child)
    print_childs(start_node)

if __name__ == "__main__":
    merge_brackens(workdir)
#run_all_bracken(list, workdir)
#run_all_kraken(list, inputdir, workdir)

