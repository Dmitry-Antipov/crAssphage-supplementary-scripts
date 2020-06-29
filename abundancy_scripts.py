#!/usr/bin/python

#add some  contamination for meta-genomic simulation

import sys
import os
import random
import gzip
from os import listdir
from os.path import isfile, join
infile = sys.argv[1]

def count_kmer_cov(line):
#    print line
    arr = line.strip().split('_')
    cov = float(arr[-1].split(',')[0])
    len = float (arr[3])
    total = cov* len   
    return total

def count_total_kmer_cov(infile):
    total = 0 
    for line in gzip.open(infile, 'r'):
        if line[0] == ">" and line.find("NODE") != -1:
            total += count_kmer_cov(line)
    print total
    return total

def count_known_kmer_cov(infile, pattern):
    for line in gzip.open(infile, 'r'):
        if line[0] == ">" and line.find("NODE") != -1:
            if line.find(pattern) >=0:
                print "found"
                return count_kmer_cov(line)
    return -1

def get_relative_abund(infile, pattern):
    total =  count_total_kmer_cov(infile)
    if total == 0:
        return -1
    else:
        return count_known_kmer_cov(infile, pattern)/total


#print (get_relative_abund(sys.argv[1], sys.argv[2]))

def find_sample_assembly(pattern):
    circular_prefix = "/Bmo/dantipov/gut_pipeline/circulars/"
    common_prefix = "/Bmo/data/projects/human_gut_meta_NCBI/"

    grep_line = "grep \"" + pattern + "\" " + circular_prefix + "*.fasta > tmp.txt"
    print grep_line
    os.system(grep_line)
    line = open ("tmp.txt",'r').readline()
    full_name = line.split(":")[0]
    print full_name
    file_circ = os.path.basename(full_name)
    file_origin = file_circ[:-15] + ".gz"
    path_origin = join(common_prefix, file_origin)
    print path_origin
#    print (get_relative_abund(path_origin, pattern))
    return path_origin

def get_all_abundancies(infile, outfile):
    out_f = open(outfile,'w')    
    for line in open(infile, 'r'):
        print line
        pattern = line.split()[0]
        sample = find_sample_assembly(pattern)
        abund = get_relative_abund(sample, pattern)
        out_f.write(pattern + "\t" + str(abund) + "\n")

def get_length_from_pattern(pattern):
    circular_prefix = "/Bmo/dantipov/gut_pipeline/circulars/"
    grep_line = "grep " + pattern + " " + circular_prefix + "*.fasta > tmp.txt"
    print grep_line
    os.system(grep_line)
    line = open ("tmp.txt",'r').readline()
    start = line.find("NODE_")
    if start != -1:
        arr = line[start:].split('_')
        len = int (arr[3])
        return len
    else:
        return -1

def add_all_len(infile, outfile):
    outf = open(outfile, "w")
    for line in open(infile, 'r'):
        print line
        pattern = line.split()[0]
        unnormalized = float(line.split()[1])
        len = get_length_from_pattern(pattern)
        normalized = -1
        if len != -1 and unnormalized != -1:
            normalized = unnormalized / len
        name = find_sample_assembly(pattern)
        have_reads = False
        if name.find("SRR") != -1 or name.find("ERR") != -1:
            have_reads = True
        outf.write(line.strip() + "\t" + str(len) + "\t" + str(normalized) + "\t" + str(have_reads) + "\t" + name + "\n")
#add_all_len(sys.argv[1], sys.argv[2])
get_all_abundancies(sys.argv[1], sys.argv[2])
#find_sample_assembly(sys.argv[1])    
#def get_containing_filename(pattern):
    

