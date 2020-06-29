#!/usr/bin/python


import sys
import os
import random
import gzip
from os import listdir
from os.path import isfile, join
infile = "/Bmo/dantipov/gut_pipeline/abundancy_check/abundancies_extended.tsv" 
namesfile = "/Bmo/dantipov/metagenomes_Koonin/june_abund/names.txt"
first_line = "contig\tfraction_of_reads\tlength\tnormalized_fraction_of_reads\treads_available\tsample_name"

lines = {}
for line in open (infile, 'r'):
    arr = line.split()
    lines[arr[0].split('.')[0]] = line.strip()
#    print arr[0]
print first_line
for line in open (namesfile, 'r'):
    if line[0] == ">":
        name = line.split()[0][1:]
#        print name
        if name in lines:
#            print name        
            print lines[name]


#get_all_abundancies(sys.argv[1], sys.argv[2])
#find_sample_assembly(sys.argv[1])    
#def get_containing_filename(pattern):
    

