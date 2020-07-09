!/usr/bin/python3

import sys
import csv
import os
from Bio import SeqIO

infile = "/Bmo/dantipov/gut_pipeline/june_abund/table_1.tsv"
genomes = "/Bmo/dantipov/gut_pipeline/june_abund/all_genomes.fa"


with open (infile, "r") as csv_input:
    parsed = list(csv.reader(csv_input, delimiter="\t"))
    
phages_dict = {}
types = set()
for i in parsed[1:]: 
    phages_dict[i[0]] = i[4]
    types.add(i[4])

d = dict((element,[]) for element in types) 

for record in SeqIO.parse(genomes, "fasta"):
    d[phages_dict[record.id]].append(record)


unique_kmers = {}
for i in d:
    if len(d[i])>2: 
        SeqIO.write(d[i], i+".fasta", "fasta")
        os.system("jellyfish count -m 100 -s 100M  -t 2 -C "+i+".fasta -o "+i+".jf")
        os.system("jellyfish stats "+i+".jf -o "+i+".stats")
        with open(i+".stats") as f:
            kmers = f.readlines()
            kmers = [i.strip().split(" ") for i in kmers] 
            print(kmers)
            unique_kmers[i] = [kmers[0][-1],kmers[1][-1],kmers[2][-1],str(len(d[i]))]
            total_kmers = 0
            for line in kmers:
                total_kmers+=int(line[1])
            unique_kmers[i] = [int(kmers[0][1]), total_kmers] 
            
print("Name Unique Distinct Total Sequences")
for i in unique_kmers:
    print(i, " ".join(unique_kmers[i]))
