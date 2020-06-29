#!/usr/bin/python3

from Bio import Entrez
import sys

sra_list=[]
for line in open(sys.argv[1], "r"):
     sra_list.append(line.strip().split()[0])

Entrez.email = "mike.rayko@gmail.com"


retmax = 10
for i in range (0, len(sra_list),retmax):
    id_list = ",".join(sra_list[i: min(i+retmax, len(sra_list))])

    try:
        handle = Entrez.efetch(db="sra", id = id_list, rettype="runinfo", retmode="text")
    except:
        print (f'Failed in {str(i)} - {min(i+retmax, len(sra_list))}')
        continue

    meta = handle.readlines()
    for i in range (1, len(meta)-1):
        meta_dict = dict(zip(meta[0].split(","), meta[i].split(",")))
        if meta_dict["Platform"] != "ILLUMINA":
            print(meta_dict["Run"] + "\t "+ "Non_illumina")
        else:
            if meta_dict["LibraryLayout"] == "PAIRED":
                print(meta_dict["Run"] + "\t "+ str(int(meta_dict["avgLength"])/2)) 
            else:
                print(meta_dict["Run"] + "\t "+ str(int(meta_dict["avgLength"])))


    handle.close()

