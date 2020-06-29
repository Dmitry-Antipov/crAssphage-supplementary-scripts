from Bio import Entrez
import sys

sra_list=[]
for line in open(sys.argv[1], "r"):
     sra_list.append(line.strip())

Entrez.email = "mike.rayko@gmail.com"



for i in sra_list:
    handle = Entrez.efetch(db="sra", id = i, rettype="runinfo", retmode="text") # or esearch, efetch, ...
    meta = handle.readlines()
    meta_dict = dict(zip(meta[0].split(","), meta[1].split(",")))
    if meta_dict["Platform"] != "ILLUMINA":
        print(meta_dict["Run"] + "\t "+ "Non_illumina")
    else:
        if meta_dict["LibraryLayout"] == "PAIRED":
            print(meta_dict["Run"] + "\t "+ str(int(meta_dict["avgLength"])/2)) 
        else:
            print(meta_dict["Run"] + "\t "+ str(int(meta_dict["avgLength"])))


    handle.close()

