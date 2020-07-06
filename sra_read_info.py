#!/usr/bin/python3
from Bio import Entrez
import sys

sra_list=[]
for line in open(sys.argv[1], "r"):
     sra_list.append(line.strip().split()[0])



def get_metadata(sra_list):
    Entrez.email = "mike.rayko@gmail.com"
    retmax = 200

    out_dict = {}
    for i in range (0, len(sra_list),retmax):
        id_list = ",".join(sra_list[i: min(i+retmax, len(sra_list))])
    
        try:
            handle = Entrez.efetch(db="sra", id = id_list, rettype="runinfo", retmode="text")
        except:
            print (f'Failed in {str(i)} - {min(i+retmax, len(sra_list))}')
            continue
        meta = []
        for line in handle.readlines():
            if len(line.strip().split(","))>1:
                if '"' in line:
                    heid = line.split('"')
                    meta.append((heid[0] + heid[1].replace(","," ") + heid[2]).strip().split(","))
                else:
                    meta.append(line.strip().split(","))

        for k in range (1, len(meta)):
            if len(meta[0]) != len(meta[k]) or meta[0][0] == meta[k][0]:
                continue
            else:
                meta_dict = dict(zip(meta[0], meta[k]))
                out_dict[meta[k][0]]=meta_dict
        handle.close()

    for i in sra_list:
        if i not in out_dict:
            meta_dict = {x: "NA" for x in meta_dict}
            out_dict[i]=meta_dict


    return (out_dict)
    

meta = get_metadata(sra_list)

for i in sra_list:
    if meta[i]["Platform"] != "ILLUMINA":
        print(i + "\t "+ "Non_illumina")
    else:
        if meta[i]["LibraryLayout"] == "PAIRED":
            print(f'{i} \t {str(int(meta[i]["avgLength"])/2)} \t {meta[i]["LibraryStrategy"]}') 
        else:
            print(f'{i} \t {str(int(meta[i]["avgLength"]))} \t {meta[i]["LibraryStrategy"]}')
            