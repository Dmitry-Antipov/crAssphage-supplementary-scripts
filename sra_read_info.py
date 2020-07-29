#!/usr/bin/python3
from Bio import Entrez
import sys
import csv
def get_sra_info(infile):
    sra_list=[]
    for line in open(infile, "r"):
         sra_list.append(line.strip().split()[0])
    sra_set = set(sra_list)
    res = {}
    Entrez.email = "d.antipov@spbu.ru"


    retmax = 200
    for i in range (0, len(sra_list),retmax):
        id_list = ",".join(sra_list[i: min(i+retmax, len(sra_list))])

        try:
            handle = Entrez.efetch(db="sra", id = id_list, rettype="runinfo", retmode="text")
        except:
            print (f'Failed in {str(i)} - {min(i+retmax, len(sra_list))}')
            exit(1)

        meta = handle.readlines()
        print (meta)
        parsed = list(csv.reader(meta))
        names = parsed[0]
        for i in range(1, len(parsed)) :
#            names = list(csv.reader(StringIO(meta[0])))[0]
            values = parsed[i]
            if len(names) != len (values):
                continue
            else:
                meta_dict = dict(zip(names, values))
                if meta_dict["Run"] not in sra_set:
                    continue
                
                if meta_dict["LibraryLayout"] == "PAIRED":
                    avg_len = str(int(meta_dict["avgLength"])/2) 
                else:
                    avg_len = str(int(meta_dict["avgLength"]))
#ID strategy platform avrg_read_length number_of_reads
                res[meta_dict["Run"]] = (meta_dict["Run"] + " "+ meta_dict["LibraryStrategy"] + " "+ meta_dict["Platform"] + " " + avg_len +" "+ meta_dict["spots"])
        handle.close()
    for sra in sra_list:
        print(res[sra])

if __name__ == "__main__":
    get_sra_info(sys.argv[1])

