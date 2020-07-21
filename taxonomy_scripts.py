#!/usr/bin/python3
import sys
import os
import random
from os import listdir
from os.path import isfile, join

from joblib import Parallel, delayed
import csv
import subprocess

#end of used taxids
start_taxid = 2748050
taxonomy_table = "/Bmo/dantipov/gut_pipeline/june_abund/table_1.tsv"
crassfamily_id = 1978007
poroviridae_id = 10744
unclass_podoviridae_id = 196895
input_dir = "/Bmo/dantipov/gut_pipeline/june_abund/kraken_genomes_saved/"
patched_dir = "/Bmo/dantipov/gut_pipeline/june_abund/kraken_genomes_splitted/"
kraken_build = "/home/dantipov/other_tools/kraken2/kraken/kraken2-build"
#kraken_db = "/Bmo/dantipov/gut_pipeline/kraken_viral_db/"
kraken_db = "/Bmo/dantipov/gut_pipeline/kraken_viral_db_clusters/"
#kraken_db = "/Bmo/dantipov/gut_pipeline/standard_db_updated/"
kraken_names = kraken_db + "/taxonomy/names.dmp"
kraken_nodes = kraken_db + "/taxonomy/nodes.dmp"
kraken_original_names = kraken_db + "/taxonomy/names_save.dmp"
kraken_original_nodes = kraken_db + "/taxonomy/nodes_save.dmp"
kraken_cleared_names = kraken_db + "/taxonomy/names_cleared.dmp"
kraken_cleared_nodes = kraken_db + "/taxonomy/nodes_cleared.dmp"
kraken_additional_names =  kraken_db + "/taxonomy/names_additional.dmp"
kraken_additional_nodes =  kraken_db + "/taxonomy/nodes_additional.dmp"
fake_subgroups ={"partial genome", "phi17_2_NC_021798", "phi13_2_NC_021803", "MethylophagaPAJP01000008"}
clusters_f = "/Bmo/dantipov/gut_pipeline/june_abund/clusters.cls"


all_names_file = ""
class reference_stats:
    def __init__(self, name, sample_id):
        self.name = name
        self.sample_id = sample_id
        self.group = ""
        self.subgroup = ""
        self.cluster = ""
    def __str__(self):
        return(self.name +" " + str(self.sample_id) + " " + self.group)

class node:
    def __init__(self, name, sample_id, parent_id, isleaf):
        self.name = name
        self.sample_id = sample_id
        self.parent_id = parent_id
        self.isleaf = isleaf
        self.depth = -1
        self.iscluster = False
#used to remove subtree
        self.bad = "unknown"
    def __str__(self):
#(f"{taxid}\t|\t1978007\t|\tspecies\t|\tCS\t|\t3\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t1\t|\t0\t|\t|\t")
#print(f"{taxid}\t|\t{line.split()[0]}\t|\t \t|\tscientific name\t|\t")

        return(nodes_str(self))

    def nodes_str(self):
        if  self.isleaf == False:
            if self.iscluster:
                rank = "genus"
#more precise bullshit?
            else:
                rank = "clade"
        else:
            rank = "species"
        return (f"{self.sample_id}\t|\t{self.parent_id}\t|\t{rank}\t|\tCS\t|\t3\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t1\t|\t0\t|\t\t|")

    def names_str(self):
        return(f"{self.sample_id}\t|\t{self.name}\t|\t \t|\tscientific name\t|")

def add_names (order, nodes):
    outf = open(kraken_additional_names, "w")
    for ref in order:
        outf.write(nodes[ref].names_str() + "\n")

def add_nodes (order, nodes):
    outf = open(kraken_additional_nodes, "w")
    for ref in order:
        outf.write(nodes[ref].nodes_str() + "\n")

def patch_fasta (order, nodes):
    os.system ("rm " + patched_dir + "/*.fasta")
    outf = open(join(patched_dir,   "concat.fasta"),"w")
    for ref in order:
        inf = join(input_dir, ref + ".fasta")
        if os.path.isfile(inf):
           # outf = open(join(patched_dir, ref + ".fasta"),"w")
            for line in open (inf, "r"):
                if line[0] ==">":
                    line = f">kraken:taxid|{nodes[ref].sample_id}|{nodes[ref].name}|" + line[1:]
                outf.write(line)
   


def create_taxonomy():
    refs = {}
    sample_id = start_taxid
    groups = {}
    subgroups = {}

    clusters = []   
    name_to_clust = {}
    check = {}
    clust_id = 0
    clusters_to_sub = {}
    for line in open(clusters_f, 'r'):
        arr = line.strip().split()
        if len(arr) == 0:
            break
        clusters.append(set(arr[1:]))
        clust_id +=1
        clust_name = clust_id 
        for name in arr[1:]:
            name_to_clust[name] = str(clust_name)
    subs_to_clust = {}
    with open (taxonomy_table, "r") as infile:#
        for line in infile:
            #print (line.split())
            arr= line.split("\t")
            if arr[2] == "translation":
                continue
            sample_id += 1
            name = arr[0].strip()
            group = arr[3].strip()
            subgroup = arr[4].strip()
            if subgroup in fake_subgroups:
                subgroup = ""
            subgroups[subgroup] = group 
            if group not in groups:
                groups[group] = set()            
            groups[group].add(subgroup)
#nt ID  len, nt translation     group   subgroup        representative? source           
            refs[sample_id] =  reference_stats(name,sample_id)
            refs[sample_id].group = group
            refs[sample_id].subgroup =  subgroup
            if name in name_to_clust and subgroup != "":
#                old_name = subgroup + "__" + name_to_clust[name]
                if subgroup not in subs_to_clust:
                    subs_to_clust[subgroup] = {}
                old_name = name_to_clust[name]
                if old_name in subs_to_clust[subgroup]:
                    new_name = subs_to_clust[subgroup][old_name]
                else:
                    new_name = subgroup+ "__" + str(len(subs_to_clust[subgroup]) + 1)
                    subs_to_clust[subgroup][old_name] = new_name
                clust = new_name
                refs[sample_id].cluster = clust
                clusters_to_sub[clust] = subgroup
    nodes = {}
    order =[]
    for group in groups:
        if group.strip() == "outgroup":
#            parent_id = poroviridae_id
            continue
        sample_id += 1
        parent_id = unclass_podoviridae_id
        nodes[group] = node(group, sample_id, parent_id , False)
#        for supgroups in groups[group]:
#            subgroups[subgroup] = sample_id
        order.append(group)
    for subgroup in subgroups:
        if subgroup == "":
            continue
        if subgroups[subgroup] == "outgroup":
            continue
        sample_id += 1
        nodes[subgroup] = node(subgroup, sample_id, nodes[subgroups[subgroup]].sample_id, False)
        order.append(subgroup)
    for cl in clusters_to_sub:
        sample_id += 1
        nodes[cl] = node(cl, sample_id, nodes[clusters_to_sub[cl]].sample_id, False)
        nodes[cl].iscluster = True
        order.append(cl)

#    for cl in clusters:
#        sample_id += 1
            
    leafs = []
    for ref_id in refs:
        ref = refs[ref_id]
        parent = -239
        print (ref)
        if ref.group == "outgroup":
            parent = poroviridae_id
        else:
            parent = nodes[ref.group].sample_id
            if ref.subgroup != '':
                parent = nodes[ref.subgroup].sample_id
            if ref.cluster != '':
                parent = nodes [ref.cluster].sample_id
        nodes[ref.name] = node (ref.name, ref.sample_id, parent, True)
        order.append(ref.name)

    add_names(order, nodes)
    add_nodes(order, nodes)  
    print("additions for taxonomy created..")  
    patch_fasta (order, nodes)
    print ("fasta patched..")                
#    for ref in refs:
#        print (refs[ref])


def add_to_lib():
# --add-to-library /Bmo/dantipov/gut_pipeline/june_abund/kraken_genomes_splitted/crAssphage.fasta --db  /Bmo/dantipov/gut_pipeline/june_abund/kraken_genomes_splitted/
    for ref in listdir(patched_dir):
        os.system(kraken_build + " --add-to-library " + join(patched_dir,ref) +" --db " + kraken_db)    


def read_nodes(names_file, nodes_file):
    id_to_names = {}
    depth = {}
    nodes = {}
    for line in open (names_file, "r"):
        elems = line.strip().split("\t|\t")
        if elems[-1].split("\t")[0] == "scientific name":
            sample_id = int(elems[0])
            id_to_names[sample_id] = elems[1]
    with open(nodes_file, 'r') as f:
        for line in f.readlines():
            elems = line.strip().split("\t|\t")
            sample_id = int(elems[0])
            name = id_to_names[sample_id]
            parent_id = int(elems[1])
            nodes[sample_id] = node (name, sample_id, parent_id, True)
            if sample_id == 1:
                nodes[sample_id].depth = 0
    def calc_depth(node):
        if node.depth >= 0:
            return node.depth
        else:
            return calc_depth(nodes[node.parent_id]) + 1

    for f in nodes:
        if nodes[f].depth == -1:
            nodes[f].depth = calc_depth(nodes[f]) 

    return nodes

def add_names_for_braken():
    nodes = read_nodes(kraken_names, kraken_nodes)
#    for f in nodes:
#        print (nodes[f].name + " " + str (nodes[f].depth))
    sp = " "
    report_file ="/home/dantipov/other_tools/kraken2/ERR688574.report"
    for line in open(report_file, "r"):     
        if (len(line.split("\t")[5].strip())) != 0:
            print (line.strip())
        else:
            sample_id = int(line.split("\t")[4])
            print (line.strip() + "\t" +  nodes[sample_id].depth* 2 * sp +nodes[sample_id].name)

def clear_extra_nodes():
    nodes = read_nodes (kraken_original_names, kraken_original_nodes)
    nodes[1].bad = "no"
    nodes[unclass_podoviridae_id].bad = "yes"
    def calc_depth(node):
        if node.depth >= 0:
            return node.depth
        else:
            return calc_depth(nodes[node.parent_id]) + 1

    def calc_bad(node):
        if node.bad == "unknown":
            calc_bad(nodes[node.parent_id])
            node.bad = nodes[node.parent_id].bad 

    for f in nodes:
        calc_bad(nodes[f]) 
    bad_ids = set()
    for node in nodes:
        if nodes[node].bad == "yes":
            bad_ids.add(node)
    nodes[unclass_podoviridae_id].name = "Extended crAssphage family"
    patched_namesf = open(kraken_cleared_names, 'w')
    patched_nodesf = open(kraken_cleared_nodes, 'w')
    root_patched = False
    
    for line in open(kraken_original_names, 'r'):
        id = int(line.split()[0])
        if not id in bad_ids:
            patched_namesf.write(line.strip() + "\n")
        elif id == unclass_podoviridae_id and not root_patched:
            patched_namesf.write(nodes[id].names_str() + "\n")
            root_patched = True
    for line in open(kraken_original_nodes, 'r'):
        id = int(line.split()[0])
        if id == unclass_podoviridae_id or not id in bad_ids:
            patched_nodesf.write(line.strip() + "\n")
    print("cleared internal podoviridae nodes..")            

def prepare_db():
    create_taxonomy()
    clear_extra_nodes()
    os.system("cat "+ kraken_cleared_names +" "+ kraken_additional_names + " > " + kraken_names)
    os.system("cat "+ kraken_cleared_nodes +" "+ kraken_additional_nodes + " > " + kraken_nodes) 
#    os.system("cat "+ kraken_cleared_names + " > " + kraken_names)
#    os.system("cat "+ kraken_cleared_nodes + "  > " + kraken_nodes)
    os.system("rm " + kraken_db + "/*.k2d")
    os.system("rm " + kraken_db + "/seqid2taxid.map")
    os.system("rm "+ kraken_db + "/library/added/*")
    print ("cleared old files...")
    add_to_lib()
    print ("added new sequences to lib...")
    os.system(kraken_build + " --threads 20 --build --db " + kraken_db)

#create_taxonomy()
prepare_db()
#add_to_lib()
#add_names_for_braken()
