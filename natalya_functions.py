import sys 
import os 
import string 
import re 
import subprocess 
import datetime 
import fastaparser 
from os.path import join, basename
from genericpath import isdir, exists
from joblib import Parallel, delayed
from operator import itemgetter

data_dir = "/Nancy/mrayko/Libs/virsorter-data-hallmarks/"
mash_bin = "/home/dantipov/other_tools/mash/mash"
seqtk_bin = "/home/dantipov/other_tools/seqtk/seqtk"



'''
Takes dir with spades result as an input, 
reads saffols.fasta, for each:
opens dir with k-mer assembly (from 129 down to 50), select largest.
Checks if it overlaps in start and end.
If yes, prints out contig name with "is circular"
'''

def extract_circular_from_file(file, indir, outdir):
    out_file = join(outdir, os.path.splitext(file)[0] + ".circular.fasta")
    contigs = fastaparser.read_fasta(join(sys.argv[1], file))
    circulars = []
    count = []
    for contig in contigs:
        arr = contig[0].strip(';').split('_')
#      if float(arr[3]) > 500:
        if len(contig[1]) < 500: continue
        for kval in range (200, 50, -1):
#            kval = 55
            if kval >= len(contig[1]) or len(contig[1]) < 500:
                continue
            start = contig[1][:kval]
            end = contig[1][-kval:]
       
            if start == end:
#               print (">" + contig[0][1:])
#               print (contig[1])
#                print (" k equal " + str(kval))
                print (contig[0] + " is circular " + str(kval))
#                contig[0] = contig[0] + " k: " + str(kval)
                circulars.append(contig)
                break   
    fastaparser.write_fasta_to_file(out_file, circulars)

def run_virsorter_one(infile, outcommondir):
    outdir = join(outcommondir, os.path.splitext(os.path.basename(infile))[0])

    #wrapper_phage_contigs_sorter_iPlant.pl -f ///Bmo/dantipov/gut_pipeline/circulars/GCA_900276575.1_SRR761713_genomic.fna.circular.fasta  --db 2 --wdir /Bmo/dantipov/tmp/virs/ --virome --ncpu 10 --data-dir   /Nancy/mrayko/Libs/virsorter-data-hallmarks/ 
    run_string = "wrapper_phage_contigs_sorter_iPlant.pl -f  " +infile + " --db 2 --wdir " + outdir + " --virome --ncpu 50 --data-dir " + data_dir
    print (run_string)
    os.system(run_string)

def run_virsorter_all(indir, outcommondir):
#    os.system("activate virsorter")
#    file = "/Bmo/dantipov/gut_pipeline/circulars/GCA_002924505.1_ASM292450v1_genomic.fna.circular.fasta"
#    file = "/Bmo/dantipov/gut_pipeline/circulars/GCA_900284195.1_SRR3131891_genomic.fna.circular.fasta"
#    file = "/Bmo/dantipov/gut_pipeline/circulars/GCA_900270585.1_SRR1214756_genomic.fna.circular.fasta"
#     run_virsorter_one(file, outcommondir)
#    Parallel(n_jobs=7)(delayed(run_virsorter_one)(join(indir, file), outcommondir)
     for file in os.listdir(indir):
        run_virsorter_one(join(indir,file), outcommondir)


#GCA*.*.fna
def glue_and_rename (indir, outfile):
    for file in os.listdir(indir):
        arr = file.split('.')
        if len(arr) < 4: 
            continue
        contigs = fastaparser.read_fasta(join(indir, file))
        for contig in contigs:
            new_name = contig[0] + " " + arr[0] + "." +arr[1]
            print new_name
            fastaparser.write_fasta_to_file(outfile, zip([new_name], [contig[1]]))

def split_and_rename(infile, outdir):
    contigs = fastaparser.read_fasta(infile)
    for contig in contigs:
        filename = contig[0].split()[0][1:]
        filename = join(outdir, filename + ".fasta")
        print filename
        fastaparser.write_fasta_to_file(filename, [contig])


#directory with virsorter output, all contigs
def extract_prophages(indir, contigs_file):
#    /Bmo/dantipov/gut_pipeline/concatenated_circulars/all_circulars_glued/Predicted_viral_sequences/VIRSorter_prophages_cat-
    os.system("grep \'>\' " + indir+ "VIRSorter_prophages_cat*fasta > tmp1.grep")
    ids = open("tmp2.grep", 'w')
    for line in open("tmp1.grep", 'r'):
        start_line = line.split('>')[1][10:]
        arr = start_line.split('_')
        seq_id = arr[0] + "." + arr[1]
#        process = subprocess.Popen(['grep', seq_id, contigs_file], stdout=subprocess.PIPE)
#        stdout = process.communicate()[0]        
#        ids.write(stdout)
        ids.write(seq_id + "\n")
#    os.system(seqtk_bin + " subseq " + contigs_file + " tmp2.grep " )
    os.system(seqtk_bin + " subseq " + contigs_file + " tmp2.grep  > " + join(indir, "Prophages_extracted.fasta"))

def process_extracted(indir):
    all_viruses = join(indir, "all_viruses.fasta")
    os.system ("rm " + all_viruses)
    os.system("cat " + indir + "/VIRSorter_cat*fasta " + indir + "/Prophages_extracted.fasta > " + all_viruses)

def extract_not_listed(infasta, list):
    listed = set()
    for line in open (list, 'r'):
        listed.add(">" + line.split()[0])
    print len(listed)
    contigs = fastaparser.read_fasta(infasta)
    print len (contigs)
    outcontigs = []
    for contig in contigs:
        if not contig[0].split()[0] in listed:
#            print contig[0]
            outcontigs.append(contig)
#        else:
#            listed.remove(contig[0])
#    for c in listed:
#        print c
    print len (outcontigs)
    outfasta = infasta[:-6]+".unknown.fasta"
    os.system("rm "+ outfasta)
    fastaparser.write_fasta_to_file(outfasta, outcontigs)
       

def run_mash(all_viruses):
    os.system(mash_bin +" sketch -i " + all_viruses)
    all_sketches = all_viruses + ".msh"
    os.system(mash_bin +" dist " + all_sketches + " " + all_sketches + " -p 50 > " +all_viruses+".dist.table")
    

#need fastas and sketch
def found_most_similar(work_dir):
    contigs_info = []
    for file in os.listdir(work_dir):
        arr = file.split('.')
        if arr[-1] == "fasta":
            contigs = fastaparser.read_fasta(join(work_dir, file))
            contigs_info.append([file, len(contigs[0][1])])

    all_sorted = sorted(contigs_info, key = itemgetter(1))
    max_ind = len(all_sorted)
    low_ind = 0
    high_ind = 0
    similar_list = []
    used = []
    for i in range(0, max_ind):
        used.append(False)
        cur_len = all_sorted[i][1]
        first_mash = join(work_dir, all_sorted[i][0] + ".msh")
        while all_sorted[low_ind][1] < cur_len * 0.8 and low_ind < max_ind -1:
            low_ind +=1
        while all_sorted[high_ind][1] < cur_len * 1.2 and high_ind < max_ind :
            high_ind +=1
        if i % 10 == 0:
            print "processing... " + str(i) + " range: " + str(low_ind) + "-" + str(high_ind) 
        sim = []
        for j in range(low_ind, high_ind):
            second_mash =  join(work_dir, all_sorted[j][0] + ".msh")
            process = subprocess.Popen([mash_bin, 'dist', first_mash, second_mash], stdout=subprocess.PIPE)
            stdout = process.communicate()[0]
            arr = stdout.split()
            dist = float(arr[2])
            if dist < 0.2:
                sim.append(j)
        similar_list.append([i, len(sim), sim])   
        if i % 10 == 0:
            print(len(sim))
    most_similar = sorted(similar_list, key = itemgetter(1), reverse = True)
    for k in most_similar:
        print k
    for contigs in most_similar:
        print all_sorted[contigs[0]][0] + " " + str(contigs[1]) +" " + str(used(contigs[0]))
        for j in contigs[2]:
            used[j] = True

def get_short_name(fa_name):
    return fa_name[1:].split()[0]

def clean_table(table_file):
    for line in open(table_file, 'r'):
        arr = line.split()
        if float(arr[2]) <= 0.1:
            print line.strip()

def get_name(str):
    return str[2:-2]

def find_closest(old_to_new_file, max_dist):
    old_to_new = {}
    for line in open (old_to_new_file, 'r'):
        arr = line.split()
        old_name = arr[0]
        new_name = arr[1]
        dist = float(arr[2])
        if old_name not in old_to_new and dist <= float(max_dist):
            old_to_new[old_name] = new_name
    for old in old_to_new.keys():
        print old + " " + old_to_new[old]

def extract_interesting(old_to_new_file, max_dist):
    old_to_new = {}
    for line in open (old_to_new_file, 'r'):
        arr = line.split()
        old_name = arr[0]
        new_name = arr[1]
        dist = float(arr[2])
        if  dist <= float(max_dist):
            if old_name not in old_to_new:
                old_to_new[old_name] = []
            old_to_new[old_name].append(new_name)
    for old in old_to_new.keys():
        old_f = open(old+".list", "w")
        for new in old_to_new[old]:
            old_f.write(new.strip()+"\n")


def parse_mash(contig_file, table):
    contigs = fastaparser.read_fasta(contig_file)
    similar_lists = {}
    for contig in contigs:
        similar_lists[get_short_name(contig[0])]=[]
    for line in open(table, 'r'):
        arr = line.split()
        dist = float(arr[2])
        if dist < 0.1:
            similar_lists[arr[0]].append(arr[1])
    print "processed input"
    to_sort = []
    for l in similar_lists:
        to_sort.append([l, len(similar_lists[l])])
    sorted_similar = sorted(to_sort, key = itemgetter(1), reverse = True)
    outcontigs = []
    used = set()
    for contig_info in sorted_similar:
        if contig_info[0] not in used:
            for similar in similar_lists[contig_info[0]]:
                used.add(similar)
            if contig_info[1] > 10:
                print contig_info
                print similar_lists[contig_info[0]]
#far from optimal but whynot                
                for contig in contigs:
                    if get_short_name(contig[0]) == contig_info[0]:
                        outcontigs.append(contig)
                        break
    result_f = join(os.path.dirname(contig_file), "interesting.fasta")
    os.system("rm " + result_f)
    fastaparser.write_fasta_to_file(result_f, outcontigs)


def check_mash (best_list, table_file):
    centers = {}
    centers_list = []
    for line in open(best_list, "r"):
        x = eval(line)
        centers[x[0]] = x[1]
        centers_list.append(x[0])
    used = set()
    next_pairs = []
    for line in open (table_file, "r"):
        arr = line.split()
        first = arr[0]
        second = arr[1]
        dist = float(arr[2])
        if dist <= 0.1:
            if len (next_pairs) == 0 or next_pairs[-1][0] != second:
                next_pairs.append([second, []])
            next_pairs[-1][1].append(first)
#    for line in open (table_file, "r"):
#        arr = line.split()
#        first = arr[0]
#        second = arr[1]
#        dist = float(arr[2])
#        if dist <= 0.1:
#            if not (first in next.keys()):
#                next[first] = []
#            next[first].append(second)
    next = dict(next_pairs)
    for l in centers_list:
        neighbours = len(next[l])
        count = 0
        members = []
        for x in next[l]:
            if x not in used:
                count += 1
                used.add(x)
                members.append(x)
        print "cluster info: " + l + "\t" + str(neighbours) + "\t" + str(count)
        for m in members:
            print m
    to_sort = []
    for l in next.keys():
        if l not in used:
            to_sort.append([l, len(next[l])])
    sorted_similar = sorted(to_sort, key = itemgetter(1), reverse = True)
    outcontigs = []
#    used = set()
    for l in sorted_similar:
        if l[0] not in used:
            count = 0
            members = []
            for similar in next[l[0]]:
                if similar not in used:
                    used.add(similar)
                    count += 1
                    members.append(similar)
            print "small_cluster info: " + str(l[0]) + "\t" + str(l[1]) + "\t" + str(count)
            for m in members:
                print m

def remove_extra(infile):
    clusts = []
    count = 0
    name = ""
    n_size = 0
    cur_list = []
    for line in open (infile, 'r'):
        arr = line.split()        
        if len(arr) < 5:
            if count == 0:
                print line.strip()
            else:
                cur_list.append(line.strip())
                count -= 1
        else:
            if len(clusts) >= 1 and name == clusts[-1][0]:
                clusts[-1][2] = cur_list 
                clusts[-1][1] = n_size
            else:
                clusts.append([name, n_size, cur_list])
            name = arr[2]
            n_size = int(arr[3])
            count = int(arr[4])
            cur_list = []
    for x in clusts:
        intro = "cluster info: "
        if int(x[1]) <= 10:
            intro = "small_" + intro
        print intro + x[0] + "\t" + str(x[1]) + "\t" + str(len(x[2]))
        for next in x[2]:
            print next
        

def get_clust_size(clust_size_file, interesting_list):
    clust_sizes = {}
    for line in open(clust_size_file, "r"):
        arr = line.split()
        if len(arr) > 0:
            name = arr[2]
            size = arr[4]
            clust_sizes[name] = size
    for line in open (interesting_list, "r"):
        name = line.strip()
        print clust_sizes[name]

def run_quast_all (assemblies_dir, ref, out_dir):
    for ass in os.listdir(assemblies_dir):
        if len(ass.split('.')) > 2 and ass.split('.')[-2] == "fna":
            string = "quast.py -t 30 --min-identity 90 --fast " + join(assemblies_dir, ass) + " -R " + ref + " -o " + join (out_dir, basename(ass[:-7]))
            os.system(string)

def parse_quast_grep(infile):
    counts = []

    for i in range(0, 10): 
        counts.append(0)
    for line in open(infile, "r"):
        arr = line.split()
        if len(arr) >=3:
            frac = int(float(arr[3])/10)
            counts[frac] += 1
    print counts
#blastn  -query /Bmo/dantipov/gut_pipeline/concatenated_circulars/all_circulars_glued/Predicted_viral_sequences/all_viruses.fasta  -db /Bmo/ncbi_nt_database/nt  -evalue 0.0001 -outfmt 5 -out /Bmo/dantipov/gut_pipeline/concatenated_circulars/all_circulars_glued/Predicted_viral_sequences/all_viruses.xml -num_threads 40 -num_alignments 50
#python parse_blast_xml.py -i /Bmo/dantipov/gut_pipeline/concatenated_circulars/all_circulars_glued/Predicted_viral_sequences/all_viruses.xml -o . > /Bmo/dantipov/gut_pipeline/concatenated_circulars/all_circulars_glued/Predicted_viral_sequences/all_viruses.xml.parsed
#grep "KNOWN" /Bmo/dantipov/gut_pipeline/concatenated_circulars/all_circulars_glued/Predicted_viral_sequences/all_viruses.xml.parsed > /Bmo/dantipov/gut_pipeline/concatenated_circulars/all_circulars_glued/Predicted_viral_sequences/all_viruses.xml.parsed.known
#/home/dantipov/other_tools/seqtk/seqtk subseq /Bmo/dantipov/gut_pipeline/concatenated_circulars/all_circulars_glued/Predicted_viral_sequences/all_viruses.fasta /Bmo/dantipov/gut_pipeline/concatenated_circulars/all_circulars_glued/Predicted_viral_sequences/all_viruses.xml.parsed.known > /Bmo/dantipov/gut_pipeline/concatenated_circulars/all_circulars_glued/Predicted_viral_sequences/all_viruses.xml.parsed.known.fasta


if __name__ == "__main__":
#    run_virsorter_all(sys.argv[1], sys.argv[2])
#    dir = sys.argv[1]
#    for file in os.listdir(dir):
#        extract_circular_from_file(file, dir, sys.argv[2])
#    glue_and_rename(sys.argv[1], sys.argv[2])
#    run_virsorter_one(sys.argv[1], sys.argv[2])
#    extract_prophages(sys.argv[1], sys.argv[2])
#    process_extracted(sys.argv[1])
#    extract_not_listed(sys.argv[1], sys.argv[2])   
#    run_mash(sys.argv[1])
#    parse_mash(sys.argv[1], sys.argv[2])
#    clean_table(sys.argv[1])
#    find_closest(sys.argv[1], sys.argv[2])
#    extract_interesting(sys.argv[1], sys.argv[2])
    check_mash(sys.argv[1], sys.argv[2])
#    get_clust_size(sys.argv[1], sys.argv[2])
#    run_quast_all(sys.argv[1],sys.argv[2], sys.argv[3])
#    parse_quast_grep(sys.argv[1])
#    remove_extra(sys.argv[1])
