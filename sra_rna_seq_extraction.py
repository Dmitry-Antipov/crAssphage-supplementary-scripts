
from Bio import Entrez
import sys
from xml.etree import cElementTree as ElementTree

Entrez.email = "mike.rayko@gmail.com"


sra_list=[]
for line in open(sys.argv[1], "r"):
    if len(line)>0:
      sra_list.append(line.strip().split()[0])


print(sra_list)

id_list = []
id_to_sra = {}
for i in sra_list:
    try: 
        handle = Entrez.esearch(db="sra", term = i)
        record = Entrez.read(handle)
        for id in record["IdList"]:
            id_list.append(id)
            id_to_sra[id] = i
        handle.close()
    except:
        print("No id in SRA database")
        continue

handle = Entrez.efetch(db="sra", id = id_list, retmode = "xml") # or esearch, efetch, ...

records = ElementTree.parse(handle)
root = records.getroot()

strategy = []

for library in root.iter('LIBRARY_STRATEGY'):
    strategy.append(library.text)


print(sra_list)
print(id_list)
print(strategy)
for i in range(0,len(id_list)):
    print(id_to_sra[id_list[i]], id_list[i], strategy[i])

handle.close()


