#!/usr/bin/python
## input is like as follows:
# NODE_11062_length_3179_cov_4.856914__Array72__9 742738  Flavonifractor plautii 1_3_50AFAA       species SRR6000869      left    S001388

## Abstract : each contig like 'NODE_11062_length_3179_cov_4.856914' should have only one match of host 
from collections import Counter
import os
from sys import argv
#usage: python crispr_PC2host.py crispr_array_taxa_interaction.txt
infile = argv[1]
contig2host = {}   ## store contig-> host_taxid
taxid2taxa = {}  ## taxaid -> taxaname
taxid2level = {}  ## taxaid -> taxalevel
lines = open(infile,'r')
next(lines)

for line in lines:
        infos = line.split('\t')
        array = infos[0]  ## NODE_11062_length_3179_cov_4.856914__Array72__7
        contig = array.split('__')[0] ## NODE_11062_length_3179_cov_4.856914
        taxid = infos[1].strip()
        contig2host.setdefault(contig,[]).append(taxid)
        taxaname,taxalevel = infos[2].strip(),infos[3].strip()
        taxid2taxa[taxid] = taxaname
        taxid2level[taxid] = taxalevel

for key in contig2host.keys():
        taxids = contig2host[key]
        taxa = Counter(taxids).most_common(1)[0][0]
        contig2host[key] = taxa


clusters='/home/jfh/data/Virome/Gajewski/03_cdhit/cluster_results.txt'
seq2clu = {}
cluster_result = open(clusters,'r')
next(cluster_result)
for result in cluster_result:
        infos = result.split('\t')
        if infos[2].strip()=='1':
                seq = infos[3].strip()
        elif infos[2].strip()=='0':
                seq = infos[4].strip()
        clu = 'PC'+infos[1].strip()
        seq2clu[seq] = clu


outFile = open('crispr_all_repeat_interaction.txt','w')
#outFile.write("PC"+'\t'+"host_taxa_id"+'\t'+"host_taxa_level"+'\t'+"host_taxa_name"+'\t'+"sample_source\n")
lines = open(infile,'r')
next(lines)
for line in lines:
        infos = line.split('\t')
        array = infos[0]  ## NODE_11062_length_3179_cov_4.856914__Array72__7
        contig = array.split('__')[0] ## NODE_11062_length_3179_cov_4.856914
        host_taxa_id = contig2host[contig]
        host_taxa_level = taxid2level[host_taxa_id]
        host_taxa_name = taxid2taxa[host_taxa_id]
        seq = infos[6].strip()
        PC = seq2clu[seq]
        outFile.write(PC+'\t'+host_taxa_id+'\t'+host_taxa_level+'\t'+host_taxa_name+'\t'+infos[4].strip()+'\n')
outFile.close()
cmd = "cat crispr_all_repeat_interaction.txt|sort|uniq > crispr_all_interaction.txt;rm crispr_all_repeat_interaction.txt"
os.popen(cmd,'r')
with open('crispr_all_interaction.txt','r+') as f:
    content = f.read()        
    f.seek(0, 0)
    f.write("PC"+'\t'+"host_taxa_id"+'\t'+"host_taxa_level"+'\t'+"host_taxa_name"+'\t'+"sample_source\n"+content)
f.close()