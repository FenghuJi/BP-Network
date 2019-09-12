#!/usr/bin/python
## Our goal is to create a reader friendly interaction table for porphage like that :
## phage_cluster_id || host_taxa_id ||host_taxa_level||host_taxa_name||sample_source
## Also, a interaction matrix for each sample

import os
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-p","--prophage",dest="prophage",help="set the prophage interaction text file",metavar="prophage interaction")  ##/home/jfh/data/Virome/Gajewski/07_interaction/gut_prophage_interaction.txt
parser.add_option("-o","--output",dest="outfile",help="set the output filename",metavar="output")  ## /home/jfh/data/Virome/Gajewski/07_interaction/prophage_result.txt
parser.add_option("-m","--namemap",dest="name2seq",help="set the name2seq mapping text file",,metavar="name mapping file") # /home/jfh/data/Virome/Gajewski/03_cdhit/prophage_virus.names.new2old.txt
(options,args) = parser.parse_args()
prophages = options.prophage
output = options.outfile
nameseq = options.name2seq


##Some variable to modify for different project
clusters='/home/jfh/data/Virome/Gajewski/03_cdhit/cluster_results.txt'


## Create dictionary of key(contig name) -> value (cluser id )
## contig name -> seq_id --------- seq_id -> cluster ------------contig name -> clutser
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

## contig name to new seq name
contig2seq = {}
contigmaps = open(nameseq,'r')
for map in contigmaps:
        infos = map.split('\t')
        contig = infos[1]
        seq = infos[0]
        contig2seq[contig] = seq


### generate interaction file 
outFile = open(output,'w')
outFile.write("PC"+'\t'+"host_taxa_id"+'\t'+"host_taxa_level"+'\t'+"host_taxa_name"+'\t'+"sample_source\n")
interacts = open(prophages,'r')
next(interacts)
for interact in interacts:
        infos = interact.split('\t')
        contig = infos[0]
        try:
                PC = seq2clu[contig2seq[contig]]
        except:
                print("Can't find matched PC ID!")
        outFile.write(PC+'\t'+infos[1]+'\t'+infos[4]+'\t'+infos[5].strip()+'\t'+infos[2]+'\n')
outFile.close()
