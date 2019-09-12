#!/usr/bin/python
import os
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-t","--taxa_bacteria",dest="taxa",help="set the bacteria taxa file",metavar="taxa_bacteria")
parser.add_option("-w","--array_path",dest="array_path",help="set the virus blast path",metavar="virus blast path")
parser.add_option("-l","--filelist",dest="filelist",help="set the filelist of samples",metavar="filelist")
(options,args) = parser.parse_args()
taxa_bacteria = options.taxa
array_path = options.array_path
filelist = options.filelist

files = open(filelist,'r')
taxa_results = open(taxa_bacteria,'r')

## Create dictionary  array-> [host taxid,host taxa_name,level,sample id,position】
taxa_map = {}
next(taxa_results)  
for result in taxa_results:
        infos = result.split('\t')
        arrays = infos[0].split('__')
        array = arrays[0]+'__'+arrays[1]
        taxid = infos[1]
        level,tax_name = infos[4],infos[5].strip()
        sample_name,position = sample_position = infos[3].split('/')[8].split('_')
        taxa_map[array]=[taxid,tax_name,level,sample_name,position]


outfile = open('crispr_array_taxa_interaction.txt','w')
outfile.write('array_id\ttaxid\ttaxa_name\tlevel\tsample\tposition\tmapped_cluster\n')
for f in files:
    f = f.strip()
    blast_file = array_path + '/' + f
    blast = open(blast_file,'r')
    for line in blast.readlines()[4:]:
        blast_info = line.split('\t')
        array = blast_info[0]
        array_keys = array.split('__')
        array_key = array_keys[0]+'__'+array_keys[1]
        cluster = blast_info[1]
        length = blast_info[3]
        if taxa_map.has_key(array_key):
            if (length >= 35) and (taxa_map[array_key][3] == f) :
                array_result = '\t'.join(str(n) for n in taxa_map[array_key])
                outfile.write(array+'\t'+array_result+'\t'+cluster+'\n')

outfile.close()
