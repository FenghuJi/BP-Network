#!/usr/bin/python
import os
import re
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i", "--input",dest="in_file", help="input cut blast and species results",metavar="blast_species_file")
parser.add_option("-o","--output",dest="out_file",help="set the output filename",metavar="output file name")
(options, args) = parser.parse_args()
infile = options.in_file
output = options.out_file

out = open(output,'w')
out.write('contig\ttaxa_id\tsample\tposition\tlevel\tspecies\n')
records = open(infile,'r')
next(records)
for record in records:
        infos = record.split('\t')
        phage = infos[0]

        if bool(re.search('_left',phage)):
                phage = phage.replace('_left','')
        elif bool(re.search('_right',phage)):
                phage = phage.replace('_right','')
        taxid = infos[1]
        sample,position = infos[3].split('/')[8].split('_')
        level,species = infos[4],infos[5].strip()

        out.write(phage+'\t'+taxid+'\t'+sample+'\t'+position+'\t'+level+'\t'+species+'\n')
out.close()
