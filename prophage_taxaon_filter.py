#!/usr/bin/python
# Usage: python get_unique_taxa.py -i ncbi_taxon.txt
## For use , it's needed to sort the ncbi_taxon file at first
## bash : {---- less prophage_cut_flanking_to_bacterial_ncbi_taxon_id.txt |sort -k1|uniq|cut -f1,2,4,5,6|sed 's/_left//g'|sed 's/_right//g' > filter-ncbi-id.txt -----} 
## python ../../scripts/prophage_taxaon_filter.py -i filter-ncbi-id.txt
import os
from collections import Counter
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i", "--input",dest="infile", help="set ncbi_taxon",metavar="ncbi taxon assignment")
(options, args) = parser.parse_args()
infile = options.infile


if __name__ == '__main__':
	outfile = 'filter_taxon_result.txt'
	OUT = open(outfile,'w')
	most_common_result = {}

	taxa_info = open(infile,'r')
	next(taxa_info)
	for infos in taxa_info:
		mapping = infos.split('\t')
		contig = mapping[0].strip()
		taxid = mapping[1].strip()
		sample = mapping[2].strip().replace('/mnt/raid1/wangteng/phage-bacteria_interactions/interations/jfh/03_cut/','')
		level = mapping[3].strip()
		taxa = mapping[4].strip()
		result = "-".join([taxid,sample,level,taxa])
		most_common_result.setdefault(contig,[]).append(result)
		
	for key in most_common_result.keys():
		taxids = most_common_result[key]
		taxa_result = Counter(taxids).most_common(1)[0][0]
		taxid = taxa_result.split('-')[0]
		sample = taxa_result.split('-')[1]
		level = taxa_result.split('-')[2]
		final_assign = taxa_result.split('-')[3]
		OUT.write(key+'\t'+ taxid +'\t'+ sample+'\t'+level+'\t'+final_assign+'\n')
	OUT.close()