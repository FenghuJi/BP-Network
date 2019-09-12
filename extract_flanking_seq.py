#!/usr/bin/python
##Extract flanking sequence for prophage
import os
from optparse import OptionParser
from Bio import SeqIO
parser = OptionParser()
parser.add_option("-d", "--work_diretory",dest="work_dir", help="set working directory",metavar="work_path")
parser.add_option("-l","--filelist",dest="filelist",help="set filelist of the project",metavar="project_filelist")
(options, args) = parser.parse_args()
work_dir = options.work_dir
filelist = options.filelist

def get_position(seq_name):   ##From prophage name,we can get its contig name,start and stop
	temp_contig = seq_name[10:].split('_')
	temp_contig[6] = temp_contig[6].split('-')[0]
	contig=temp_contig[0]
	for i in range(1,6):
		contig=contig+'_'+ temp_contig[i]
	contig= contig +'.'+ temp_contig[6]
	start=int(temp_contig[-2].split('-')[1])
	stop=int(temp_contig[-2].split('-')[2])
	return contig,start,stop




def get_flanking_seq(seq,start,stop):
	seq_len = len(seq)
	left=''
	right=''
	if start >= 250:
		left = seq[49:249]
	elif ((start < 250) and (start > 200)):
		left = seq[0:start - 51]
	if (stop + 250) < seq_len:
		right = seq[seq_len-251:seq_len-51]
	elif (stop + 250) > seq_len:
		right = seq[stop+49:seq_len-1]
	return left,right


with open(filelist,'r') as files:
	for file in files:  ## Loop for every sample
		file=file.strip()

		contig_fasta=work_dir+'/01_fasta/'+file+'.fasta'  ###Here we need to read the fasta for extract flanking sequence
		
        out_dir = work_dir + '/04_flank_seq'
		pro_path = work_dir + "/03_prophage_seq/" + file
		leftFile = open(out_dir + '/' + file + '_left','w')
		rightFile = open(out_dir + '/' + file + '_right','w') 

		seqs_dict = {} ## ids mapping seqs
		contig_fin = open(contig_fasta,'r')
		for record in SeqIO.parse(contig_fin,'fasta'):
			fa_id = record.id
			fa_seq = record.seq
			seqs_dict[fa_id] = fa_seq

        with open(pro_path,'r') as prophages:
            for prophage in prophages:
                prophage = prophage.strip()
                contig,start,stop = get_position(prophage)
                seq = seqs_dict[contig]
                left,right = get_flanking_seq(seq,start,stop)
                if len(left)!=0:
                	leftFile.write(">"+ prophage + 'left' +'\n')
                	leftFile.write(str(left + '\n'))
                if len(right)!=0:
                	rightFile.write(">" + prophage + 'right' +'\n')
                	rightFile.write(str(right + '\n'))

        seqs_dict.clear()
        leftFile.close()
        rightFile.close()
