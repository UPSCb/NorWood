#!/usr/bin/python
#from __future__ import print_function
'''
	Name: get_orfs.py
	Author: David Sundell
	Contact: davidsundell84@gmail.com
	Modified: 2015-03-17
'''

info=""" 
		This script will count open reading frames in a fasta source file, 
			count ORFs in that file and write ORFs and positions to a file, 
			if requested the orfs can be outputed as a fasta or gff3 file.
"""


''' Import nessesary libraries '''

from argparse import ArgumentParser

parser = ArgumentParser(description = """
		count open reading frames in fasta sequences
	""")

parser.add_argument('input',type=str, help="Fasta file with sequences where ORFs should be counted")
parser.add_argument('--gff3', action="store_true", help="output ORFs as a gff3 file")
parser.add_argument('--fasta', action="store_true", help="output ORFs as a fasta file")
parser.add_argument('-s','--sort', action="store_true", help="output the result sorted, by transcript name")
parser.add_argument('--sort_file', type=str, default=False, help="File from which print order is read")
parser.add_argument('--reverse', action="store_true", help="Add count on reverse strand to output")
parser.add_argument('-o', '--output', type=str, help="If --gff3 or --fasta or --merge is specified specify output file name")

args = parser.parse_args()

from pprint import pprint
import sys
from BCBio import GFF
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein

sequences = {"seq":"ATGCATCTAAGAATTTAATTGCTGCATTAGTTTAATATGCATCGAATGCTGCATTAG"}

#seq = "".join(seq.split("\n"))

#seq = seq[0:3381]

def get_sequence(f):
	'''Function to get sequences from a fasta file'''
	sequences = {}
	order = []
	with open(f,"r") as fin:
		for row in fin:
			if row.startswith(">"):
				seq_name = row.split(">")[1].split(" ")[0]
				sequences[seq_name] = ""
			else:
				sequences[seq_name]+=row.strip("\n").upper()
			order.append(seq_name)
	return sequences,order

def get_sort(f):
	order = []
	with open(f, "r") as fin:
		fin.readline()
		for row in fin:
			name = row.strip("\n").split("\t")[0]
			order.append(name.lower())
	return order

def find_orfs(sequence,positions):
	pass

def find_codon(my_dna,codon):
	count = my_dna.count(codon)
	start = 0
	all_pos = []
	for x in range(count):
		pos = my_dna.find(codon,start)
		all_pos.append(pos)
		start = pos+1
	return all_pos

def countORFS(start_c,stop_c):
	ORFS = []
	for start in start_c:
		index = 0
		if len(stop_c) <= index:
			break
		while True:
			if len(stop_c) <= index:
				break
			stop = stop_c[index]
			'''Check so that the start and the stop is in the same frame otw continue'''
			if stop > start and (stop-start) % 3 == 0:	
				frame = (start % 3)
				if frame == 0:
					frame = 3
				## the script does not count the stop codon, which is normally counted, therefore the number will be +1
				length = ((stop-start)/3)+1
				length = stop-start
				## The cutoff for the min length ORF is the shortest know protein of 21 amino acids.
				if length > 20:
					ORFS.append((start,stop,frame,length))
				#stop_f = stop
				break
			index +=1
			#else:
			#	print "not same frame"
	return ORFS

def count_CDS(sequences,order,args):
	count = 0
	print "name\tORFs/100bp\tlongestORF\tORF\\%mrna\tnORFs"
	if args.sort:
		loop = sequences.keys()
		if args.sort_file:
			loop = get_sort(args.sort_file)
	else:
		loop = order
	for seq in loop:
		my_dna = Seq(sequences[seq],generic_dna)
		#print my_dna
		start_c = find_codon(my_dna,"ATG")
		stop_c = find_codon(my_dna,"TAA")
		stop_c = find_codon(my_dna,"TGA") + stop_c
		stop_c = find_codon(my_dna,"TAG") + stop_c
		stop_c.sort()
		if args.reverse:
			my_dna_r = my_dna.reverse_complement()

			start_c_r = find_codon(my_dna_r,"ATG")
			stop_c_r = find_codon(my_dna_r,"TAA")
			stop_c_r = find_codon(my_dna_r,"TGA") + stop_c
			stop_c_r = find_codon(my_dna_r,"TAG") + stop_c
			stop_c_r.sort()
		#print start_c
		#print stop_c
		#exit()

		ORFS = []
		#stop_f = 0
		ORFS = countORFS(start_c,stop_c)
		if args.reverse:
			ORFS_r = countORFS(start_c_r,stop_c_r)
			ORFS += ORFS_r
			#for orf_i in range(len(ORFS)):
			#	print ORFS[orf_i][3]," + ",ORFS_r[orf_i][3]," = ", (ORFS[orf_i][3] + ORFS_r[orf_i][3])
			#	exit()
			#	ORFS[orf_i][3] = ORFS[orf_i][3] + ORFS_r[orf_i][3]
		#print ORFS
		import operator
		ORFS = sorted(ORFS,key=operator.itemgetter(3))
		#print ORFS
		if len(ORFS) > 0:
			seq_name = seq
			norfs = round((len(ORFS)*100.0)/(len(my_dna)*3*6),2)
			norfs2 = len(ORFS)
			longestORF = int(ORFS[-1][3])+3
			quote = round(float(longestORF)/len(my_dna),2)
			if quote >= 0.0:
				print seq.strip("\n"),"\t" ,norfs,"\t",longestORF,"\t", quote, "\t",norfs2
		else:
			print seq.strip("\n"),"\t0\t0\t0\t0"
			count +=1
		#for x in ORFS:
		#	print "ORF Length: ", (x[1]-x[0])/3, " Pos: ", x
	#print count

### merge file command
# paste /mnt/picea/home/david/wood/cpat/Potri_non_novel_genes.tsv temp_res/all_result_no_sort_plus_reverse.txt | cut -f1,2,3,4,5,6,8,9,10,11,12 > /mnt/picea/home/david/wood/metatable/Potri_meta_table.tsv

if __name__=="__main__":
	seqs,order = get_sequence(args.input)
	count_CDS(seqs,order,args)
