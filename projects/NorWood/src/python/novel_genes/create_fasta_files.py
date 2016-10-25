#!/usr/bin/python

'''
	Author: David Sundell
	Contact: david.sundell@umu.se
	Created: 2015-06-29
	Modified: 2015-06-29
	Description:
		Script to create fasta files from the PASA frameDP merged output
		will output three files, a nucleotide fasta file for both coding 
		and non coding genes as well as a amino acid fasta file for the 
		protein coding genes.
'''

import sys

if len(sys.argv) == 1:
	print "usage python create_fasta_files.py filter_file frameDP_peptide_source frameDP_sequence_source"
	exit()
'''Input files'''

# filter
filter_file = sys.argv[1]

# datafiles
pepSource = sys.argv[2]
nukeSource = sys.argv[3]
model_to_potri = sys.argv[4]

def read_filter(infile):
	''' This function takes a file with the coding and non coding genes that are to be saved, and will filter the output
		input file is a list with coding and non coding genes separated by a line of "ncoding" where coding genes are first'''
	coding = {}
	ncoding = {}
	c = True
	with open(infile,"r") as f:
		for row in f:
			gene = row.strip("\n")
			if gene == "ncoding":
				c = False
				continue
			if c:
				coding[gene] = True
			else:
				ncoding[gene] = True
	return coding,ncoding

def read_source(source_file,tp="peptide"):
	if tp == "peptide":
		pass
	elif tp == "sequence":
		pass
	else:
		print "source type not handled please input peptide or sequence source from frameDP"

def read_seq_fasta(in_file):
	'''Read a fasta file split by headers'''
	sequences = {}
	sequence = ""
	with open(in_file,"r") as f:
		for row in f:
			data = row.strip("\n")
			if data.startswith(">"):
				data = data.strip(">")
				if sequence != "":
					sequences[headers[0]] = sequence
				headers = data.split(" ")
				tmp = []
				for x in range(len(headers)):
					if "temp_gene" not in headers[x]:
						if "novel" in headers[x]:
							if "novel_gene" not in headers[x]:
								tmp.append(x)
						else:
							tmp.append(x)
				for x in reversed(tmp):
					del headers[x]
				sequence = ""
			else:
				sequence += data
	sequences[headers[0]] = sequence
	return sequences

def write_protein_file(in_file,coding,gdic=False):
	output = open("output.txt","w")
	with open(in_file,"r") as f:
		for row in f:
			data = row.strip("\n")
			if data.startswith(">"):
				write = False
				data = data.strip(">")
				headers = data.split(":")
				if headers[0] in coding.keys():
					write = True
					gName = headers.pop(0)
					output.write(">"+gName+" "+":".join(headers)+"\n")
			else:
				if write:
					output.write(row)
					

def write_seq_file(coding,ncoding,prot_seq=False,nuke_seq=False,output="output.txt",gdic=False):
	with open(output, "w") as f:
		coding = coding.keys()
		coding.sort()
		ncoding = ncoding.keys()
		ncoding.sort()
		#for gene in coding:
		#	f.write(">"+gdic[gene]+"\n"+nuke_seq[gene]+"\n")
		for gene in ncoding:
			f.write(">"+gdic[gene]+"\n"+nuke_seq[gene]+"\n")

def get_name_dic(f_in):
	gdic = {}
	with open(f_in) as f:
		for row in f:
			data = row.strip("\n").split("\t")
			gene = data[1]
			potri = data[2]
			gdic[gene] = potri
	return gdic

def get_protein_seq(sequence,start=0,reverse=False):
	''' walk through all the possible orfs and pick the longest one '''
	if reverse:
		start+=1
		my_seq = Seq(sequence[:-int(start)],generic_dna)
		my_seq = my_seq.reverse_complement()
	else:	
		my_seq = Seq(sequence[int(start):], generic_dna)
	my_prot = my_seq.translate()
	return (my_prot)
			
''' Main script '''

coding,noncoding = read_filter(filter_file)
#prot_seq = read_seq_fasta(pepSource)
nuke_seq = read_seq_fasta(nukeSource)

gdic = get_name_dic(model_to_potri)
#print nuke_seq.keys()
#write_seq_file(coding,noncoding,nuke_seq=nuke_seq,gdic=gdic)
write_protein_file(pepSource,coding,gdic)






#Old script

"""
from __future__ import print_function

''' Import nessesary libraries '''

from pprint import pprint
import sys
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein

usage=''' python create_fasta_seq.py <gff_file> <fasta_file> <output prefix>
	optional: if output is "examine" the script will 
				output a overview of the features in 
				the gff file and exit.
'''

if len(sys.argv) == 0:
	print(usage)
	sys.exit()

''' Help functions '''
def getSequence(sub_feature, source,strand):
	''' The function takes a sub_feature of exon type, a source file for sequence and the strand (+/-) 
	It returns the bases of the exon matched to the source'''
	start = sub_feature.location.start
	end = sub_feature.location.end
	sequence = ""
	if strand == -1:
		''' reverse complement '''
		sequence = source[start:end].reverse_complement().seq
	else:
		sequence = source[start:end].seq
	return sequence

def read_name_dictionary(infile):
	''' This function takes a name of a dictionary file and returns a temp to final name dictionary
	metatable/model_to_potri_table.txt
	'''
	# novel_model_1_54f5ee74	novel_gene_1_54f5ee74	Potri001C000100	2021876	2024153	01
	dic = {}
	with open(infile,"r") as f:
		for row in f:
			data = row.strip("\n").split("\t")
			gene_id=data[2]
			start = data[3]
			stop = data[4]
			chrom = data[5]
			tmp_id=data[1]
			miRNA_id=data[0]
			dic[gene_id] = {
				"gene_id":gene_id,
				"start":start,
				"stop":stop,
				"chr":chrom,
				"miRNA_id":miRNA_id,
				"tmp_id":tmp_id
			}
	return dic


def coding_and_non_coding_filter(infile):
	''' This function takes a file with the coding and non coding genes that are to be saved, and will filter the output
		input file is a list with coding and non coding genes separated by a line of "ncoding" where coding genes are first'''
	coding = {}
	ncoding = {}
	c = True
	with open(infile,"r") as f:
		for row in f:
			gene = row.strip("\n")
			if gene == "ncoding":
				c = False
			if c:
				coding[gene] = True
			else:
				ncoding[gene] = True
	return coding,ncoding

def get_best_protein_seq(sequence):
	''' walk through all the possible orfs and pick the longest one '''
	ORF = [0,1,2,3,1,2]
	my_seq = Seq(seq[1:], generic_dna)
	mini = 100
	best_prot = ""
	reverse = False
	for x in ORF:
		if x == 3:
			reverse = True
			x = 0
		my_seq = Seq(sequence[x:], generic_dna)
		if reverse:
			my_seq = my_seq.reverse_complement()
		my_prot = my_seq.translate()
		print(my_prot.count("*")),
		if my_prot.count("*") < mini:
			mini = my_prot.count("*")
			best_prot = my_prot
	print(cgene,": "),
	print(best_prot)


def read_merged_fasta(in_file):
	'''Read a fasta file split by headers'''
	sequences = {}
	sequence = ""
	with open(in_file,"r") as f:
		for row in f:
			data = row.strip("\n")
			if data.startswith(">"):
				data = data.strip(">")
				if sequence != "":
					sequences[headers[0]] = sequence
				headers = data.split(" ")
				tmp = []
				for x in range(len(headers)):
					if "temp_gene" not in headers[x]:
						if "novel" in headers[x]:
							if "novel_gene" not in headers[x]:
								tmp.append(x)
						else:
							tmp.append(x)
				for x in reversed(tmp):
					del headers[x]
				sequence = ""
			else:
				sequence += data
	sequences[headers[0]] = sequence
	return sequences

''' Global variables '''
## Set print on and off
pr = True
## Feature to extract
fture = "exon"
###
examine = False

''' Run script '''
if __name__=="__main__":
	if len(sys.argv) == 1:
		print(usage)
		exit()
	in_file = sys.argv[1]
	fasta_file = sys.argv[2]
	gene_dic_file = sys.argv[3]
	filter_file = sys.argv[4]
	if pr:
		print(sys.argv[5])
		#out_file = open(sys.argv[5],"w")
	name_dic = read_name_dictionary(gene_dic_file)
	coding,ncoding = coding_and_non_coding_filter(filter_file)	
	sequences = read_merged_fasta(fasta_file)
	
	for cgene in coding:
		cgenex = cgene.replace(".","")
		seq = sequences[name_dic[cgenex]["tmp_id"]].upper()
		header = ">"+cgene
		get_best_protein_seq(seq)
		
	exit("test row 114")
	













	''' For each chromosome in gff file '''
	for index,record in enumerate(gff):
		#print(record)
		#print(record.id)
		scaffold = records[record.id]
		#print scaffold[37264471:37264847].reverse_complement().seq
		#print "number of records: ",len(record.features)
		''' For each feature in Chromosome'''
		for rec in record.features:
			''' For all mRNA in gene'''
			gene_id = rec.id
			for sf in rec.sub_features:
				''' for each feature in mRNA '''
				#print sf
				mRNA_id = sf.id
				mRNA = []
				strand = 1
				nExon = 0
				''' for all exons within the mRNA'''
				for sf2 in sf.sub_features:
					if sf2.type == fture:
						strand = sf2.strand
						if strand == -1:
							count2+=1
						nExon += 1
						mRNA.append(str(getSequence(sf2,scaffold,strand)))
					else:
						continue
				''' Join sequences for exon and print'''
				if strand == -1:
					'''Build exons in reverse order'''
					string = "".join(reversed(mRNA))
				else:
					string = "".join(mRNA)
				if string != "":
					count += 1
					temp = ""
					#" "+str(strand)+" "
					#header = ">"+mRNA_id+" "+gene_id+" exons: "+str(nExon)+" "+str(record.id)
					if gene_id in coding:
						header = ">"+coding[gene_id]["gene_id"]
						protein = True
					elif gene_id in ncoding:
						protein = False
						header = ">"+ncoding[gene_id]["gene_id"]
					else:
						print("Error occured gene not founc",gene_id)
						exit()
					if pr:
						out_file.write(header)
						print(sf.location,file=out_file)
						out_file.write(string+"\n")
					#print ">",gene_id,mRNA_id,str(nExon),"\n",string
				#else:
					#print(">",gene_id,mRNA_id,str(nExon),"\n")
		print("features found (",record.id,"): ",count)
		print("features on (-) strand: ", count2)
		count2 = 0
		count = 0
	out_file.close()
"""
