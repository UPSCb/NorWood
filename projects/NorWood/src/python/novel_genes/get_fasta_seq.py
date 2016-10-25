#!/usr/bin/python
from __future__ import print_function
'''
	Name: get_fasta_seq.py
	Author: David Sundell
	Contact: davidsundell84@gmail.com
	Modified: 2015-03-17
'''

''' Import nessesary libraries '''

from pprint import pprint
import sys
from BCBio import GFF
from Bio import SeqIO

usage=""" python get_fasta_seq.py <gff_file> <fasta_file> <output>
	optional: if output is "examine" the script will 
				output a overview of the features in 
				the gff file and exit.
"""

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
	in_file = sys.argv[1]
	fasta_file = sys.argv[2]
	if sys.argv[3] == "examine":
		examine = True
	if pr:
		out_file = open(sys.argv[3],"w")
	if examine:
		examiner = GFF.GFFExaminer()
		gff_file = open(in_file, "r")
		pprint(examiner.parent_child_map(gff_file))
		pprint(examiner.available_limits(gff_file))
		gff_file.close()
		exit(0)

	records = SeqIO.to_dict(SeqIO.parse(open(fasta_file,",rU"),"fasta"))
	gffDict = {}
	f = open(in_file,"r")
	#print(f.readline())
	gff = GFF.parse(f,base_dict=records)
	count = 0
	count2 = 0
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
				#print(gene_id)
				coding = True
				if "N" in gene_id:
					#print("Non coding")
					coding = False
					strand = sf.strand
					if strand == -1:
						count2+=1
					mRNA.append(str(getSequence(sf,scaffold,strand)))
					#print(mRNA[-1])
					nExon = 0
				if coding:
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
					header = ">"+mRNA_id+" "+gene_id+" exons: "+str(nExon)+" "+str(record.id)
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
