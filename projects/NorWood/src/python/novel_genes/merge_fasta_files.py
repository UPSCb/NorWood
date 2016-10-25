#/usr/bin/python

import sys
from os import walk
import os

if len(sys.argv) == 1 or sys.argv[1] == "-h":
	print "Usage: python merge_fasta_files.py workdir outputfile"


mypath = sys.argv[1]
output = sys.argv[2]
## emtpy output file
open(output,"w").close()

def parseFrameDPFasta(path):
	read = False
	fasta = False
	'''function to read the sequence from a frameDP gff file'''
	outfile = open(output,"a")
	with open(path) as tmp_f:
		#for row in tmp_f:
			# if read and fasta and row.startswith(">"):
			# 	read = False
			# 	fasta = False
			# x	#exit()
			# if fasta:
			# 	outfile.write(row)
			# if row.startswith("##FASTA"):
			# 	read = True
			# if read and row.startswith(">"):
			# 	fasta = True
			# 	outfile.write(row)
		outfile.write(tmp_f.read())
	outfile.close()
	return


f=[]
mainDir = next(walk(mypath))[1]
for d in mainDir:
	subDir = next(walk(mypath+"/"+d))[1]
	for sub_d in subDir:
		files = next(walk("/".join([mypath,d,sub_d])))[2]
		for f in files:
			if f.endswith("gff3"):
				parseFrameDPFasta("/".join([mypath,d,sub_d,f]))