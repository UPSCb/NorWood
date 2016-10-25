#!/usr/bin/python
'''
	Name: merge_fa_and_frameDP_gff.py
	Author: David Sundell
	Contact: davidsundell84@gmail.com
	Modified: 2015-03-17
'''

''' Import nessesary libraries '''
import sys

usage=""" 
python merge_fa_and_frameDP_gff.py <fasta file> <merge file from frameDP> <out>
	or any file with features starting with ">"
"""

from argparse import ArgumentParser

parser = ArgumentParser(description = """
	python merge_fa_and_frameDP_gff.py <fasta file> <merge file from frameDP> <out>
			 or any file with features starting with ">"

	run python merge_fa_and_frameDP_gff.py -h to get detailed help
	""")
parser.add_argument('file', metavar="", type=str, help="The source fasta file")
parser.add_argument('merge', metavar="", type=str, help="Choose the number of files to split fasta file in, the script itself determine the since of each file")
parser.add_argument('out', metavar="", type=str, help="Choose the number of sequences for each fasta subfile, default is 1000 lines per file")
parser.add_argument('--gff3', metavar="", default=False, type=str, help="Add one additonal output file, for a merged gff3 file to be written")

if len(sys.argv) == 1:
	print usage
	exit()

args = parser.parse_args()
	
'''Parse file'''

if args.gff3:
	gff_file=open(args.gff3,"w")

def read_fasta(file):
	'''Reads a fasta file where new sequence starts with >'''
	### Maybe not needed, file can be written while it is walked through

def _split_seq(string,pices=50):
	string = string.strip("\n")
	while string:
		yield string[:pices]
		string = string[pices:]

def read_gff_fasta(fin):
	'''Reads a frameDP fasta file with gff headers, preferably merge gff files with cat *.gff3 > all.gff3'''
	read = False
	name=False
	gff=False
	seq_dic = {}
	string = ""
	with open(fin, "r") as f:
		for row in f:
			if row.startswith("##FASTA"):
				read = True
				start = True
				gff = False
				# store previous feature to its name
				if name:
					seq_dic[name] = string
				## Continue is used to not include the line with ##FASTA in the output
				continue
			if row.startswith("##seq"):
				gff = True
			elif gff and args.gff3:
				gff_file.write(row)
			if read and row.startswith(">"):
				## Fetch the next > rowstart and end read, this is where the protein is written
				if row.startswith(">") and not start:
					read = False
					continue
				start = False
				### Split row, fetch the name and store the comming lines in array
				#print row
				name = row.strip(">").split(" ")[0]
				string = row
			elif read:
				if len(row) > 51:
					row = "\n".join(_split_seq(row))+"\n"
				string += row
	seq_dic[name] = string
	keys = seq_dic.keys()
	return seq_dic
	

if __name__=="__main__":
	fin = args.file
	fout = args.out
	merge = args.merge
	updated = 0
	updated_seq = read_gff_fasta(merge)
	with open(fin, "r") as f:
		fout = open(fout,"w")
		for row in f:
			## Set update to false, if updated it will be turned to True
			update = False
			if row.startswith(">"):
				name = row.strip(">").split(" ")[0]
				try:
					string = updated_seq[name]
					#string = row.strip("\n") + "\t" +string
					fout.write(string)
					update = True
					write = False
					updated += 1
				except KeyError:
					## No update exist continue
					write = True
					fout.write(row)
			elif not update and write:
				if len(row) > 50:
					row = "\n".join(_split_seq(row.upper()))
				fout.write(row.upper().strip("\n")+"\n")
	print "Updated: ", str(updated)
	fout.close()

if args.gff3:
	gff_file.close()
