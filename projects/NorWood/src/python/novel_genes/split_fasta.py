#!/usr/bin/python
'''
	Name: split_fasta.py
	Author: David Sundell
	Contact: davidsundell84@gmail.com
	Modified: 2015-03-17
	Latest: 2015-10-09 - output folder option added
'''

''' Import nessesary libraries '''
import sys

usage=""" 
python split_fasta.py <fasta file> <outfile_prefix> -n nfiles -l nseq -e suffix
	or any file with features starting with ">"
"""

from argparse import ArgumentParser

parser = ArgumentParser(description = """
	python split_fasta.py <fasta file> -s sfiles -n nseq -e suffix
			 or any file with features starting with ">"

	run python split_fasta.py -h to get detailed help
	""")
parser.add_argument('file', metavar="", type=str, help="The source fasta file")
parser.add_argument('prefix', metavar="", type=str, help="The prefix for output fasta files")
parser.add_argument('-s','--sfiles', metavar="", default=False, type=int, help="Choose the number of files to split fasta file in, the script itself determine the since of each file")
parser.add_argument('-n','--nseq', metavar="", default=1000, type=int, help="Choose the number of sequences for each fasta subfile, default is 1000 seq per file")
parser.add_argument('-e','--end', metavar="", default="fna", type=str, help="Change ending of output files default prefix_n_to_n2.(fna)")
parser.add_argument('-o',"--output", metavar="",default=".",type=str, help="Select output folder default current folder")

if len(sys.argv) == 1:
	print usage
	exit()

'''Parse file'''

if __name__=="__main__":
	args = parser.parse_args()

	fin = args.file
	out_pref = args.prefix
	n_seq = False
	if args.sfiles:
		n_seq = sum(1 for line in open(fin))
		split_lines = (int(n_seq/args.sfiles)/2)+1
	elif args.nseq:
		split_lines = args.nseq
		n_seq = sum(1 for line in open(fin))
	file_n = 1
	count = 0
	with open(fin, "r") as f:
		outF = "%s/%s_%i_to_%i.fna" % (args.output,out_pref,file_n,file_n+split_lines-1)
		fout = open(outF,"w")
		for row in f:
			if row.startswith(">"):
				count += 1
			if count % (split_lines+1) == 0:
				fout.close()
				file_n += split_lines
				if n_seq and (n_seq/2 <file_n+split_lines):
					outF =  "%s_%i_to_%i.%s" % (out_pref,file_n,n_seq/2,args.end)
				else:
					outF = "%s_%i_to_%i.%s" % (out_pref,file_n,file_n+split_lines-1,args.end)
				fout = open(outF,"w")
				count +=1
			fout.write(row)
	fout.close()
