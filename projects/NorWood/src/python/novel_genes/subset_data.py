#!/usr/bin/python
'''
	Name: subset_byname.py
	Author: David Sundell
	Contact: davidsundell84@gmail.com
	Modified: 2015-03-17
'''

''' Import nessesary libraries '''
import sys

usage=""" 
python subset_byname.py <infile> <key_file> <outfile> delimiter="\\t" column=1 
		Observe key file to subset by must have only one column
"""


''' remporary key_match due to mRNA vs gene id usage'''
with open("/mnt/picea/home/david/wood/PASA/20150211/Potri03_update.12239_features.fasta","r") as f:
	translate_dic = {}
	for row in f:
		if row.startswith(">"):
			data = row.strip(">").split(" ")
			k1 = data[0]
			k2 = data[1]
			translate_dic[k1] = k2
			translate_dic[k2] = k1


''' Retrieve system inputs '''
if len(sys.argv) == 1:
	print usage
	exit(0)

try:
	infile = sys.argv[1]
	index_file = sys.argv[2]
	outfile = sys.argv[3]
except:
	print usage
	exit(0)

if len(sys.argv) > 4:
	delimiter = sys.argv[4]
else:
	delimiter = "\t"
if len(sys.argv) > 5:
	column = int(sys.argv[5])-1
else:
	column = 0 ## first column in python has index 0

def read_file(f,delimiter,column):
	key_file = {}
	with open(f, "rU") as fl:
		for row in fl:
			data = row.strip("\n").split(delimiter)
			key = data[column]
			try:
				key_file[key].append(row)
			except KeyError:
				key_file[key] = [row]
	return key_file

if __name__=="__main__":
	## Key file is requested to only have one column
	''' Retrieve keys and datafile '''
	keys = read_file(index_file,delimiter,0).keys()
	data = read_file(infile,delimiter,column)
	dkeys = data.keys()
	dkeys.sort()
	''' Subset file and write to output '''
	out = open(outfile,"w")
	keys.sort()
	not_found = []
	for key in keys:
		try:
			out.write("".join(data[translate_dic[key]]))
		except KeyError:
			not_found.append(key)

	print "Ngenes not found", len(not_found)