#!/usr/bin/python
'''
	Name: novel_feature_names.py
	Author: David Sundell
	Contact: davidsundell84@gmail.com
	Modified: 2015-04-08
'''
import sys
''' Import nessesary libraries '''
from argparse import ArgumentParser

parser = ArgumentParser(description = """
		Name novel features with proper names
		python novel_feature_names.py 
			-g Potri03_update.12239_orig_features.fasta (to get information about chr and pos)
			-d output directory
			-o output_file
			-ab Potri_C,Potri_N
			-f all_novel_coding_genes.txt,novel_non_coding_genes.txt
	""")
parser.add_argument('-g','--gff',type=str, help="The source frameDP merged gff file")
parser.add_argument('-d','--dir', type=str, help="The output directory")
parser.add_argument('-o','--output', type=str, help="File to output results")
parser.add_argument('-ab','--abbreviations', type=str, help="abbreviations separated by comma, supply --files with same order as abbreviations")
parser.add_argument('-f','--files', type=str, help="Input files with gene names in rows, files separated by comma")
parser.add_argument('-v', '--verbose', type=bool, default=False)

def __parse_fasta(fasta,name_dict):
	''' Parse fasta file created in the long non coding RNA pipeline, it should look something like
		>temp_gene_310 temp_model_4211 exons: 1 Chr01[4259299:4259634](+)
	'''
	dic = {}
	scaffold = {}
	novel_features = name_dict
	with open(fasta, "r") as f:
		for row in f:
			if row.startswith(">"):
				info = row.strip(">").split(" ")
				gname = info[0]
				mname = info[1]
				Chr = info[4].split("[")[0]
				start = int(info[4].split(":")[0].split("[")[1])
				end = int(info[4].split(":")[1].split("]")[0])
				length = max([start,end]) - min([start,end])
				inf = {"chr":Chr,"start":start,"end":end,"length":length}
				dic[gname] = inf
				dic[mname] = inf
				if gname in novel_features or mname in novel_features:
					try:
						scaffold[Chr].append([mname,start,end,length,gname])
					except KeyError:
						''' Add new scaffold '''
						scaffold[Chr] = [[mname,start,end,length,gname]]

	for sc in scaffold.keys(): 
		scaffold[sc] = sorted(scaffold[sc],key=lambda x: x[1])
	return scaffold

def __read_genes(gene_file,directory=""):
	'''Simply read a file with genes by row'''
	genes = {}
	with open(directory+gene_file, "r") as f:
		for row in f:
			genes[row.strip("\n")] = True
	return genes

def __create_name_table(name_dicts,abbreviations):
	'''Create table with new gene names'''
	pass

def __sort_name_dict(name_dict):
	'''Sort name dict by chromosome, start, length'''

def __print_gene(feature,abbreviations,Chr_int,g_int):
	if True:
		print "%s\t%s\t%s%03d%s%06d\t%d\t%d\t%s" % (feature[0],feature[-1],abbreviations[k][0],Chr_int,abbreviations[k][1],g_int,feature[1],feature[2],"Chr_"sc_i) 

def __print_scaffold(feature,abbreviations,Chr_int,g_int,sc_i):
	if True:
		print "%s\t%s\t%s%s%06d\t%d\t%d\t%s" % (feature[0],feature[-1],abbreviations[k][0],abbreviations[k][1],g_int,feature[1],feature[2],"scaffold_"+sc_i) 


''' Main script '''
import re
			
if __name__=="__main__":
	args = parser.parse_args()
	#print args
	name_dicts = []
	for gfile in args.files.split(","):
		name_dicts.append(__read_genes(gfile,args.dir))
	if args.verbose:
		print "Dictionaries read:"
		for x in range(len(name_dicts)):
			print "\t",args.files.split(",")[x],len(name_dicts[x])
	abbreviations = args.abbreviations.split(",")
	if args.verbose:
		print "Abbreviations given: "
	for x in range(len(abbreviations)):
		abbreviations[x] = abbreviations[x].split("_")
		if args.verbose:
			print "\t",abbreviations[x]
	'''Read gff file'''
	#print len(scaffolds["Chr01"])
	Chr_int = 0
	for k in range(len(abbreviations)):
		scaffolds =__parse_fasta(args.gff,name_dicts[k])
		scaff = scaffolds.keys()
		scaff.sort()
		chromosome = True
		for sc in scaff:
			sc_i = re.findall(r'\d+',sc)[0]
			Chr_int = int(sc_i)

			if Chr_int > 19 and chromosome:
				chromosome = False
				g_int = 100

			if chromosome:
				g_int = 100	
			x = 0
			for feature in scaffolds[sc]:
				try:
					if chromosome:
						__print_gene(feature,abbreviations,Chr_int,g_int)
					else:
						__print_scaffold(feature,abbreviations,Chr_int,g_int,sc_i)
					x += 1
					g_int += 100
				except KeyError:
					print "---"
					sys.exit("Error, gene not found")
		Chr_int += 1

	'''Print coding gene file'''
	#for genes in 

	'''Print non coding gene file'''

	
