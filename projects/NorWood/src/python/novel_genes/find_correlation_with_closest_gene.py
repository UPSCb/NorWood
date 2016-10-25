#!/usr/bin/python

'''

Script is made to take a file with two fasta attribute columns and then find the correlation associated with those two genes
has currently to be run at server with the databases, might be written to be run with source files.

'''

import argparse
import GFF.Attributes as _Attributes
#import db as db

parser = argparse.ArgumentParser(description='Find correlation',add_help=False)
req = parser.add_argument_group('required arguments')
req.add_argument('-i','--input', type=str, help='The closest distance file', metavar='')
req.add_argument('-o','--output',type=str, help='The output file',metavar='')

opt = parser.add_argument_group('optional arguments')
opt.add_argument('-t','--table', default="wood",type=str,metavar='', help="Select correlation table default = Correlation_wood_clr")
opt.add_argument('-db','--database', default="popgenieDB_V3", type=str ,metavar='', help="Select database")
opt.add_argument('-h','--help', action='help', help='show this help message and exit')
args = parser.parse_args()

#db_in = db.dbConnect(args.database)

def _get_corr(source,target):
	corr = 0
	query = '''SELECT corr FROM Correlation_%s_clr WHERE gene_i1 = %s AND gene_i2 = %s''' % (args.table,source,target)
	print query
	return corr

def _get_dic():
	dic = {}
	query = "SELECT * FROM GeneTable"

	return dic

def _get_avg_expression(gene):
	query = '''SELECT AVG(expression) FROM Expression_%s''' % (args.table)
	exprs = 0

	return exprs

def Run(corr_file,output_file):
	'''read input file'''
	output = open(output_file,"w")
	print >> output,"\t".join(["source_id","target_id","distance","corr","s_avg_exprs","t_avg_exprs"])
	Attributes = _Attributes.Attributes()
	dic = _get_dic()
	with open(corr_file, "r") as f:
		for row in f:
			info =row.strip("\n").split("\t")
			## If first column is chr remove
			if info[0][0:3] == "Chr":
				info.pop(0)
			source = info[0]
			target = info[1]
			distance = int(info[2])
			source_id = Attributes.read_attributes(source)
			target_id = Attributes.read_attributes(target)
			print source_id
			print target_id
			#source_i,target_i = dic[source_id],dic[target_id]
			corr = _get_corr(source_id,target_id)
			s_exprs = _get_avg_expression(source_id)
			t_exprs = _get_avg_expression(target_id)
			if corr:
				output.write("\t".join([source_id,target_id,str(distance),str(corr),str(s_exprs),str(t_exprs)])+"\n")
			exit()

if __name__=="__main__":
	Run(args.input,args.output)