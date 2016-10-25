#!/usr/bin/python

usage="""fragmentNames.py input output"""

description="""
	This script is taking a file with gene names of kind Potri.001[NC]000100
	and will translate them to fragment names Potri.001F000100, any possible 
	name match will be corrected from the coding number
"""
import sys
if len(sys.argv) == 1:
	print usage
	print description
	exit()

def checkDup(duplicates, name):
	integer = int(name[-6:])
	if integer <= 100:
		integer+=1000
	if name[:-6]+str(integer-100).zfill(6) not in duplicates:
		return name[:-6]+str(integer-100).zfill(6)
	else:
		return checkDup(duplicates,name[:-6]+str(integer-100).zfill(6))

outfile = open(sys.argv[2],"w")
dupl = []
with open(sys.argv[1],"r") as f:
	columns = f.readline().strip("\n").split("\t")
	columns.append("fragment.id")
	print >> outfile,"\t".join(columns)
	for row in f:
		names = row.strip("\n").split("\t")
		if "Potri" in names[0] and "TRUE" in names[1] and "FALSE" in names[2] and "fragment" in names[-1]:
			print "found"
			names.append(names[0].replace("N","F").replace("C","F"))
			if names[-1] not in dupl:
				dupl.append(names[-1])
			else:
				print "WARNING: duplicate found",names
				newName = checkDup(dupl,names[-1])
				if newName:
					print names,newName
					names[-1] = newName
					dupl.append(newName)
				else:
					names[-1] = "NA"
		else:
			names.append("NA")
		print >> outfile,"\t".join(names)
				
