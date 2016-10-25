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
python countCDS.py <infile> <outfile>
"""

lim = 0.6
outfile = open(sys.argv[2],"w")
f = open(sys.argv[1],"r")
while True:
	row = f.readline()
	row2 = f.readline()
	if not row or not row2:
		break
	mRNA = row.strip("\n").split("\t")
	CDS = row2.strip("\n").split("\t")
	mlen = 1+int(mRNA[4])-int(mRNA[3])
	Clen = 1+int(CDS[4])-int(CDS[3])
	#print row,
	#print row2,
	quote = float(Clen)/mlen
	#print mlen,"/",Clen, " = ",quote
	if quote > lim:
		outfile.write(row)
		outfile.write(row2)
	