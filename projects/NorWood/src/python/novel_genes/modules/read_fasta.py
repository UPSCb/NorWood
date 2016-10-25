#!/usr/bin/python

'''
remember to make all separations in the file tab separated, use following command in unix commandline

awk '{print $0}' FrameDPoutput.gff3 OFS="\t" > tab-sep_output.gff3

'''

class ReadFasta(object):
	"""docstring for ReadFasta"""
	def __init__(self):
		super(ReadFasta, self).__init__()
		
	def readFasta(self,infile):
		''' Read a fasta formatted file and return a dictionary with > name and sequence '''
		sequences = {}
		with open(infile,"r") as f:
			for row in f:
				if row_in.startswith(">"):
					name = row.strip(">").strip("\n").split(" ")[0]
					name = name.split("|")[0]
					sequences[name] = ""
				else:
					sequences[name] += row.strip("\n")
		return sequences


if __name__ == "__main__":
	out = open("Potri.count", "w")
	for gene in cds.keys():
		try:
			l_tran = len(tran[gene])
			l_cds = len(cds[gene])
		except KeyError:
			pass
		quote = float(l_cds)/l_tran
		out.write(gene+"\t"+str(quote)+"\t"+str(l_cds)+"\n")

	out.close()




class GeneFeature(object):
	"""docstring for GeneFeature"""
	def __init__(self,file_obj):
		super(GeneFeature, self).__init__()
		self.file_obj = open(file_obj,"r")
		self.read_seq = False
		self.FASTA = False
		self.read_prot = False
		self.read_gff = False
		#print self.file_obj.readline()

	def read_feature(self):
		'''Read a feature from a frameDP source file ,including gff prot seq and genom seq'''
		#### Get the sequence region
		feature = {}
		row_in = self.file_obj.readline()
		while row_in:
			row = row_in.strip("\n")
			if row_in.startswith("##sequence-region"):
				self.read_gff = True
				row = row.split(" ")
				feature["sequence-region"] = {"id":row[-3],"start":int(row[-2]),"end":int(row[-1]),"mRNA":[],"CDS":[]}
			elif row_in.startswith("##FASTA"):
				self.FASTA = True
				self.read_gff = False
			elif self.read_gff:
				row = row.split("\t")
				row[3] = int(row[3])
				row[4] = int(row[4])
				#if row[0] == "temp_gene_3995":
				#	print row
				feature["sequence-region"][row[2]].append(row)
			elif row_in.startswith(">") and self.FASTA:
				self.read_seq = True
				self.FASTA = False
				row = row.strip(">").split(" ")
				_id = row.pop(0) ## collect id
				try:
					_len = int(row.pop(0).split("=")[-1]) ## store length
				except ValueError:
					_len = None
				row.pop(0) ## remove def=
				row.pop(0) ## remove ;

				'''Collect the frameDP comment, if it was modified'''
				frameDP = ""
				if len(row) > 0:
					frameDP+= " ".join(row)
				feature["FASTA"] = {"id":_id,"len":_len,"frameDP":frameDP,"sequence":""}
			
			elif self.read_seq:
				while not row.startswith(">"):
					feature["FASTA"]["sequence"]+=row
					row = self.file_obj.readline().strip("\n")
				self.read_prot = True
				self.read_seq = False
				row = row.split(":")
				_id = row[0][1:]
				_start = int(row[1])
				_end = int(row[2])
				_strand = row[3]
				feature["prot"] = {"id":_id,"start":_start,"end":_end,"strand":_strand,"sequence":""}

			elif self.read_prot and not row_in.startswith(">"):
				feature["prot"]["sequence"]+=row
				self.read_prot = False
			
			## read new row
			row_in = self.file_obj.readline()
			if row_in.startswith("##gff-version"):
				break
		return feature
	
class ReadFrameDP(object):
	"""docstring for ReadFrameDP"""
	def __init__(self, gff_source):
		super(ReadFrameDP, self).__init__()
		self.GF = GeneFeature(gff_source)
		
	def __iter__(self):
		while True:
			feature = self.GF.read_feature()
			if not feature:
				break
			yield feature

'''
testing

from read_fasta import ReadFrameDP
testObj = ReadFrameDP("/mnt/picea/home/david/wood/PASA/20150211/frameDP/gff_files/Wood-frameDP-tab-sep.gff3")
for x in testObj: print x

from read_fasta import GeneFeature
g = GeneFeature("/mnt/picea/home/david/wood/PASA/20150211/frameDP/gff_files/Wood-frameDP-tab-sep.gff3")
g.read_feature()

'''



		