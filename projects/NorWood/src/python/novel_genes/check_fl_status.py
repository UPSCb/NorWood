#!/usr/bin/python
'''
	Name: check_fl_status.py
	Author: David Sundell
	Contact: davidsundell84@gmail.com
	Modified: 2015-10-01
	Latest update: 2015-12-01
		bug fix, only the first CDS was counted now all CDS are summed up to be coding
'''
## Last modifictaion -% mRNA option added

	

class CheckFullLengthStatus(object):
	"""docstring for CheckFullLengthStatus"""
	def __init__(self):
		super(CheckFullLengthStatus, self).__init__()
		
	def _split_seq(self,string,pices=50):
		string = string.strip("\n").upper()
		while string:
			yield string[:pices]
			string = string[pices:]

	def _merge(self,string):
		return "".join(string.split("\n"))

	def _startCodon(self,codon):
		''' Function check if a codon sequence is the start codon ATG
			Returns True if it is otherwise False
		'''
		if codon.upper() == "ATG":  	## (start)
			return True
		else:
			return False

	def _stopCodon(self,codon):
		'''	Funciton checks if the codon supplied is any of the stop codons
			They are checked in the order of frequence TAA 63% TGA 29% TAG 8% to speed up the check, 
			if not a stop codon False is returned
		'''
		if codon.upper() == "TAA":  	## (amber)
			return True
		elif codon.upper() == "TGA": 	## (ocra)
			return True
		elif codon.upper() == "TAG": 	## (opal)
			return True
		else:
			return False

	def check_fl(self,sequence,CDS,strand="+",pr=False,getStatus=False):
		'''	This function requires a string with a sequence the strand, and a CDS tuple with start and stop'''
		## Initiate variable fl_status, default is False until otherwise proven
		fl_status = False
		## Merge sequence in case the sequence is fasta formatted (with 51 charachters per line)
		sequence = self._merge(sequence)
		## Revers complement in case of minus strand
		bDic = {"A":"T","T":"A","G":"C","C":"G","N":"N"}
		if strand == "-":
			sequence = sequence[::-1]
			sequence = "".join([bDic[x.upper()] for x in sequence])
			start = len(sequence) - CDS[1]
		 	stop = start + (CDS[1]-CDS[0])
		 	startCodon = sequence[start:start+3]
		 	stopCodon = sequence[stop-1:stop+2]
		else:
			start,stop = CDS
			## Fetches the start and end of the CDS -1 is because python counts arrays from 0 so the start index is +1 from the char index
			startCodon = sequence[start-1:start+2]
			stopCodon = sequence[stop-3:stop]

		# if strand == "-":
		# 	sequence = sequence[::-1]
		# 	sequence = "".join([bDic[x.upper()] for x in sequence])
		# 	start,stop = CDS
		# 	startCodon = sequence[start-1:start+2]
		# 	stopCodon = sequence[stop-3:stop]
		# else: 
		# 	start = len(sequence) - CDS[1]
		# 	stop = start + (CDS[1]-CDS[0])
		# 	startCodon = sequence[start:start+3]
		# 	stopCodon = sequence[stop-1:stop+2]

		if strand == "-":
			startCodon,stopCodon = stopCodon,startCodon
		if (strand == "-" and False) or pr:
			print "<----start--->"
			#print sequence
			#if strand == "-":
				#sequence = sequence[::-1]
				#print sequence
			## Get the start and the stop Codon with help of the CDS
			print "CDS: ",[start,stop]
			print "Start: ",startCodon, " status: ", self._startCodon(startCodon), " Stop: ", stopCodon, "status: ",self._stopCodon(stopCodon)
			print sequence
			print len(sequence)
			print sequence[start:start+20]
			print sequence[start-1:stop]
			print "<---end--->"
		

		## Check if there is a start and a stop codon in the start and end respectively if so change fl status to True
		if getStatus:
			return self._startCodon(startCodon),self._stopCodon(stopCodon)
		if self._startCodon(startCodon) and self._stopCodon(stopCodon):
			fl_status = True
			#print startCodon, stopCodon
		## Chop up the sequence in 51 chars per line, to return a fasta formatted sequence. Also return fl status.
		#sequence = "\n".join(_split_seq(full_string,pices=51))
		return fl_status


	def _percent_CDS(self,mRNA,CDS):
		mlen = 1+int(mRNA[1])-int(mRNA[0])
		Clen = 0
		'''Count all CDS and sum length'''
		for i in range(len(CDS)/2):
			Clen += 1+int(CDS[i+1])-int(CDS[i])
		#print "mRNA length: ",mlen
		#print "CDS length:  ",Clen
		quote = float(Clen)/mlen
		#print quote
		return quote





class ReadFrameDP(object):
	"""Read a frameDP merged file, retrieve CDS start:end and sequence"""
	def __init__(self):
		super(ReadFrameDP, self).__init__()

	def _split_seq(self,string,pices=50):
		string = string.strip("\n").upper()
		while string:
			yield string[:pices]
			string = string[pices:]

	def _parse_GFF_line(self,line):
		''' pase a line with GFF info, return a dictionary with info'''
		#temp_gene_10009	FrameD	mRNA	1	253	.	.	.	ID=temp_gene_10009 ; Name=temp_gene_10009
		info = line.strip("#").strip("\n").split("\t")
		_type = info[2] ## mRNA or CDS in this case
		name = info[0]
		strand = info[6]
		pos = (int(info[3]),int(info[4]))
		return name,_type,{"strand":strand, "pos":pos}

	def read_file(self,f):
		''' Reads a sequence file and return a dictionary with sequence info and the sequence
			{
				seq_name: {	seq_head: "> info"
							seq_info: {"mRNA": {"strand":"+/-", "pos":(start:end) }
										"CDS": {"strand":"+/-", "pos":(start:end) }
									}
							sequence: "ATCG"
				}
			}
		'''
		sequences = {}
		_type = None
		name = None
		i = 0
		with open(f, "r") as _f:
			for row in _f:
				#print row
				## If row starts with ## Seq the sequence info will follow
				if row.startswith("##seq"):
					_type = "seq"
					continue
				if row.startswith("##FA"):
					_type = "fa"
					continue
				if _type == "fa_read" and row.startswith(">"):
					_type == "prot"
				#print _type

				## Read fasta sequence when type is "fa" or "fa_read"
				if row.startswith(">") and _type == "fa":
					## The type is change because the next appearance of a > will be the protein sequence
					_type = "fa_read"
					info = row.strip(">").split(" ")
					name = info[0]
					sequences[name]["seq_head"] = row
				elif _type == "fa_read":
					## Append lines with sequence.
					if row.startswith(">"):
						_type = "prot"
						continue
					sequences[name]["sequence"]+="\n".join(self._split_seq(row))+"\n"
				## Read sequence info with GFF format
				if _type == "seq":
					name,type1,info = self._parse_GFF_line(row)
					try:
						sequences[name]["seq_info"][type1] = info
					except KeyError:
						sequences[name] = {
											"seq_head":"",
											"sequence":"",
											"seq_info":{type1:info}
									}
		return sequences


''' Main script '''

if __name__=="__main__":
	''' Import nessesary libraries '''
	from argparse import ArgumentParser

	parser = ArgumentParser(description = """
			Check full length status in frameDP output
		""")
	parser.add_argument('input',type=str, help="The source frameDP merged gff file")
	parser.add_argument('output', type=str, help="File to output results")
	parser.add_argument('-g ', '--gene_file', type=str, help="Add a file to output seq names in")
	parser.add_argument('-fl','--full_length', metavar="", default=False, type=bool, help="Add argument to filter output file by the once that are full length")
	#parser.add_argument('-pc','--pc_mrna', metavar="",default=0.4,type=float,help="Add a percentage of mRNA that is minimum to be considered true fl, 0.4 = 40%")
	parser.add_argument('--dump', metavar="", type=str, help="The dump argument adds a possibility to dump non fl sequences to a separate file if -fl is true")
	parser.add_argument('-pc','--pc_mrna',metavar="",type=float,default=0.4,help="Add a percentage of mRNA that is minimum to be considered true fl, 0.4 = 40 percent")

	args = parser.parse_args()

	rFDP = ReadFrameDP()
	sequences = rFDP.read_file(args.input)
	fl_status = CheckFullLengthStatus()
	seqs = sequences.keys()
	seqs.sort()
	flstat = 0
	fout = open(args.output, "w")
	if args.gene_file:
		gene_out = open(args.gene_file, "w")
		gene_out.write("gene\n")
	for seq in seqs:
		CDS_pc = 0
		seq_obj = sequences[seq]
		c_seq = seq_obj["sequence"]
		pos = seq_obj["seq_info"]["CDS"]["pos"]
		strand = seq_obj["seq_info"]["CDS"]["strand"]
		#if seq_obj["seq_head"].split(">")[1].split(" ")[0] == "temp_gene_6273":
		#	print fl_status.check_fl(c_seq,pos,strand,pr=False)
		if fl_status.check_fl(c_seq,pos,strand):
			if fl_status._percent_CDS(seq_obj["seq_info"]["mRNA"]["pos"],seq_obj["seq_info"]["CDS"]["pos"]) > args.pc_mrna:
				#fout.write("above_40.seqNames")
				#print seq_obj
				_full_length = True
				CDS_pc+=1
			else:
				_full_length = False
				#fout.write("below_40.seqNames")
			#print seq_obj
			flstat += 1
			if args.gene_file:
				if args.full_length:
					if _full_length:
						gene_out.write(seq+"\n")
				else:
					gene_out.write(seq+"\n")
			fout.write(seq_obj["seq_head"])
			fout.write(seq_obj["sequence"])
	print "nr fl status; ", flstat
	print "nr fl > 40% CDS", CDS_pc
	
