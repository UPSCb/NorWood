import sys
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from modules.read_fasta import ReadFrameDP
from argparse import ArgumentParser
import os
import subprocess

parser = ArgumentParser(description = """
	python create_gff_final_file.py <fasta file> <merge file from frameDP> <out>
			 or any file with features starting with ">"

	run python create_gff_final_file.py -h to get detailed help
	""")
parser.add_argument('input', metavar="", type=str, help="a PASA? gff3")
parser.add_argument('-g','--gene_annotation', metavar="", default=False,type=str, help="Chose a file if indexes needs to be updated")
parser.add_argument('-gf','--gene_filter', metavar="",default=[], type=str,help="Filter the final set of genes from the gene annotation file")
parser.add_argument('-f', '--frameDP', type=str,help="The source frameDP merged in_file")
parser.add_argument('-o','--output', metavar="", default=False, type=str, help="Add one additonal output file, for a merged gff3 file to be written")
parser.add_argument('--create_genes',action='store_true',help="Add to create new gene features for frameDP updates instead of adding them as an mRNA to existing")
parser.add_argument('--replace_genes',action='store_true',help="Add to replace frameDP genes instead of adding them, use in combination with --create_genes")

args = parser.parse_args()
'''Parse features from file '''

### filter to remove the 38 overlapping genes  temporary
#ovlp = [gene.strip("\n") for gene in open("/mnt/picea/home/david/wood/2015_12_02_final_annotations_38_genes_removed/removed_genes/overlapping_genes.txt")]
#ovlp = False:
from check_fl_status import CheckFullLengthStatus
def checkStartStop(sequence,CDS,gene_id=False):
	'''
		Take a sequence and check if it contains a start or a stop codon  
		code nothing:0, start:1, stop:2, rstart:3, rstop:4
	'''
	checkFL = CheckFullLengthStatus()

	if gene_id == "":
		pstart,pstop = checkFL.check_fl(sequence,CDS,"+",getStatus=True,pr=True)
		rstart,rstop = checkFL.check_fl(sequence,CDS,"-",getStatus=True,pr=True)
		print pstart,pstop
		print rstart,rstop
		exit()
	else:
		pstart,pstop = checkFL.check_fl(sequence,CDS,"+",getStatus=True)
		rstart,rstop = checkFL.check_fl(sequence,CDS,"-",getStatus=True)
	#print pstart,pstop,rstart,rstop
	if pstart: code = 1
	elif pstop: code = 2
	elif rstart: code = 3
	elif rstop: code = 4
	else: code = 5
	return code

def create_sub_feature(mRNAs,CDSs,model_to_id,genes=[],create_gene=False,sequence=False):
	#print genes
	print mRNAs
	CDS_change = {'+':1,'-':-1}
	top_feature = []
	sub_features = []
	sub_qualifiers_base = {"source":"frameDP"}
	gene_id_end = 20
	for i in range(len(mRNAs)):
		mRNA = mRNAs[i]
		gene_id = mRNA[0]
		if model_to_id:
			if not genes:
				genes = model_to_id.keys()
			if gene_id in genes:
				gene_id = model_to_id[gene_id]
		else:
			gene_id = mRNA[0]
		if create_gene:
			#gene_id = gene_id[:-2]+str(gene_id_end)
			#gene_id_end += 20
			gene_qualifiers = sub_qualifiers_base.copy()
			gene_qualifiers["ID"] = [gene_id]
			gene_qualifiers["Name"] = [gene_id]
			rec_info = {"start":mRNA[3],"end":mRNA[4],"type":"gene","strand":CDS_change[mRNA[6]],"id":mRNA[0]}
			gene_feature = new_feature(rec_info,gene_qualifiers)
		mRNA_id = gene_id+"."+str(i+1)
		rec_info = {"start":mRNA[3],"end":mRNA[4],"type":"mRNA","strand":CDS_change[mRNA[6]],"id":mRNA_id}
		sub_qualifiers = sub_qualifiers_base.copy()
		sub_qualifiers["ID"] = [mRNA_id]
		sub_qualifiers["Parent"] = [gene_id]
		top_feature = new_feature(rec_info,sub_qualifiers)
		#sub_feature =new_feature(rec_info,sub_qualifiers)
		index = 1
		#print mRNAs
		for x in range(len(CDSs)):
			fragment = False
			CDS = CDSs[x]
			#print CDS
			'''Function to create sub features from mRNA and CDS'''
			rec_info = {"start":CDS[3],"end":CDS[4],"type":"exon","strand":CDS_change[CDS[6]],"id":mRNA_id+".exon"+str(index)}
			sub_qualifiers1 = sub_qualifiers_base.copy()
			sub_qualifiers1["ID"] = [mRNA_id+".exon"+str(index)]
			sub_qualifiers1["Parent"] = [mRNA_id]
			if CDS[-1].split(" ; ")[-1].split("=")[-1] == "fragment":
				fragment=True
				sub_qualifiers1["Note"] = "fragment"
			sub_features.append(new_feature(rec_info,sub_qualifiers1))
			if model_to_id:
				if CDS[0] in genes:
					CDS[0] = model_to_id[CDS[0]]
			CDS_val = CDS_change[CDS[6]]			
			if fragment:
				code = checkStartStop(sequence,[CDS[3],CDS[4]],gene_id=gene_id)
			else:
				code = 10
			if code == 3 or code == 4:
				CDS_val = -1*CDS_val
			rec_info = {"start":CDS[3],"end":CDS[4],"type":"CDS","strand":CDS_val,"id":mRNA_id+".cds"+str(index)}
			sub_qualifiers2 = sub_qualifiers_base.copy()
			sub_qualifiers2["ID"] = [mRNA_id+".cds"+str(index)]
			sub_qualifiers2["Parent"] = [mRNA_id]
			sub_features.append(new_feature(rec_info,sub_qualifiers2))
			'''create UTRs'''
			
			'''First find out if a UTR is supposed to be added (for fragments)
				if there is a start or a stop codon at CDS start/end then add UTR
			'''
			if mRNA[3] != CDS[3] or mRNA[4] != CDS[4]:
				'''create 3' UTR'''
				if mRNA[4] > CDS[4] and (code == 2 or code == 4 or code == 10):
					sub_qualifiers3 = sub_qualifiers_base.copy()
					sub_qualifiers3["ID"] = [mRNA_id+".utr3p"+str(index)]
					sub_qualifiers3["Parent"] = [mRNA_id]
					UTR_start = CDS[4]+1
					UTR_end = mRNA[4]
					rec_info = {"start":UTR_start,"end":UTR_end,"type":"three_prime_UTR","strand":CDS_change[CDS[6]],"id":gene_id+".utr3p"+str(index)}
					sub_features.append(new_feature(rec_info,sub_qualifiers3))
				'''create 5' UTR'''
				if mRNA[3] < CDS[3] and (code == 1 or code == 3 or code == 10):
					sub_qualifiers4 = sub_qualifiers_base.copy()
					sub_qualifiers4["ID"] = [mRNA_id+".utr5p"+str(index)]
					sub_qualifiers4["Parent"] = [mRNA_id]
					UTR_start = mRNA[3]
					UTR_end = CDS[3]-1
					rec_info = {"start":UTR_start,"end":UTR_end,"type":"five_prime_UTR","strand":CDS_change[CDS[6]],"id":gene_id+".utr5p"+str(index)}
					sub_features.insert(0,(new_feature(rec_info,sub_qualifiers4)))
			index += 1
	top_feature.sub_features = sub_features
	if create_gene:
		gene_feature.sub_features = [top_feature]
		return gene_feature
	return top_feature

def change_gene_name(feature,geneName=False):
	if geneName:
		feature.qualifiers["Name"] = [geneName]
		feature.qualifiers["ID"] = [geneName]
		feature.id = geneName
	else:
		feature.qualifiers["Name"] = [""]
	return feature

def change_mRNA_name(feature,geneName,transcript_i):
	if transcript_i and geneName:
		transcript_id = geneName+"."+str(transcript_i)
		#print transcript_id
	else:
		transcript_id = geneName
	if geneName:
		feature.qualifiers["Name"] = [transcript_id]
		feature.qualifiers["ID"] = [transcript_id]
		feature.qualifiers["Parent"] = [geneName]
	else:
		feature.qualifiers["Name"] = [""]
	#print "mRNA",feature
	return feature

def change_exon_name(feature,geneName,exon_i):
	if exon_i:
		exon_id = geneName+".exon"+str(exon_i)
	else:
		exon_id = geneName+".exon1"
	feature.id = exon_id
	feature.qualifiers["ID"] = [exon_id]
	feature.qualifiers["Parent"] = [geneName]
	#print "exon",feature
	return feature

def change_CDS_name(feature,geneName,cds_i):
	cds_id = geneName+".cds"+str(cds_i)
	feature.id = cds_id
	feature.qualifiers["ID"] = [cds_id]
	feature.qualifiers["Parent"] = [geneName]
	#print "CDS",feature
	return feature

def change_UTR_name(feature,geneName,utr_i,utr_type):
	if utr_type.startswith("three"):
		utr_id = geneName+".utr3p"+str(utr_i)
	else:
		utr_id = geneName+".utr5p"+str(utr_i)
	feature.id = utr_id
	feature.qualifiers["ID"] = [utr_id]
	feature.qualifiers["Parent"] = [geneName]
	#print "UTR",feature
	return feature

def update_annotation(feature,geneName=False):
	'''Update the names ID attributes to proper gene names, otherwise remove NO Annotation'''
	
	### First update name of the gene 
	feature = change_gene_name(feature,geneName)
	transcript_i = len(feature.sub_features)-1
	top_f = []
	## Update name of all subfeatures in gene
	for subf in feature.sub_features:
		updated = []
		_type = subf.type
		subf.id = geneName
		if _type == "mRNA":
			exon_i = 0
			cds_i = 0
			transcript_i += 1
			up = change_mRNA_name(subf,geneName,transcript_i)
			#if geneName == "Potri.001N012400":
			#	print top_f
			#updated.append(change_mRNA_name(subf,geneName,transcript_i))
			geneName+="."+str(transcript_i)
			for subf2 in up.sub_features:
				sub_type=subf2.type
				if sub_type == "exon":
					exon_i += 1
					updated.append(change_exon_name(subf2,geneName,exon_i))
				if sub_type == "CDS":
					cds_i += 1
					updated.append(change_CDS_name(subf2,geneName,cds_i))
				if sub_type.endswith("UTR"):
					updated.append(change_UTR_name(subf2,geneName,transcript_i,_type))
			up.sub_features = updated
			top_f.append(up)
		else:
			'''There might be a few CDS etc annotated without mRNA they have to be filtered'''
			return feature
		#updated.append(subf)
	#feature.sub_features = updated
	try:
		feature.sub_features = top_f
	except UnboundLocalError:
		'''There might be a few CDS etc annotated without mRNA they have to be filtered'''
		return feature
	return feature

''' Write features when parsed '''

def new_feature(rec_info,sub_qualifiers):
	'''Create new feature'''
	feature = 	SeqFeature(
					FeatureLocation(
						rec_info["start"],
						rec_info["end"]
					),
					id=rec_info["id"],
					type=rec_info["type"],
					strand=rec_info["strand"],
					qualifiers=sub_qualifiers
				)
	return feature

def temp_to_annotation():

	pass

def update_feature_location(feature,start):
	return feature
	'''Fuction takes a feature and update its feature location relative to a supplied start value'''
	s = start + feature.location.start-1
	e = start + feature.location.end
	if s > e and feature.location.strand != -1:
		feature.location.strand *= -1
		s,e = e,s ## Swap the two values, start should always be lower than end
	#print "############################################"
	#print feature
	feature.location = FeatureLocation(s,e,feature.location.strand)
	#print "-------"
	#print feature
	#print "############################################"
	return feature

def create_sub_features(features,sub_qualifiers):
	'''Create a list of sub features'''
	sub_features = []
	for feat in features:
		new_feature(feat,sub_qualifiers)

def new_gene(input,qualifiers,sub_qualifiers):
	rec = SeqRecord(sequence,"Ptricocarpa_lincRNA_update_2015-07-31")
	top_feature = new_feature(rec_info,sub_qualifiers)
	top_feature.sub_features = create_sub_features(sub_f)
	return top_feature

def read_gene_annotation(in_file):
	'''return dictionary that translates either temp gene or temp model to its correct Potri annotation'''
	dic = {}
	with open(in_file,"r") as f:
		for row in f:
			data = row.strip("\n").split("\t")
			dic[data[0]] = "Potri"+data[2].strip("Potri")
			dic[data[1]] = "Potri"+data[2].strip("Potri")
	return dic

def parse_id(gene_id,features,gene,model_to_id=False,create_gene=False,filt_genes=False):
	if model_to_id:
		if args.gene_filter:
			if gene_id in model_to_id.keys():
				gene_id = model_to_id[gene_id]
				if gene_id in filt_genes: # and gene_id not in ovlp:
					gene = update_annotation(gene,gene_id)
				else: return False
			else:
				return False
		else:
			if gene_id in model_to_id.keys(): # and gene_id not in ovlp:
				gene_id = model_to_id[gene_id]
			gene = update_annotation(gene,gene_id)
	'''Get frame DP updated sub features, if they exists'''
	try:
		#sub_features = features[gene_id]
		if create_gene:	
			new_gene = features[gene_id]
			mRNA = new_gene.sub_features[0]
			#new_gene = update_feature_location(new_gene,start=gene.location.start)
		else:
			mRNA = features[gene_id]
		#mRNA = update_feature_location(mRNA,start=gene.location.start)
		mRNA_features = mRNA.sub_features
		print mRNA_features
		#for sub_feature in range(len(mRNA_features)):
		#	'''For every feature, update its gene start and end'''
		#	mRNA_features[sub_feature] = update_feature_location(mRNA_features[sub_feature],start=gene.location.start)
		if create_gene and args.replace_genes:
			return [new_gene]
		elif create_gene:
			return [gene,new_gene]
		gene.sub_features.append(mRNA)
	except KeyError:
		'''If no such update exists only print back to GFF3'''
		if create_gene:
			return [gene]
		return gene
	return gene

def print_tmp(feature):
	'''pring gff tmp file'''
	with open("tmp.fa", "w") as f:
		protein_fasta = ">%s\n%s" % (feature["FASTA"]["id"],feature["prot"]["sequence"])
		#protein_fasta = ">%s\n%s" % (feature["FASTA"]["id"],feature["FASTA"]["sequence"])
		print >> f,protein_fasta
	return "tmp.fa"

def exonerate_correct(exonerate_info,feature):
	'''Function to update and correct the alignment to fit the genome, as well as add back the UTR info from the mRNA'''
	#print feature["sequence-region"]["mRNA"]
	#print feature["sequence-region"]["CDS"]
	fpUTR = feature["sequence-region"]["mRNA"][0][3] - feature["sequence-region"]["CDS"][0][3]
	tpUTR = feature["sequence-region"]["mRNA"][0][4] - feature["sequence-region"]["CDS"][-1][4]
	## update exonerate info to include UTRs
	#print int(exonerate_info["gene"][0]["start"])-fpUTR
	#print int(exonerate_info["gene"][0]["end"])+tpUTR
	#print exonerate_info["gene"][0]["strand"]
	exonerate_info["mRNA"] = [{"start":int(exonerate_info["gene"][0]["start"])-fpUTR,"end":int(exonerate_info["gene"][0]["end"])+tpUTR,"strand":exonerate_info["gene"][0]["strand"]}]
	exonerate_info["gene"][0]["start"] = int(exonerate_info["gene"][0]["start"])-fpUTR
	exonerate_info["gene"][0]["end"] = int(exonerate_info["gene"][0]["end"])+tpUTR
	## update feature info with corrected positions.
	'''update gene'''

	'''update mRNA'''
	#print feature["sequence-region"]["mRNA"]
	feature["sequence-region"]["mRNA"][0][3] = exonerate_info["mRNA"][0]["start"]
	feature["sequence-region"]["mRNA"][0][4] = exonerate_info["mRNA"][0]["end"]
	feature["sequence-region"]["mRNA"][0][6] = exonerate_info["mRNA"][0]["strand"]
	print feature["sequence-region"]["CDS"]
	print range(len(feature["sequence-region"]["CDS"]))
	print exonerate_info
	for i in range(len(feature["sequence-region"]["CDS"])):
		'''update CDS'''
		try:
			feature["sequence-region"]["CDS"][i][3] = int(exonerate_info["cds"][i]["start"])
			feature["sequence-region"]["CDS"][i][4] = int(exonerate_info["cds"][i]["end"])
			feature["sequence-region"]["CDS"][i][6] = exonerate_info["cds"][i]["strand"]
		except IndexError:
			print i
			pass
	return exonerate_info,feature


def run(model_to_id=[],create_gene=False):
	sub_qualifiers = {"source":"RNAseq"}
	features = {}
	#frameDP_gff = ReadFrameDP("/mnt/picea/home/david/wood/PASA/20150211/frameDP/gff_files/Wood-frameDP-tab-sep.gff3")
	frameDP_gff = ReadFrameDP(args.frameDP)
	i = 0
	#x = False
	count = 0
	source = "/mnt/picea/storage/reference/Populus-trichocarpa/v3.0/fasta/Ptrichocarpa_v3.0_210.fa"
	for feature in frameDP_gff:
		if model_to_id:
			if feature["FASTA"]["id"] in model_to_id.keys():
		#		x = True
				feature["FASTA"]["id"] = model_to_id[feature["prot"]["id"]]
				#print feature["FASTA"]
				#'''Re reverse complement sequence for correct annotations, also update where to find the CDS'''
				#if "reverse" in feature["FASTA"]["frameDP"]:
				#	start = len(sequence) - CDS[1]
				 #	stop = start + (CDS[1]-CDS[0])
				#### instead of my create sub feature, use exonerate.
				target = print_tmp(feature)
				cmd = '''exonerate --forcegtag yes -m p2g -n 1 --showtargetgff yes --showvulgar no --verbose no --showalignment no %s %s''' % (target,source)
				#print cmd
				cmd = cmd.split(" ")
				res = subprocess.check_output(cmd,subprocess.STDOUT)
				print feature["sequence-region"]["mRNA"]
				#print res
				algn_info = {}
				
				for row in res.strip("\n").split("\n"):
					if not row.startswith("#"):
						row = row.split("\t")
						print row
						if len(row) < 3:
							continue
						if not row[2].startswith("sim"):
							try:
								algn_info[row[2]].append({"start":row[3],"end":row[4],"strand":row[6]})
							except:
								algn_info[row[2]] = [{"start":row[3],"end":row[4],"strand":row[6]}]
				#print feature["sequence-region"]["mRNA"]
				#feature["sequence-region"]["mRNA"][0][6] = algn_info["gene"]["strand"]
				#feature["sequence-region"]["mRNA"][0][3] = int(algn_info["gene"]["start"])
				#feature["sequence-region"]["mRNA"][0][4] = int(algn_info["gene"]["end"])+(int(feature["sequence-region"]["mRNA"][0][4])-int(feature["sequence-region"]["CDS"][-1][4]))
				#for cds_info in range(len(algn_info["cds"])):
				#	feature["sequence-region"]["mRNA"][cds_info][6] = algn_info["cds"][cds_info]["strand"]
				#	feature["sequence-region"]["mRNA"][cds_info][3] = int(algn_info["cds"][cds_info]["start"])
				#	feature["sequence-region"]["mRNA"][cds_info][4] = int(algn_info["cds"][cds_info]["end"])+(int(feature["sequence-region"]["CDS"][cds_info][4])-int(feature["sequence-region"]["CDS"][cds_info][4]))
				if "gene" in exonerate_info[0].keys():
					exonerate_info,feature = exonerate_correct(algn_info,feature)
				else:
					print "-----------------Warning-----------------------"
					print "Warning, exonerate info not correct",
					print res
					print "-----------------Warning-----------------------"
				#print exonerate_info
				#print feature
				#print feature
				#print algn_info
				features[feature["FASTA"]["id"]] = create_sub_feature(feature["sequence-region"]["mRNA"],feature["sequence-region"]["CDS"],model_to_id,create_gene=create_gene,sequence=feature["FASTA"]["sequence"])
				#features[feature["FASTA"]["id"]] = create_sub_feature(feature["sequence-region"]["mRNA"],feature["sequence-region"]["CDS"],model_to_id)
				print features[feature["FASTA"]["id"]]
				print "\n=========================\n================================"
		#if x:
		
		#count+=1
		#print count
		#if count == 5:
		#	break
		#	exit()
	return features

if __name__=="__main__":
	if args.gene_annotation:
		print "add gene annotations"
		model_to_id = read_gene_annotation(args.gene_annotation)
	else:
		model_to_id = False
	#print model_to_id
	#exit()
	features = run(model_to_id,args.create_genes)
	'''Get gene filter if required'''
	if args.gene_filter:
		filt_genes = [ x.strip("\n") for x in open(args.gene_filter)]
	else:
		filt_genes = False
	print len(features), "frameDP updates"
	in_file = args.input
	in_handle = open(in_file)
	rec_list = []
	for rec in GFF.parse(in_handle):
		#print rec
		new_features = []
		for gene in rec.features:
			if gene.id == "temp_gene_198":
				pass#print "orig",gene
			gene = parse_id(gene.id,features,gene,model_to_id,args.create_genes,filt_genes)	
			if gene:
				for g in gene:
					if g.id == "Potri.007N006900":
						print g
					if g.id == "temp_gene_198":
						print g
						exit()
			if args.create_genes and gene:
				for g in gene:
					new_features.append(g)
			elif gene:
				new_features.append(gene)

			#print gene
		#	print gene
		#GFF.write([gene],out_handle)		
		#print len(new_features)
		rec.features=new_features
		if len(new_features)>0:
			rec_list.append(rec)
	if args.output:
		if os.path.exists(args.output):
			res = raw_input("Output "+args.output+" file already exists, overwrite? (y/n): ")
			if res != "y" or res == "Y":
				exit("Script aborted")
		out_handle = open(args.output,"w")
	else:
		print "No output file was supplied writing to standard output file output.txt"
		out_handle = open("output.txt","w")
	GFF.write(rec_list,out_handle)
	in_handle.close()


'''
in_file = sys.argv[1]
in_handle = open(in_file)
for rec in GFF.parse(in_handle):
	print rec
	print "features in chromosome"
	for gene in rec.features:
		parse_id(gene.id)

in_handle.close()
'''