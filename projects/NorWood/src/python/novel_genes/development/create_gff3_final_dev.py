import sys
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from modules.read_fasta import ReadFrameDP
from argparse import ArgumentParser
import os

parser = ArgumentParser(description = """
	python merge_fa_and_frameDP_gff.py <fasta file> <merge file from frameDP> <out>
			 or any file with features starting with ">"

	run python merge_fa_and_frameDP_gff.py -h to get detailed help
	""")
parser.add_argument('input', metavar="", type=str, help="")
parser.add_argument('-g','--gene_annotation', metavar="", default=False,type=str, help="Chose a file if indexes needs to be updated")
parser.add_argument('-gf','--gene_filter', metavar="",default=[], type=str,help="Filter the final set of genes from the gene annotation file")
parser.add_argument('-f', '--frameDP', type=str,help="The source frameDP merged in_file")
parser.add_argument('-o','--output', metavar="", default=False, type=str, help="Add one additonal output file, for a merged gff3 file to be written")
parser.add_argument('--create_genes',action='store_true',help="Add to create new gene features for frameDP updates instead of adding them as an mRNA to existing")

args = parser.parse_args()
'''Parse features from file '''

def set_qualifiers(_id,base_qualifiers,_type,rec_info=False,parent=False,mRNA_i=False,mRNA=False,CDS=False,index=False,start=False,end=False):
	'''functino to set qualifiers for a new feature'''
	CDS_change = {'+':1,'-':-1}
	gene_qualifiers = base_qualifiers.copy()
	if parent:
		gene_qualifiers["Parent"] = [_id]
		_id = _id+"."+str(mRNA_i)
	gene_qualifiers["ID"] = [_id]
	if _type == "gene":
		gene_qualifiers["Name"] = [_id]
	if mRNA:
		rec_info = {"start":mRNA[3],"end":mRNA[4],"type":_type,"strand":CDS_change[CDS[6]],"id":_id}
	elif start:
		UTR_ends = {"three_prime_UTR":"utr3p","five_prime_UTR":"utr5p"}
		rec_info = {"start":start,"end":end,"type":_type,"strand":CDS_change[CDS[6]],"id":_id+"."+UTR_ends[_type]+str(index)}
	else:
		rec_info = {"start":CDS[3],"end":CDS[4],"type":_type,"strand":CDS_change[CDS[6]],"id":_id+"."+_type+str(index)}
	gene_feature = new_feature(rec_info,gene_qualifiers)

	return gene_feature


def create_sub_feature(mRNAs,CDSs,model_to_id,genes=[],create_gene=False):
	'''Main function to create a new gene or mRNA feature from frameDP outputs'''
	#print genes
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
			gene_id = gene_id[:-2]+str(gene_id_end)
			gene_id_end += 20
			gene_feature=set_qualifiers(gene_id,sub_qualifiers_base,_type="gene",mRNA=mRNA,CDS=CDSs[0])
		mRNA_id = gene_id+"."+str(i+1)
		top_feature = set_qualifiers(mRNA_id,sub_qualifiers_base,_type="mRNA",mRNA=mRNA,CDS=CDSs[0])
		index = 1
		for x in range(len(CDSs)):
			CDS = CDSs[x]
			'''Function to create sub features from mRNA and CDS'''
			sub_features.append(set_qualifiers(mRNA_id,sub_qualifiers_base,_type="exon",index=index,CDS=CDS))
			if model_to_id:
				if CDS[0] in genes:
					CDS[0] = model_to_id[CDS[0]]			
			sub_features.append(set_qualifiers(mRNA_id,sub_qualifiers_base,_type="cds",index=index,CDS=CDS))
			'''create UTRs'''
			if mRNA[3] != CDS[3] or mRNA[4] != CDS[4]:
				'''create 3' UTR'''
				if mRNA[4] > CDS[4]:
					UTR_start = CDS[4]+1
					UTR_end = mRNA[4]
					sub_features.append(set_qualifiers(mRNA_id,sub_qualifiers_base,_type="three_prime_UTR",index=index,CDS=CDS,start=UTR_start,end=UTR_end))
				elif mRNA[3] < CDS[3]:
					UTR_start = mRNA[3]
					UTR_end = CDS[3]-1
					sub_features.insert(0,set_qualifiers(mRNA_id,sub_qualifiers_base,_type="five_prime_UTR",index=index,CDS=CDS,start=UTR_start,end=UTR_end))
			index += 1
	top_feature.sub_features = sub_features
	if create_gene:
		gene_feature.sub_features = [top_feature]
		return gene_feature
	return top_feature

'''Help Functions to create new genes '''
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
			dic[data[0]] = "Potri."+data[2].strip("Potri")
			dic[data[1]] = "Potri."+data[2].strip("Potri")
	return dic


'''Functions for updating name, could be merged to fewer functions, but it is also clear with them separated'''

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

def set_feature_name(feature,geneName=False,id_i=False,id_type=False,startwith=False):
	if not geneName:
		feature.qualifiers["ID"]=[""]
		return feature
	name=geneName
	if id_type=="gene":
		feature.qualifiers["Name"] = [name]
	if id_type!="gene":
		id_name=name+"."+id_type+str(id_i)
	feature.id=name
	feature.qualifiers["ID"] = [name]
	


'''Main functions'''



def update_annotation(feature,geneName=False):
	'''Update the names ID attributes to proper gene names, otherwise remove NO Annotation'''
	
	### First update name of the gene 
	feature = change_gene_name(feature,geneName)
	transcript_i = len(feature.sub_features)
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
					updated.append(change_UTR_name(subf2,geneName,transcript_i,sub_type))
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

def update_feature_location(feature,start):
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


def parse_id(gene_id,features,gene,model_to_id=False,create_gene=False):
	if model_to_id:
		if args.gene_filter:
			filt_genes = [ x.strip("\n") for x in open(args.gene_filter)]
			if gene_id in model_to_id.keys():
				gene_id = model_to_id[gene_id]
				if gene_id in filt_genes:
					gene = update_annotation(gene,gene_id)
				else: return False
			else:
				return False
		else:
			if gene_id in model_to_id.keys():
				gene_id = model_to_id[gene_id]
			gene = update_annotation(gene,gene_id)
		#print gene
		'''Get frame DP updated sub features, if they exists'''
	try:
		#sub_features = features[gene_id]
		if create_gene:	
			new_gene = features[gene_id]
			mRNA = new_gene.sub_features[0]
			new_gene = update_feature_location(new_gene,start=gene.location.start)
		else:
			mRNA = features[gene_id]
		mRNA = update_feature_location(mRNA,start=gene.location.start)
		mRNA_features = mRNA.sub_features
		for sub_feature in range(len(mRNA_features)):
			'''For every feature, update its gene start and end'''
			mRNA_features[sub_feature] = update_feature_location(mRNA_features[sub_feature],start=gene.location.start)
		if create_gene:
			return [gene,new_gene]
		gene.sub_features.append(mRNA)
	except KeyError:
		'''If no such update exists only print back to GFF3'''
		if create_gene:
			return [gene]
		return gene
	return gene

def run(model_to_id=[],create_gene=False):
	sub_qualifiers = {"source":"RNAseq"}
	features = {}
	#frameDP_gff = ReadFrameDP("/mnt/picea/home/david/wood/PASA/20150211/frameDP/gff_files/Wood-frameDP-tab-sep.gff3")
	frameDP_gff = ReadFrameDP(args.frameDP)
	i = 0
	#x = False
	for feature in frameDP_gff:
		if model_to_id:
			if feature["FASTA"]["id"] in model_to_id.keys():
		#		x = True
				feature["FASTA"]["id"] = model_to_id[feature["FASTA"]["id"]]
		features[feature["FASTA"]["id"]] = create_sub_feature(feature["sequence-region"]["mRNA"],feature["sequence-region"]["CDS"],model_to_id,create_gene=create_gene)
		#features[feature["FASTA"]["id"]] = create_sub_feature(feature["sequence-region"]["mRNA"],feature["sequence-region"]["CDS"],model_to_id)
		#print features[feature["FASTA"]["id"]]
		#if x:
		#	exit()
	return features

if __name__=="__main__":
	if args.gene_annotation:
		print "add gene annotations"
		model_to_id = read_gene_annotation(args.gene_annotation)
	else:
		model_to_id = False
	features = run(model_to_id,args.create_genes)
	print len(features), "frameDP updates"
	in_file = args.input
	in_handle = open(in_file)
	rec_list = []
	for rec in GFF.parse(in_handle):
		#print rec
		new_features = []
		for gene in rec.features:
			#print gene
			gene = parse_id(gene.id,features,gene,model_to_id,args.create_genes)
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