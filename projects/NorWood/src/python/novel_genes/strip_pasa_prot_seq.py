#!/usr/bin/python

'''
	Strip certain features from a PASA gff3, but keep the pasa gff3 file structure.
	Author: David Sundell (c) 2015
	Modified: 2015-03-09
'''

####
usage= """strip_pasa_prot_seq.py pasa.gff3 output_file
		all filter options are optional, but if non is set, there will be no filtering.
		optional: 	status_filter:int 
				  	valid_filter:True
				  	filter_id:Potri,PAC
		"""

''' 
	Import and system arguments
'''

import sys

if len(sys.argv) == 1 or sys.argv[1] == "-h":
	print usage
	sys.exit(0)
sys.argv.pop(0)
infile = sys.argv.pop(0)
outfile = open(sys.argv.pop(0),"w")
pasa_out=""

options = {}
for arg in sys.argv:
	arg = arg.split(":")
	options[arg[0]] = arg[1]

print options

'''
	Functions
'''

def parse_pasa_update_string(string):
	# PASA_UPDATE: PAC:27014049, single gene model update, valid-1, status:[pasa:asmbl_168722,status:12], valid-1
	## also make it work for
	# PASA_UPDATE: novel_model_89_54f5ee74, new gene addition, valid-1, status:[pasa:asmbl_168362,status:10], valid-1
	# PASA_UPDATE: temp_model_93274.1.54f72083, alt-splice addition, valid-0, status:[pasa:asmbl_168363,status:28], valid-0

	_row = string.strip().split(" ")
	_id = _row[2].strip(",")
	if _row[-1] == "valid-1":
		_valid = True
	else:
		_valid = False
	_asmbl,_status = _row[-2].split("pasa:")[1][:-2].split(",")
	_status = _status.split(":")[1]
	#print _id,_valid,_asmbl,_status
	return {"id":_id,"valid":_valid,"asmbl":_asmbl,"status":_status,"string":string}

def parse_prot_string(string):
	_row = string.split("\t")
	_prot = _row.pop()
	_row = _row[0].split(" ")
	## Remove ##PROT
	_name = _row.pop(0)

	_names = ""
	for _name in _row:
		_name = _row.pop(0)
		_names+=_name+" "
	_names.strip(" ")
	#print _names,_prot
	return [_names,_prot]
	##PROT PAC:27014049 Potri.T103200        MDHNSFEGALPSEIGNMKNLEI*

def parse_gff_string(string):
	''' Function that parses the gff3 string and fetches the ID and Name features, 
		returns the first part of the string so that ID and Name manipulations are possible
	'''
	_row = string.strip("#").split("ID=")
	_names = _row[-1].strip("\n").split(";Name=")
	if len(_names) < 2:
		_names = _row[-1].strip("\n").split(";Parent=")
	return {"ID":_names[0],"Name":_names[1],"comp_str":_names[0]+" "+_names[1],"row_pt1":_row[0]}
	#Chr01   .       gene    27088030        27097991        .       -       .       ID=Potri.001G261600;Name=%2A%2A%20NO%20NAME%20ASSIGNED%20%2A%2A

def filter_ids(string, options):
	''' Filter for id, or remove id

	'''
	pass

'''
	Main script
'''
testSet = set()

filter_status = False
filter_ids = False
# Filter the file according to the different update status ids
if "status_filter" in options.keys(): filter_status = True

# Filter ids gives an option to filter the file by ids found by pasa update or that is in the ID field of the gff3 file, 
# several ID keys are possible to add separated by ,
if "filter_id" in options.keys(): 
	filter_ids = True
	options["filter_id"] = options["filter_id"].split(",")
	first_row = True

valid = False
if "valid_filter" in options.keys(): 
	if "valid_filter" == "True" : valid = True

write = False
reached_prot = False
next = False
found_status=False

if not valid and not filter_status and not filter_ids:
	print usage
	sys.exit("No filter option was set, please add a filter.")


'''
	Walk through input file line by line and write complete pasa output unless filtered.
'''

with open(infile, "r") as f:
	for row in f:
		if row.startswith("# PASA_UPDATE"):
			if next and not reached_prot:
				continue
			else:
				next = False
			
			res = parse_pasa_update_string(row)
			
			## If a set of updates are read write to file
			## write previous update and reset parameters
			#print write
			if write:
				#print "line printed"
				outfile.write(pasa_out)
				#print pasa_out
				pasa_out=""
				write = False
				found_status=False
				if filter_ids:
					first_row = True

			if valid:
				if res["valid"]:
					next = True
					continue
			if filter_ids:
				for test in options["filter_id"]:
					if test in res["id"]:
						next = True
						write = False
						continue
			if filter_status:
				if options["status_filter"] == res["status"]:
					pasa_out+=row
					found_status = True
				elif found_status:
					pasa_out+=row
				else:
					pasa_out = ""
					next = True
					write = False
					continue
			else:
				pasa_out+=row
		elif row.startswith("#PROT"):
			reached_prot=True
			if next:
				pasa_out=""
				write = False
				continue
			parse_prot_string(row)
			pasa_out+=row
			write = True
		else:
			if next:
				continue
			if filter_ids and first_row and row != "\n":
				res = parse_gff_string(row)
				for test in options["filter_id"]:
					if False:
						print "Is ",test, " in ",res["comp_str"]," -> ",
						#print row
					if test in res["comp_str"]:
						#print "True"
						#print "next set write = false"
						next = True
						write = False
						continue
				#first_row = False
			pasa_out+=row
		
	if write:
		#print pasa_out
		outfile.write(pasa_out)
	outfile.close()

print testSet
print "Done"
