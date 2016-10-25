#!/usr/bin/python2.7

'''Class Attibutes'''

import warnings

class Attributes(object):
	"""docstring for Attributes"""
	def __init__(self, _id=False, name=False,parent=False,alias=False,target=False,longest=False,pacid=False):
		super(Attributes, self).__init__()
		self.id = _id
		self.name = name
		self.alias = alias
		self.target = target
		self.parent = parent
		self.longest = longest
		self.pacid = pacid

	'''Get function to get all variables'''

	def get_id(self):
		return self.id

	def get_name(self):
		return ",".join(self.name)

	def get_alias(self):
		return ",".join(self.alias)

	def get_target(self):
		return ",".join(self.target)

	def get_parent(self):
		return ",".join(self.parent)

	def get_longest(self):
		return self.longest

	def get_pacid(self):
		return ",".join(self.pacid)

	def get_info(self):
		'''return all information in the class as a dictionary'''
		return {
			"id":self.id,
			"name":self.name,
			"target":self.target,
			"parent":self.parent,
			"alias":self.alias,
			"longest":self.longest,
			"pacid":self.pacid
		}

	'''Set function in case variables needs to be changed after initiation'''

	def set_id(self,id):
		if id:
			self.id = id
			return self.get_id()
		else:
			self.id = False

	def set_name(self,name):
		if name:
			self.name = name.split(',')
		else:
			self.name=False

	def set_alias(self,alias):
		if alias:
			self.alias = alias.split(',')
		else:
			self.alias=False

	def set_target(self,target):
		if target:
			self.target = target.split(',')
		else:
			self.target=False

	def set_parent(self,parent):
		if parent:
			self.parent = parent.split(',')
		else:
			self.parent = False

	def set_longest(self,longest):
		if longest:
			self.longest = longest
		else:
			self.longest = False

	def set_pacid(self,pacid):
		if pacid:
			self.pacid = pacid.split(',')
		else:
			self.pacid = False

	'''Function to append a second name alias etc to an attribute'''

	def add_name(self,name):
		self.name.append(name)

	def add_alias(self,alias):
		self.alias.append(alias)

	def add_target(self,target):
		self.target.append(target)

	def add_parent(self,parent):
		self.parent.append(parent)

	'''Read gff file attributes'''

	def read_attributes(self,string):
		attributes = string.split(";")
		'''Some gff files have whitespace around their ; they will be removed
			function return the id attribute
		'''
		attributes = [x.strip(" ") for x in attributes]
		for attribute in attributes:
			id,val = attribute.split("=")
			id = id.lower()
			if id == "id":
				self.set_id(val)
			elif id == "name":
				self.set_name(val)
			elif id == "parent":
				self.set_parent(val)
			elif id == "alias":
				self.set_alias(val)
			elif id == "target":
				self.set_target(val)
			elif id == "longest":
				self.set_longest(val)
			elif id == "pacid":
				self.set_pacid(val)
			else:
				#print string
				warnings.warn("The "+id+" attribute is not stored by the attributes object",Warning)
		return self.get_id()

	'''Internal functions for Attributes'''

	def _build_print(self):
		''' Support funciton to build the print for a correct attribute output '''
		string = "ID=%s" % (self.get_id())
		if self.name:
			string+=";Name=%s" % (self.get_name())
		if self.longest:
			string+=";longest=%s" % (self.get_longest())
		if self.parent:
			string+=";Parent=%s" % (self.get_parent())
		if self.pacid:
			string+=";pacid=%s" % (self.get_pacid())
		if self.alias:
			string+=";Alias=%s" % (self.get_alias())
		if self.target:
			string=";Target=%s" % (self.get_target())
		return string

	def __repr__(self):
		return Attributes()

	def __str__(self):
		return self._build_print()

if __name__=="__main__":
	'''test object'''
	print "Testing object"
	ids = set()
	attribute1 = Attributes("gene00001",name=["testgene"])
	print attribute1
	attribute2 = Attributes("tfbs000001",parent=["gene00001"])
	print attribute2
	attribute3 = Attributes("mRNA00001",parent=["gene00001"])
	attribute3.set_name("Eden1")
	attribute3.add_name("testmRNA")
	print attribute3
	attribute4 = Attributes()
	ids.add(attribute4.read_attributes("ID=cds00002 ; Parent=mRNA00002 ; Name=edenprotein.2"))
	print attribute4
	attribute5 = Attributes()
	ids.add(attribute5.read_attributes("ID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003;Invalid=test"))
	print attribute5
	attribute6 = Attributes()
	ids.add(attribute6.read_attributes("ID=PAC:27043735;Name=Potri.001G000100.1;pacid=27043735;longest=1;Parent=Potri.001G000100"))
	print attribute6
	print ids
