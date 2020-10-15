# Anthony Androulakis
'''
from analyzeiupac import Molecule
name = "6,6-dimethyl-4-methylenebicyclo[3.1.1]hept-2-en-1-ol"
data = Molecule(name).data
print(data)
'''

class Molecule:
	def __init__(self, name):
		self.name = name
		self.substitutentmultipliers = ["di","tri","tetra","penta","hexa","hepta","octa","nona","deca"]
		self.cyclomultipliers = ["bi","tri","tetra","penta","hexa","hepta","octa","nona","deca"] #the only difference is using bi instead of di
		self.parents = ["meth","eth","prop","but","pent","hex","hept","oct","non","dec","undec","dodec", "tridec", "tetradec", "pentadec", "hexadec", "heptadec", "octadec", "nonadec", "eicos"]
		self.findparent = [i[::-1] for i in self.parents[::-1]]
		self.unsaturation = ["en", "yn"] #"an" excluded since it's so easy to check for
		self.functionalgroups = ["oic", "oate", "al", "one", "ol", "amine"] #these are not all the functional groups but are most
		self.data = self.analyzename()
	#returns the number of carbons in the parent chain (input is name split by -), index of splitname where parent name is found
	def parentfinder(self, splitname):
		searching = [i[::-1] for i in splitname[::-1]]
		for index,i in enumerate(searching):
			tempvar = []
			for j in self.findparent:
				if j in i:
					tempvar.append((i.find(j), self.findparent.index(j)))
			if len(tempvar)>0:
				return len(self.findparent)-min(tempvar)[1], len(searching)-index-1, min(tempvar)[0]
		return None, None, None
	#for finding cyclic/polycyclic molecules (up to decacyclo). for example, decacyclo would return 10 as the first value
	def cyclofinder(self, splitname, parentnum, indexofparentinsplitname, indexofparentinfocus): 
		focus = splitname[indexofparentinsplitname] #other terms don't matter for finding if cyclic
		reversedfocus = focus[::-1]
		indexofparentinfocus = focus.find(self.parents[parentnum-1])
		cyclopossibilities = ["cyclo"]+[i+"cyclo" for i in self.cyclomultipliers]
		for index, i in enumerate(cyclopossibilities[::-1]):
			findcyclo = reversedfocus.find(i[::-1])
			if findcyclo != -1:
				if i=="cyclo" and len(focus)-findcyclo == indexofparentinfocus:
					return len(cyclopossibilities)-index, len(focus)-(findcyclo+len(i))
				elif i!="cyclo" and focus.find(']')+1 == indexofparentinfocus:
					return len(cyclopossibilities)-index, len(focus)-(findcyclo+len(i))
		return None, None
	def functionalgroupfinder(self, splitname): #returns which functional group was used
		for index,i in enumerate(self.functionalgroups):
			if i == splitname[-1]:
				return index
		return None
	#central stuff?
	def analyzename(self):
		cutname = self.name.split(' ')[0]
		splitname = cutname.split("-")
		if splitname[0] == "(trans)" or splitname[0] == "(cis)":
			splitname = splitname[1:]
		parentnum, indexofparentinsplitname, focusindexfromend = self.parentfinder(splitname)
		indexofparentinfocus = len(splitname[indexofparentinsplitname])-focusindexfromend-len(self.parents[parentnum-1])
		#print(parentnum, indexofparentinsplitname, indexofparentinfocus)
		cyclotype, indexofcycloinfocus = self.cyclofinder(splitname, parentnum, indexofparentinsplitname, indexofparentinfocus)
		#print(cyclotype, indexofcycloinfocus)
		functionalgroup = self.functionalgroupfinder(splitname)
		theunsaturation = True #it doesn't really matter what you set theunsaturation to on this line, it just needs to be "initialized"
		if functionalgroup == None:
			splitname[-1] = splitname[-1][:-1] #chop off that useless e. we don't need it for this program.
		if focusindexfromend != 0:		
			unsaturation = None
			splitname[indexofparentinsplitname] = splitname[indexofparentinsplitname][:-2] #chop off that useless an.
		#time for the fun part
		thesubstituents = ('-'.join(splitname[:indexofparentinsplitname])+'-'+splitname[indexofparentinsplitname][:(indexofcycloinfocus if cyclotype != None else indexofparentinfocus)] if indexofparentinsplitname!=0 else None)
		thefunctionalgroup = (splitname[-2:] if functionalgroup != None else None)
		if thefunctionalgroup is not None:
			theunsaturation = (splitname[indexofparentinsplitname+1:-2] if theunsaturation != None else None)
		else:
			theunsaturation = (splitname[indexofparentinsplitname+1:] if theunsaturation != None else None)
		theparent = self.parents[parentnum-1]
		thecycloparent = (splitname[indexofparentinsplitname][indexofcycloinfocus:indexofparentinfocus] if cyclotype != None else None)
		print(f'parent: {theparent}\n cyclo: {thecycloparent}\n unsaturation: {theunsaturation}\n functionalgroup: {thefunctionalgroup}\n substituents: {thesubstituents}')
		return theparent, thecycloparent, theunsaturation, thefunctionalgroup, thesubstituents


