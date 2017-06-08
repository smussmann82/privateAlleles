from __future__ import print_function

from popmap import Popmap

import re
import collections
from collections import defaultdict

class VCF():
	'Class for parsing a VCF file'

	def __init__(self, infile, populations):
		#member variables
		self.header = list()
		self.preColumns = defaultdict(list)
		self.dataColumns = defaultdict(list)
		self.blacklist = set()
		self.privateAlleles = defaultdict(defaultdict(list).copy)
		self.loci = 0

		data = open(infile, 'r')
		content = data.readlines()
		self.data = [x.rstrip('\n') for x in content]
		self.removeHeader()
		self.transform()
		self.separatePops(populations)
	
	def removeHeader(self):
		counter=0
		removeList = list()
		for item in self.data:
			if(item.startswith('##')):
				self.header.append(item)
				removeList.append(counter)
			counter+=1

		for i in sorted(removeList, reverse=True):
			del self.data[i]

	def transform(self):
		dataHeader = self.data.pop(0).split()
		counter=0
		for line in self.data:
			counter+=1 #count the number of loci
			temp = line.split()
			for i,key in enumerate(dataHeader):
				if(i<9):
					self.preColumns[key].append(temp[i])
				else:
					self.dataColumns[key].append(temp[i])
		self.loci = counter

	def separatePops(self,populations):
		tempblacklist = list()
		regex = re.compile(r'./.') #compile regex to check for missing data
		for i in xrange(0, self.loci):
			tempdict = defaultdict(dict)
			missdict = defaultdict(int) #record missing data per population per locus
			allelePopCount = defaultdict(list)
			for key,value in self.dataColumns.iteritems():
				#check to see if locus is missing in either population
				if(regex.match(self.dataColumns[key][i])):
					missdict[populations.popmap[key]]+=1 #tally missing loci per population

				#if not missing, track which alleles appear in each population
				else:
					temp_list = self.dataColumns[key][i].split("|")
					for allele in temp_list:
						if allele not in allelePopCount[populations.popmap[key]]:
							allelePopCount[populations.popmap[key]].append(allele)
			#print(allelePopCount)

			#check if locus is absent for any population.  If so, blacklist it.
			for key,value in missdict.iteritems():
				if(missdict[key] == populations.popnums[key]):
					tempblacklist.append(i)

			#check to see if all alleles appear in all populations.  Record any alleles that are private, and record population in which it occurs.
			for key,value in allelePopCount.items():
				temp = value;
				#print(value)
				for allele in temp:
					#print(allele)
					counter=0 #count number of populations in which allele does not occur
					for pop,nums in populations.popnums.iteritems():
						if( allele not in allelePopCount[pop]):
							counter+=1
							#print(allele, "not in", pop)
					#check if allele is only in this population
					if(counter == populations.totalpops-1):
						#print("Private allele!")
					#add to list of private alleles
						self.privateAlleles[int(self.preColumns["#CHROM"][i])][key].append(allele)
		self.blacklist = set(tempblacklist)

	def printFile(self,outfile):
		#can probably easily generalize this to make a "printHeader" function
		f = open(outfile, 'w')
		for line in self.header:
			f.write(line)
			f.write("\n")
		f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
		od = collections.OrderedDict(sorted(self.dataColumns.items()))
		for key,value in od.iteritems():
			f.write("\t")
			f.write(key)
		f.write("\n")

		#actual printing of data begins here - can probably generalize into a "printData" function
		for i in xrange(0, self.loci):
			if( i not in self.blacklist):
				f.write(self.preColumns["#CHROM"][i])
				f.write("\t")
				f.write(self.preColumns["POS"][i])
				f.write("\t")
				f.write(self.preColumns["ID"][i])
				f.write("\t")
				f.write(self.preColumns["REF"][i])
				f.write("\t")
				f.write(self.preColumns["ALT"][i])
				f.write("\t")
				f.write(self.preColumns["QUAL"][i])
				f.write("\t")
				f.write(self.preColumns["FILTER"][i])
				f.write("\t")
				f.write(self.preColumns["INFO"][i])
				f.write("\t")
				f.write(self.preColumns["FORMAT"][i])
				for key,value in od.iteritems():
					f.write("\t")
					f.write(od[key][i])
				f.write("\n")

		f.close()


	def printPrivate(self,outfile):
		out = outfile + ".private"

		f = open(out, 'w')

		for key,value in sorted(self.privateAlleles.items()):
			for pop,allele in value.items():
				for a in allele:
					f.write(str(key))
					f.write("\t")
					f.write(pop)
					f.write("\t")
					f.write(a)
					f.write("\n")
			#print(value)
		f.close()

