from __future__ import print_function

from popmap import Popmap

import re
from collections import defaultdict


class VCF():
	'Class for parsing a VCF file'

	def __init__(self, infile, populations):
		#member variables
		self.header = list()
		self.preColumns = defaultdict(list)
		self.dataColumns = defaultdict(list)
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
		regex = re.compile(r'./.') #compile regex to check for missing data
		for i in xrange(0, self.loci):
			tempdict = defaultdict(dict)
			missdict = defaultdict(int) #record missing data per population per locus
			for key,value in self.dataColumns.iteritems():
				print(key)
				#check to see if locus is missing in either population
				if(regex.match(self.dataColumns[key][i])):
					print(self.dataColumns[key][i])
					missdict[populations.popmap[key]]+=1 #tally missing loci per population

				#if not missing, track which alleles appear in each population
				else:
					print(self.dataColumns[key][i])

			#check if locus is absent for any population.  If so, blacklist it.
			print(missdict)

			#check to see if all alleles appear in all populations.  Record any alleles that are private, and record population in which it occurs.
