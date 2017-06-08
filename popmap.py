from __future__ import print_function

from collections import defaultdict

class Popmap():
	'Class for parsing a popmap'

	def __init__(self, infile):
		#member variables
		self.popmap = dict()
		self.popnums = defaultdict(int)
		

		data = open(infile, 'r')
		content = data.readlines()
		data.close()

		content = [x.rstrip('\n') for x in content]


		for line in content:
			temp = line.split()
			self.popmap[temp[0]] = temp[1] #make map of individual->population
			self.popnums[temp[1]]+=1 #count number of individuals per population
		#print(self.popmap)
		#print(self.popnums)

		self.totalpops = len(self.popnums)
